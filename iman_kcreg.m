%IMAN_KCREG
%   Kernal Correlation Registration, for 2D point clouds (as Nx2 coordinate
%   lists). Optimizes a rigid transformation (translation and rotation) to
%   align 'Model' coordinates to fixed 'Scene' coordinates. Coordinates
%   are approximated by Kernel Density Estimates (KDEs) via ksdensity,
%   wherein sampling points are distributed in a uniform grid, and each
%   point is affected by coordinates within a range, defined by BandWidths.
%   Optimizaiton is carried out in Levels, with the first pass using sparse
%   sampling with large BandWidths to roughly align the Model. Subsequent
%   Levels use denser sampling points and smaller BandWidths for more
%   accuracy. 
%
%   [mdl, tp] = iman_kcreg(scn, mdl, 'Name', 'Value', ...)
%
%   Inputs
%       scn: Nx2 array of Scene coordinates
%       mdl: Mx2 array of Model coordinates
%
%   Outputs
%       mdl: Mx2 array of transformed Model coordinates
%       tp:  transform parameters (x-translation, y-trans, rotation-angle)
%
%   Additional input options (as Name/Value pairs)
%       NumLevels   - [num] Number of Levels of optimization to perform
%           (default: 3). More levels take smaller steps in refining the
%           optimization problem (from broad KDEs to narrow).
%       ExtraLevels - [num] Number of extra Levels to perform (default: 0).
%           Extra levels increase the number of sampling points beyond
%           typical ranges.
%       bwChange    - [num] Change to make to BandWidths (default: 0).
%           Positive values increase BandWidths, negative values decrease
%           them. 
%       Display     -
%
%   The implementation of the Kernal Correlation Registration cost function
%   is based on Tsin and Kanade (2004), doi:10.1007/978-3-540-24672-5_44 
%
%   Version 1.0

%   M. Pargett, Albeck Lab, UC Davis, MCB.  2016

function [mdl, tp, ctf] = iman_kcreg(scn, mdl, varargin)
%Version check provision
if strcmpi(scn,'version'); mdl = 'v1.0'; return; end

%Operation parameter defaults
p.numlevels = 3;        %Default Number of Levels is 3
p.extralevels = 0;      %Default Extra Levels is zero (none)
p.bwchange = 0;       	%Default is no change to BandWidths
p.screen = 0;           %Default number of screen points is zero
p.display = false;    	%Default is to not display procedure

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end


%% Preliminary definitions
%Useful range for number of sampling points per dimension
nrng = [10, 50];
%   Define number of sample points for each Level
np = linspace(nrng(1), nrng(2), p.numlevels);
%   Append extra Levels at same spacing, if applicable
if p.numlevels == 1; nps = nrng(2); else nps = np(end) - np(end-1); end
np(end + (1:p.extralevels)) = np(end) + nps*(1:p.extralevels);

%Get number of points in each coordinate set
ns = size(scn,1);       nm = size(mdl,1);

%Set region bounds by range of coordinate values
mn0 = min(scn,[],1); mx0 = max(scn,[],1);
mn = mn0 - (mx0-mn0)*0.05;  mx = mx0 + (mx0-mn0)*0.05;  %Expand each by 5%
%   Store transform parameter bounds (for optimization problem)
tpb = [mx-mn, pi];

%Define coordinate transform function for rigid transformations
ctf = @(m,p)m*[cos(p(3)),-sin(p(3));sin(p(3)),cos(p(3));p(1),p(2)];
%   p = [x_translation, y_translation, rotation_angle]

%Set optimization options
oopt = optimoptions(@fmincon, 'Algorithm', 'active-set', ...
    'Display', 'off', 'UseParallel', true);


%% Prepare Model coordinates
%Preliminary alignment - as needed to ensure initial overlap
if any(min(mdl,[],1) > mx0) || any(max(mdl,[],1) < mn0);   ald = true;
    %Moves mean of Model to mean of Scene
    mnm = mean(scn,1) - mean(mdl,1);  mdl = bsxfun(@plus, mdl, mnm);
else ald = false;
end

%Augment coordinates for transformation
mdl = [mdl, ones(size(mdl,1),1)];


%% Initial Screening (if indicated)
%FIXME: Not yet implemented.  For now, ensure that images have notable
%features (large shapes) in cell mass, and that if trimming is necessary it
%is done very carefully for optimal alignment.

%% Iterate with refining KDE scale for robust solution
tp = zeros(p.numlevels+1,3);
for s = 1:p.numlevels
    %Set BandWidth (affects region influencing each sampling point)
    bw = tpb(1:2)/(max(s - p.bwchange,1)*(np(s)-1));   bw2 = prod(bw);
    %Build uniform grid of sampling points
    [xx,xy] = meshgrid( linspace(mn(1),mx(1),np(s)), ...
                        linspace(mn(2),mx(2),np(s)) );
    xi = reshape(cat(3,xx,xy), numel(xx), 2);  %List of grid points
    
    %Get Kernel Density Estimate (KDE) for Scene (at current scale)
    [kde_s, xi] = ksdensity(scn, xi, 'Kernel', @normpdf, ...
        'Bandwidth', bw, 'Function', 'pdf');
    kde_s = kde_s.*ns.*bw2;  %Scale toward unity (avoid tiny values)
    
    %Define KDE function for Model, with parameterized transform
    cost = @(p)-dot( kde_s, nm.*bw2.*ksdensity(ctf(mdl,p), xi, ...
        'Kernel', @normpdf, 'Bandwidth', bw, 'Function', 'pdf') );
    
    %Run optimization problem
    tp(s+1,:) = fmincon(cost, tp(s,:), [],[], [],[], -tpb,tpb, [], oopt);
    %   uses optimal parameters from last level
end


%% Display procedure
if p.display
    %Set plotting parameters (subset of trace indices, color)
    ss = 1:ceil(nm/100):nm;  c = 0.75;
    %Plot progression of points
    figure;   hold on;   
    for s = 1:p.numlevels+1; 	a = ctf(mdl,tp(s,:)); 
        if s > 1 %Plot line from last to current point
            arrayfun(@(x1,y1,x2,y2)plot([x1,x2],[y1,y2], ...
                '-', 'LineWidth', 1.5, 'Color', [1,c,c]), ...
                a0(ss,1), a0(ss,2), a(ss,1), a(ss,2));
        end;    a0 = a;
    end
    plot(scn(:,1),scn(:,2),'bo', 'MarkerSize', 6);    
    a = ctf(mdl,tp(1,:));    plot(a(:,1),a(:,2),'rx', 'MarkerSize', 7);
    a = ctf(mdl,tp(end,:));  plot(a(:,1),a(:,2),'ro', 'MarkerSize', 2);
end


%% Finalize
%Register Model with final optimized parameters
mdl = ctf(mdl, tp(end,:));
%Return mean value from scene if initial alignment performed
if ald;  tp(end,1:2) = tp(end,1:2) + mnm;  end
%Restrict transformation parameters to end values to return
tp = tp(end,:);

end

%Obsolete display code for debugging:
% figure;  ksdensity(scn, xi, 'Kernel', @normpdf, ...
%     'Bandwidth', bw, 'Function', 'pdf'); view(2);
% 
% figure; ksdensity(ctf(mdl,tp(end,:)), xi,'Kernel', @normpdf, ...
%     'Bandwidth', bw, 'Function', 'pdf'); view(2);
% 
% test = ctf(mdl,tp(end,:)); figure; plot(test(:,1),test(:,2),'rx');
% hold on; plot(scn(:,1),scn(:,2),'bo');
% 
% for s = 1:nlevels+1
%     test = ctf(mdl,tp(s,:)); figure; plot(test(:,1),test(:,2),'rx');
%     hold on; plot(scn(:,1),scn(:,2),'bo');
% end

