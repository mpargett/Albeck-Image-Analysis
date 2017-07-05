%IMAN_XYSHIFT
%   Correct for shift in XY position during imaging
%
%   [xyin, ctrans] = iman_xyshift(xyin, p, shd, fname, varargin)

%   ctrans - a structure contain the shift coordinate transformations
%       tidx - time indices (of first shifted frame)
%       tmat - transformation matrices
%
%   Version 2.0, 2016.08.08

function [xyin, ctrans] = iman_xyshift(xyin, p, shd, dao, varargin)
%Version check provision
if strcmpi(xyin,'version'); xyin = 'v2.0'; return; end
    
%Assess direction of shift desired (no scaling allowed)
if ~exist('shd', 'var'); shd = 1; else shd = sign(shd); end

%Assess for bad frames to be skipped
if isstruct(xyin) && isfield(xyin,'xCoord')
    badt = find(arrayfun(@(x)isempty(x.xCoord), xyin));
else badt = [];
end
        
%% Determine necessary shift inference procedure
if all(isfield(p,{'tidx','tmat'})) %IF ctrans provided as p
    ctrans = p; %Assign ctrans and skip to deshifting
elseif ~isempty(p.dx) && ~isempty(p.dy) %IF shift values are provided
    %Define the provided translation
    ctrans = sub_defined_shift(p);
    
else    %otherwise the shift requires inference
    %FIXME Parallelize here if many shifts to be performed.
    sz = 1;  %Generally not considering z-dimension
    
    %Determine tracking location (i.e. well / XY point)
    %   Use Tracking Well if provided, or Current Well, or default to 1
    if isfield(p,'trackwell') && ~isempty(p.trackwell); twell = p.trackwell;
    elseif  isfield(p,'curwell') && ~isempty(p.curwell); twell = p.curwell;
    else    twell = 1;
    end
    
    %Determine tracking channel
    if isfield(p,'trackchan') && ~isempty(p.trackchan)
        sc = p.trackchan;    else sc = [];      %IF not present, empty
    end
    
    %IF shift is specified to be only translation, use reduced procedure
    if isfield(p,'shifttype') && strcmpi(p.shifttype, 'spot')
        %THIS OPTION TO RUN SIMPLIFIED PROCEDURE (NOT YET IMPLEMENTED)
        infer_proc = @sub_image_registration;
    else %Use Translation as default (rotation/translation)
        infer_proc = @sub_image_registration;
    end
    
    %Include any manually indicated 'bad' time points
    if isfield(p, 'badt'); badt = unique([badt(:); p.badt(:)]); end
    
    %Perform shift inference for each shift event
    ctrans = struct('tidx', [], 'tmat', [], 'stime', p.stime); cst = 1;
    for st = 1:numel(p.frame)
        %Load pre-shift image
        tpre = p.frame(st)-1;   tpst = p.frame(st);  %Frame indices
        %   Check validity of pre and post frames (skip bad first frame)
        badi = find(tpre == badt);
        if ~isempty(badi); prange = p.stime:tpre;  %#ok<ALIGN>
            prange(badt(1:badi) - p.stime + 1) = []; 
        if isempty(prange);   continue; end;    tpre = prange(end);    end
        badi = find(tpst == badt);
        if ~isempty(badi); prange = tpst:(p.stime-1 + length(xyin));  %#ok<ALIGN>
            prange(badt(badi:end)-tpst+1) = []; 
        if isempty(prange);   continue; end;    tpst = prange(end);    end
        
        if isempty(sc); %IF no channel specified, check all channels
            nc = r.getSizeC; sc = 1:nc; im1 = zeros(r.getSizeY,r.getSizeX,nc);
            for s = 1:numel(sc); im1(:,:,s) = ...
                    iman_getframe(dao, [tpre, sc(s), twell, sz]); end
            [~,sc] = max(range( reshape(im1, numel(im1)./nc, nc), 1));
            im1 = double( im1(:,:,sc) );
        else
            im1 = iman_getframe(dao, [tpre, sc, twell, sz]);
        end
        
        %Load post-shift image        
        im2 = iman_getframe(dao, [tpst, sc, twell, sz]);
        
        %Infer transformation between images
        ctrans(cst).tidx = tpst;
        ctrans(cst).tmat = infer_proc(im1, im2, lower(p.shifttype));
        %IF multiple shifts exist, propagate previous to give total shift
        if cst > 1; ctrans(cst).tmat = ctrans(cst-1).tmat*ctrans(cst).tmat; end
        cst = cst + 1;
    end
end

%Generate de-shifted coordinates or image (as needed)
xyin = sub_deshift(xyin, ctrans, shd, badt, varargin{:});

end

%% Defined shift based on provided translations
function [ctrans] = sub_defined_shift(p)
ctrans = struct('tidx', [], 'tmat', [], 'stime', p.stime);
for s = 1:numel(p.frame)
    %Store time index of shift
    ctrans(s).tidx = p.frame(s);   
    %Store translation matrix to correct provided shift values
    ctrans(s).tmat = [[eye(2), [0; 0]]; [-p.dx(s), -p.dy(s), 1]];
    %IF multiple shifts exist, propagate previous to give total shift
    if s > 1; ctrans(s).tmat = ctrans(s-1).tmat*ctrans(s).tmat; end
end

end


%% Optimizing shift inference procedure
%   Runs MATLAB image registration procedures via imregtform
%   Adds a layer to screen over optimization parameters, and regularizes to
%       prevent 'fixing' shifts that are very small (likely erroneous)
function tmat = sub_image_registration(im1, im2, ttype)
%Disable warnings likely to arise
warning('off','images:regmex:registrationOutBoundsTermination');
%   Use default of translation if no transform specified
if ~exist('ttype','var') || isempty(ttype); ttype = 'translation'; end

cost_thresh = 0.1; %Set cost threshold for satisfactory alignment
stol = [1e-4, 1e-2, 3, 1e-2, 1e-4, 3];      %Shift tolerances
rw = 0.1;         %Set pseudo-regularization weight
nstep = 10;         %Set maximum number of step size parameters to test

%Self-normalize incoming images (constrains cost values)
im1 = im1./mean(im1(:));  im2 = im2./mean(im2(:));

%Get spatial reference object, to align images
ra = imref2d(size(im1));
%Get default optimizer and metric values for registration
[tform_opti, tform_met] = imregconfig('monomodal');
%   Modify gradient optimizer parameters (for any transformation type)
tform_opti.MaximumIterations = 1e3;     %Increased max cycles for accuracy
tform_opti.GradientMagnitudeTolerance = 1;  %Higher tolerance, for speed

%Set optimizer parameters by transformation type (ttype)
switch ttype
    case 'rigid'    %IF a rotation
        tform_opti.MinimumStepLength = 1e-5;    %Low, rotation is sensitive
        tform_opti.RelaxationFactor  = 0.65;    %Optimized per step sizes
        %Consider a range of max step sizes in log space. Nominal = 1e-2
        max_step = logspace(-3,-1,nstep);  ms_nomlog = -2;
    otherwise      %IF a translation (default)
        tform_opti.MinimumStepLength = 1e-3;    %Higher, no rotation here
        tform_opti.RelaxationFactor  = 0.8;     %Higher, more robust        
        %Consider a range of max step sizes in log space. Nominal = 1e-1
        max_step = logspace(-2,0,nstep); ms_nomlog = -1;
end
%Sort max step sizes to start search near nominal and progress away
[~,si] = sort( abs(log10(max_step) - ms_nomlog) ); max_step = max_step(si);

%Perform image registration, looping through the step size parameter
cost_val = inf(1,nstep); tmat = cell(1,nstep);
for s = 1:nstep     %FOR each step size parameter, start with default
    tform_opti.MaximumStepLength = max_step(s);     %Set step size
    tform = imregtform(im2, ra, im1, ra, ...
        ttype, tform_opti, tform_met, 'PyramidLevels', 5);
    tmat{s} = tform.T;  %Store transformation matrix
    
    %Evaluate a cost function (by transforming and subtracting images)
    imt = imwarp(im2, ra, tform, 'OutputView', ra);    idx = imt > 0;
    %   Mean of squared difference, neglecting values off either frame
    cost_val(s) = mean( (im1(idx) - imt(idx)).^2 );
    
    %BREAK loop IF the threshold is satisfied (save evaluation time)
    if cost_val(s) < cost_thresh; break; end
end

%Use transformation matrix giving lowest cost value
%   Get minimum cost and apply
[cval, si] = min(cost_val);        tmat = tmat{si};

%Filter out small transformations (may be erroneous and uTrack will cover)
if  ( mean((im1(:) - im2(:)).^2) - cval )   < rw  ...
   	||	all( abs(tmat(1:6) - [1,0,0,0,1,0]) < stol);  tmat = eye(3); end

end


%% Shift correction procedure
function xyin = sub_deshift(xyin, ctrans, shd, badt, varargin)
%Procedure depends on type of input (coordinate structure, array, or image)
ns = numel(ctrans); %Number of shift events

if isempty([ctrans.tidx]);  %Check if valid shifts available
    warning('IMAN:XYshift', ['No coodinate transform matrix provided. ',...
        'No XY shift performed.']); return; 
end

%IF shift is in reverse direction, get inverse tranformation matrices
if shd < 0; for s = 1:ns; ctrans(s).tmat = inv(ctrans(s).tmat); end; end

%IF input is a coordinate structure
%--------------------------------------------------------------------------
if isstruct(xyin) && all(isfield(xyin, {'xCoord','yCoord'})) 
    tpts = [[ctrans.tidx] - ctrans(1).stime + 1, length(xyin)+1];
    for ts = 1:ns %FOR each differently shifted time segment
        t_range = setdiff(tpts(ts):(tpts(ts+1)-1), badt); %Exclude bad time
        if isempty(t_range); continue; end %Skip if no coordinates
        %Transform coordinates
        xyt = arrayfun(@(x)([x.xCoord, x.yCoord, ones(length(x.xCoord),1)]...
            *ctrans(ts).tmat), xyin(t_range), 'UniformOutput', false);
        %Get indices for non-negative coordinates (u-Track will error)
        goodi = cellfun(@(x)all(x > 0, 2), xyt,'UniformOutput',false);
        %Re-assign transformed coordinates to structure
        xyt = cellfun(@(x)num2cell(x(:,1:2),1)',xyt,'UniformOutput',false);
        xyt = cat(2,xyt{:});
        [xyin(t_range).xCoord] = deal(xyt{1,:}); 
        [xyin(t_range).yCoord] = deal(xyt{2,:}); 
        %Eliminate negative valued coordinates
        xyin(t_range) = arrayfun(@(xy1,g1) ...
            structfun(@(x)x(g1{1}), xy1, 'UniformOutput', false), ...
            xyin(t_range), goodi);
    end    
    
elseif isnumeric(xyin)     
    sz = size(xyin);  %Get size of input
    
    %IF input appears to be a coordinate array
    %----------------------------------------------------------------------
    if sz(end) == 2 && all(round(xyin(:)) == xyin(:) | isnan(xyin(:)) )
    %xyin as any array must be nCoords x nTime x 2
        tpts = [[ctrans.tidx] - ctrans(1).stime + 1, sz(2)+1];
        xyin = reshape(xyin, prod(sz)/2, 2);
    for ts = 1:ns %FOR each differently shifted time segment
        %Get relevant index range
        xyrng = sz(1)*(tpts(ts)-1)+1:sz(1)*(tpts(ts+1)-1);
        %Transform coordinates
        xyt = [xyin(xyrng,:),ones(length(xyrng),1)]*ctrans(ts).tmat;
        xyin(xyrng,:) = xyt(:,1:2);
    end
    %Reshape xyin to original condition
    xyin = reshape(xyin, sz);
    
    %IF input appears to be an image
    %--------------------------------------------------------------------------
    else    
        %Expect input varargin to provide a time index.  Get time regime
        tidx = varargin{1};    treg = find(tidx < [ctrans.tidx, Inf], 1)-1;
        if treg < 1; return; end    %IF before any shifts, return now
        %Define spatial reference based on orginal image size
        ra = imref2d(sz);
        %Define the transformation
        tform = affine2d(ctrans(treg).tmat);
        %Correct the shift
        xyin = imwarp(xyin, ra, tform, 'OutputView', ra);
    end
end
end

