%IMAN_ICELLCOUNT
%   Estimate cell count by total frame intensity, using saved processed
%   data (via the iman_celltracer pipeline).
%
%   [cce, scc] = iman_icellcount(bn, ...)
%       returns the intensity-based cell count estimate (cce) as well as
%       the segmented cell count (scc, optional), for the experiment
%       located by bn, the base path name.  bn must indicate the path and
%       base file name for the processed data, for example, 'L:\Processed
%       Data\ExperimentName\BaseName', excluding the '_xy01.mat' etc. found
%       on the processed data files.
%
%   Additional paraemeters may be specified as Name/Value pairs:
%
%   xy - XY position(s) to include (vector)
%   t - Time point(s) to include (vector)
%   c - Channel to use (scalar)
%   z - Z position(s) to include (vector, not yet supported)
%   trim - Percentiles at which to trim tails of intensity distribution 
%           (1x2 vector, e.g. [5,95])
%   nuc - Logical, TRUE if channel shows nuclei (FALSE not yet supported)
%   cthresh - Threshold number of segmented cells for a 'good' frame. Cell
%       count estimates will be set to zero for any frames with a segmented
%       cell count below the threshold (default 10).
%

function [cce, scc] = iman_icellcount(bn, varargin)
%Default parameters
p = struct('t',1,'c',1,'xy',1,'z',1,'trim', [], 'nuc', true, 'cthresh', 10);

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%Assert that image processing paths are correct
iman_assertjavapaths;

%% Process data
%Get global data for experiment
gd = load([bn,'_Global.mat']);

%Get DAO for image data
%   Ensure filename available
if exist(gd.p.fname,'file'); fn = gd.p.fname;
else  warning('Image file not found. Select via browser.'); fn = uigetfile;
end
%   Access image file
dao = iman_imageaccess(fn, gd.p.indsz);
%Get frame intensity sums
y = iman_framestats(dao, 'xy', p.xy, 't', p.t, 'c', p.c, 'z', p.z, ...
    'trim', p.trim, 'sum', true, 'mean', false, 'var', false);

%Get average cell intensities
%   Get Channel name(s)
nucyt = {'_Nuc','_Cyt'};
cn = [gd.p.bkmd.exp.Channel{p.c},nucyt{2-p.nuc}];  an = 'nArea';
dh = '0';  xyn = @(x)[dh(x<10),num2str(x)];
%   Initialize outputs
scm = nan(numel(p.xy), numel(p.t), numel(p.c), numel(p.z)); 
sca = scm; scc = scm;  bkg = scm;
for sx = 1:numel(p.xy)
    for sz = 1:numel(p.z)  %FIXME - doesn't actually handle Z slices
        %   Load processed data
        d = load([bn,'_xy', xyn(p.xy(sx)), '.mat'], ...
            'valcube','vcorder','bkg');
        %   Store background
        bkg(sx,:,:,sz) = d.bkg;
        %   Get Channel index for intensity data
        ci = strcmpi(cn, d.vcorder);  ci = ci(1:size(d.valcube,3));
        %   Take average intensity
        scm(sx,:,:,sz) = nanmean(d.valcube(:,p.t,ci),1);
        %   Get segmented cell count
        scc(sx,:,:,sz) = sum(~isnan(d.valcube(:,p.t,ci)),1);
        %Get area per cell (average)
        if p.nuc
            %   Get Channel index for area data
            ai = strcmpi(an, d.vcorder);  ai = ai(1:size(d.valcube,3));
            %   Take average area
            sca(sx,:,:,sz) = nanmean(d.valcube(:,p.t,ai),1);
        else %FIXME - Something for cyto? (Not immediately necessary)
        end
    end
end

%Calculate cell count estimates (subtracts background)
cce = (y.sum - (bkg.*gd.GMD.cam.PixSizeX.*gd.GMD.cam.PixSizeY))./(scm.*sca);
%   Adjust for low segmented cell counts
cce(scc < p.cthresh) = 0;


end