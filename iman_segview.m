%IMAN_SEGVIEW
%   Segmentation viewer for Image Processing diagnostics.  Use this viewer
%   to overlay processed masked onto images, creating a movie.  One channel
%   and XY position may be viewed per movie, and the time points to use may
%   be specified.  The viewer relies on information stored in the "_Global"
%   file saved with each processing run, which will be loaded
%   automatically.
%
%   [segstack] = iman_segview(bpath, chan, xy, ...)
%       returns the image stack, segstack, which may be converted into a
%   movie or played using implay.  
%
%   INPUTS
%       bpath - [string] Path to your processed data, plus the base name of
%                   the processed data files (do not include _xy##.mat).
%       chan  - [scalar, or string] Index or name of channel to show.
%       xy    - [scalar] Index of the the XY to view .
%   
%   OPTIONAL INPUTS
%   Optional inputs may be passed as Name/Value pairs (two sequential
%   inputs per parameter to set: first the name, then the new value).
%       dao      - Previously generated data access object for images.
%       trange   - [Vector] Time points to use (e.g. [1:50])
%       nucmask  - [Logical] Set TRUE to show nuclear masks.
%       cytmask  - [Logical] Set TRUE to show cytoplasm masks.
%       printid  - [Logical] Set TRUE to print cell IDs.
%       imcolor  - [1x3 array] RGB color for the image channel
%       nuccolor - [1x3 array] RGB color for nuclear masks
%       cytcolor - [1x3 array] RGB color for cytoplasm masks
%       idcolor  - [1x3 array] RGB color cell IDs
%       idsize   - [Scalar or 1x2] Size of numbers to print for IDs.
%           Default: 2.  1 gives numbers at 3 pixels by 5 pixels.
%       nw       - [Scalar] Number of workers to use in parallel.  Only
%           recommended if more than 200 time points, or very large images
%           (e.g. 1x1 binning).
%       keeppool - [Logical] Set TRUE to keep parallel pool open (e.g. for
%           faster future runs if generating multiple movies).
%
%
%   EXAMPLE USAGE
%   Calling for Channel 1 and XY 22, for Time Points 10-30 (Note that the
%   actual data file would be called 'testfile_xy22.mat'. Only include the
%   base name when specifying bpath.): 
%       bpath = 'L:\Processed Data\Test Folder\testfile';
%       im_mask_viewer(bpath, 1, 22, 'trange', 10:30);

%           
% -----------------------------------------------------------------------------
% VERSION 1.0


function [segstack] = iman_segview(bpath, chan, xy, varargin)

%% Parse Inputs
%   Default parameters
p.dao = []; p.trange = [];  %No Dao and all time points
p.playnow = true;           %Play when finished
p.nucmask = true;           %Plot nuclear masks
p.cytmask = false;          %Do not plot cytoplasm masks
p.printid = false;          %Do not print cell IDs
p.imcolor   = [0, 0, 0.7];      %Color to plot image channel (Blue)
p.nuccolor  = [0, 1, 0.25];     %Color to plot nuc masks (Yellow/Green)
p.cytcolor  = [1, 0.3, 0.3];    %Color to plot cyt masks (Red)
p.idcolor   = [1,1,1];          %Color to print Cell IDs (White)
p.idsize    = 2;            %Text size to print Cell IDs
p.nw        = 1;            %Number of workers for parallel
p.keeppool  = false;        %Do not keep parallel pool

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end


%% Load processed data
dh = '0';   %Allow for leading zero padded XY in loading path 
%Load processed data file 
fname = [bpath,'_xy', dh(xy<10), num2str(xy), '.mat'];
if exist(fname,'file');    pdat = load(fname, 'masks', 'valcube', 'vcorder');
    xi = strcmpi(pdat.vcorder, 'XCoord');   %Get X channel
    yi = strcmpi(pdat.vcorder, 'YCoord');   %Get Y channel
    if isempty(pdat.masks); error('IMAN:NoMask', ['No stored masks were ',...
            'found for ', fname, '.']); end
else     error(['Data File ',fname,' does not exist. Check your path.'])
end

%Load global metadata file
fname = [bpath, '_Global.mat'];
if exist(fname,'file');    glb = load(fname);
else     error(['Global file (metadata) ', ...
            fname,' does not exist. Check your path.'])
end


%% Validations
%Verify input parameters
%   Time series requested
if isempty(p.trange);   p.trange = glb.op.trng(1) : glb.op.trng(end);
elseif numel(p.trange) == 2; p.trange = p.trange(1) : p.trange(2); end
%   IF requested time is out of range, throw error
if any(p.trange > glb.op.trng(end) | p.trange < glb.op.trng(1) );
    error('IMAN:TRange',['Requested time range is greater than ',...
        'maximum time in processed data']);    
end
%   Channel requested
if ischar(chan); cid = strcmpi(glb.GMD.exp.Channel, chan);
else                cid = chan;      end
%   If bad channel, throw error
if isempty(cid) || cid > numel(glb.GMD.exp.Channel)
    error('IMAN:Channel', ['Channel ', num2str(chan), 'not found.']);
end

%Ensure java path availability for ND2 files
if ~isempty(regexpi(glb.p.fname, '.nd2$')); 
            pth = iman_assertjavapaths(); end

%% Generate frames
%Initialize variables
nt = numel(p.trange);
%   Pre-define structuring element for dilation
st2 = strel('disk',2);
%   Pre-generate number images for IDs
if p.printid; n = pixnums(p.idsize);  else n = {}; end %Shenanigans

%Check and initialize for Parallel Processing
if p.nw > 1;  nper = pp_init(p, nt, pth);  else nper = nt;  end

%List typical warnings to suppress
wid = {'MATLAB:Java:DuplicateClass', 'images:initSize:adjustingMag'};
if p.nw > 1;  %IF parallel, turn off on all workers
    for s = 1:length(wid)
        pctRunOnAll(['warning(''off'', ''',wid{s},''');']);
    end
else    for sw = 1:length(wid); warning('off', wid{sw}); end
    %   ^Or else, just loop for the client
end

%Initialize parallelized variables
segstack = cell(p.nw,1);  imname = glb.p.fname;  indsz = glb.p.indsz;
tst = glb.op.trng(1);  npix = [glb.GMD.cam.PixNumX, glb.GMD.cam.PixNumY];
dao = p.dao;
parfor ps = 1:p.nw
    %Load Data Access Object
    if isempty(dao) || ~dao.info.ThreadSafe;  
            daop = iman_imageaccess(imname, indsz);
    else    daop = dao; 
    end
    
    %Inialize output frame structure
    segstack{ps} = struct('cdata', cell(nper(ps),1), 'colormap', []);
    k = 0;  %Initialize output frame index
    for st = 1:nper(ps);
        %Get current time index (as start + already completed)
        s = p.trange( sum( nper(1:(ps-1)) ) + st );     %#ok<PFBNS>
        
        k = k + 1;  %Iterate output frame index
        %Align time if processing doesnt start at 1.
        mskframe = s - (tst-1);
        %Get the frame
        IM = iman_getframe(daop, [s,cid,xy,1]);
        IM = IM./prctile(IM(:),95);             %Scale for visibility
        
        %Nuclear Masks
        if p.nucmask
            IM = bsxfun(@times, IM, shiftdim(p.imcolor(:),-2));
            %   Get nuclear masks
            nm = bwunpack(pdat.masks(mskframe).nuc, npix(2)); %#ok<PFBNS>
            %   Dilate and subtract for a boundary
            nm = imdilate(nm,st2) & ~nm;
            %   Set mask boundaries to desired color
            IM(repmat(nm,1,1,3)) = ones(nnz(nm),1)*p.nuccolor;
        end
        
        %Cytoplasms Masks
        if p.cytmask    %Repeat plotting as with nuclear, if called
            %   Get masks, boundaries, and write to image
            nm = bwunpack(pdat.masks(mskframe).cyt, npix(2));
            nm = imdilate(nm,st2) & ~nm;
            IM(repmat(nm,1,1,3)) = ones(nnz(nm),1)*p.cytcolor;
        end
        
        %Cell IDs
        if p.printid    %Only plot IDs if requested
            gi = find(~isnan(pdat.valcube(:,s,xi))); %Get good indices
            %   Get number images for each ID, set up offsets for numbers
            y = num2pix(gi,n); nn = [-1/2, 1/2]; nh = size(y{1},1).*nn;
            %   Initialize the ID mask for the whole image
            idmask = size(IM);  idmask = false(idmask(1:2));
            for sg = 1:numel(y) %FOR each ID
                %   Get x bounds for insert region
                xb = floor(pdat.valcube(gi(sg), s, xi) + size(y{sg},2).*nn);
                %   Get y bounds for insert region
                yb = floor(pdat.valcube(gi(sg), s, yi) + nh);
                %Enforce image bounds
                if xb(1) < 1;           xb = xb + 1 - xb(1);        end
                if xb(2) > npix(1);     xb = xb + npix(1) - xb(2);     end
                if yb(1) < 1;           yb = yb + 1 - yb(1);        end
                if yb(2) > npix(2);     yb = yb + npix(2) - yb(2);     end
                %   Place number images into ID mask
                idmask((yb(1)+1):yb(2), (xb(1)+1):xb(2)) = y{sg};
            end
            %Write text to image
            IM(repmat(idmask,1,1,3)) = ones(nnz(idmask),1)*p.idcolor;
        end
        
        %Store frame
        segstack{ps}(k).cdata = IM;
    end     %END Serial Loop
end     %END Parallel Loop

%Compile parallel processed stacks
segstack = cell2mat(segstack);

%Open in MATLAB movie player if desired now
if p.playnow;    implay(segstack,5);   end

%Clear parallel pool
if ~p.keeppool;     delete(gcp('nocreate'));    end


end


%Subfucntion: Initialize Parallel Processing
function nper = pp_init(p, nt, pth)
    if isempty(gcp('nocreate')) %Initial Parallel pool, if not present
        parpool('local', p.nw, 'IdleTimeout', Inf);     end
    ipw = ceil( nt./p.nw );               %Approx. images per Worker
    nper = ipw*ones(1,p.nw);  dif = p.nw*ipw - nt - 1; %Fix overestimate
    nper(end-dif:end) = nper(end-dif:end) - 1;  %Minus 1 excess per Worker
    
    %Assert java class paths on workers
    if iscell(pth.jc); runtext = ['javaclasspath({''',pth.jc{1}];
        for s = 2:numel(pth.jc); runtext = [runtext,''',''',...
                pth.jc{s}]; end;  runtext = [runtext,'''});']; %#ok<AGROW>
    else runtext = ['javaclasspath ',pth.jc,';'];
    end;     pctRunOnAll(runtext);
end

%Subfunction: Convert number to pixel pattern
function y = num2pix(x,n)
%Ensure cell encapsulation of x
if ~iscell(x); x = num2cell(x); end
%Get digit indices
%   First split digits from each number in x
di = cellfun(@(xx) cell2mat(arrayfun(@(a)str2double(a)+1, ...
                    num2str(xx), 'Un', 0)), x, 'Un', 0);
%   Then add spacer indices in between each
di = cellfun(@(xx) reshape(cat(1, xx, ones(size(xx)).*11), ...
                    1, numel(xx).*2), di, 'Un', 0);
%Concatenate pixel patterns (neglect trailing spacer)
y = cellfun(@(yy)cat(2, n{yy(1:end-1)}), di, 'Un', 0);
end

%Subfunction:  Define pixel shapes for numbers 0-9
function n = pixnums(scl)
%Define pixel maps for numbers
n{1}  = [1,1,0;1,0,1;1,0,1;1,0,1;0,1,1]; %0
n{2}  = [0,1,0;1,1,0;0,1,0;0,1,0;0,1,0]; %1
n{3}  = [1,1,1;0,0,1;0,1,1;1,0,0;1,1,1]; %2
n{4}  = [1,1,1;0,0,1;1,1,1;0,0,1;1,1,1]; %3
n{5}  = [1,0,1;1,0,1;1,1,1;0,0,1;0,0,1]; %4
n{6}  = [1,1,1;1,0,0;1,1,0;0,0,1;1,1,0]; %5
n{7}  = [0,1,1;1,0,0;1,1,1;1,0,1;0,1,1]; %6
n{8}  = [1,1,1;0,0,1;0,1,0;1,0,0;1,0,0]; %7
n{9}  = [1,1,1;1,0,1;1,1,1;1,0,1;1,1,1]; %8
n{10} = [0,1,1;1,0,1;1,1,1;0,0,1;1,1,0]; %9
n{11} = [0;0;0;0;0];    %Space
%   Convert to logicals
n = cellfun(@logical, n, 'Un', 0);

%Apply scaling for larger numbers as needed
if exist('scl','var') && any(scl > 1); %IF block scaling
    if length(scl) == 1; scl = [scl,scl]; end
    b1 = ones(scl);  b0 = zeros(scl);
    for s = 1:numel(n)
        c = cell(size(n{s}));
        [c{n{s}}]  = deal(b1);
        [c{~n{s}}] = deal(b0);
        n{s} = cell2mat(c);
    end
end

end