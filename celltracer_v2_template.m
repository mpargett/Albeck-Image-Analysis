%CELLTRACER_V2_SCRIPT
%   


function [ip, op] = celltracer_v2_template()
%Initialize Celltracer Parameter Structures
[ip, op] = celltracer_v2_partner('initialize');


%% Define Imaging Parameters
%   --------- EDIT THIS SECTION TO DEFINE THE IMAGE DATA ---------
%Data locations
%   Target imaging data file name (include path)
ip.fname = '';  
%   Target save name for processed data files (include path)
ip.sname = '';  

%Data - basic statistics
%   Provide the expected size of imaging data
ip.indsz.xy = [];   %Number of XY positions
ip.indsz.z  = [];   %Number of Z positions
ip.indsz.t  = [];   %Number of Time points
ip.indsz.c  = [];   %Number of color Channels
ip.tsamp    = [];   %Time between samples (MUST BE FILLED)

%Data - detailed definitions
%   Provide the baseline image intensity value (from Camera)
ip.bval = 100;
%   Initialize the properly sized bkg structure
ip.bkg(ip.indsz.xy) = ip.bkg;   
%   Provide background regions or values for each XY position
[ip.bkg(1).reg]    = deal([0, 0; 0, 0]);     %Region (or fixed values)
[ip.bkg(1).fix]   = deal(false);	%Is .reg a fixed intensity value
[ip.bkg(1).dyn]    = deal(false);  	%Is this XY bkg dynamic (sample each frame)?
[ip.bkg(1).altxy]  = deal([]);     	%Index of alternate XY to use for bkg
ip.bkg(1).bkonly = false;  	%Is this XY only for background
%       ***CRITICAL:  Ensure background information is properly defined!
%       We recommend using a background well, and setting all other wells'
%       .altxy field the to the index of the background well.

%   Provide expected (or measured) frame-shift events
ip.xyshift.frame = [];      %Time indices of first shifted frames
ip.xyshift.dx    = [];      %x magnitude of shift (omit for tracking)
ip.xyshift.dy    = [];      %y magnitude of shift (omit for tracking)
ip.xyshift.trackwell = [];  %XY index for use in tracking
ip.xyshift.trackchan = [];  %Channel to use in tracking
ip.xyshift.badtime   = [];  %Time points bad for tracking (optional)
ip.xyshift.shifttype = 'translation';  %Type of shift to perform

%   Provide backup imaging MetaData:
%       (** = Critical for all runs, ++ = Only for objbias)
    %   Parameters of the objective
    %   ---------------------------------------------------------------
    ip.bkmd.obj.Desc       =   'Apo_20x';  %Description of Objective
    ip.bkmd.obj.Mag        =   20;         %Magnification ++
    ip.bkmd.obj.WkDist     =   1;          %Working Distance (mm) ++
    ip.bkmd.obj.RefIndex   =   1;          %Index of Refraction ++
    ip.bkmd.obj.NA         =   0.75;       %Numerical Aperture ++
    %   Parameters of the camera
    %   ---------------------------------------------------------------
    ip.bkmd.cam.Desc       =   'Zyla5';    %Descriptor of camera
    ip.bkmd.cam.PixSizeX   =   0.65;       %Pixel horizontal size (um) **
    ip.bkmd.cam.PixSizeY   =   0.65;       %Pixel vertical size (um) **
    ip.bkmd.cam.PixNumX    =   1280;       %Number of Pixels in X **
    ip.bkmd.cam.PixNumY    =   1080;       %Number of Pixels in Y **
    ip.bkmd.cam.BinSizeX   =   2;          %Number of Pixels per bin, X
    ip.bkmd.cam.BinSizeY   =   2;          %Number of Pixels per bin, Y
    %   Parameters of the exposures
    %   ---------------------------------------------------------------
    %   Declare here names for the color channels, in proper order.  
    %       All other entries in ip.bkmd.exp must be in the same order.
    ip.bkmd.exp.Channel    =   {'', '', ''};  %Name(s) of Channel(s)
    
    %   ----- OPTIONAL - Only for use with Spectral Unmixing -----
    %   Declare the Filters for each Channel ++
    ip.bkmd.exp.Filter     =   {'', '', ''};
    %   Declare the target Fluorophores for each Channel ++
    ip.bkmd.exp.FPhore     =   {'', '', ''};
    %       -> See iman_naming for valid Filter and Fluorophore Names    
    %   Declare the Channel names associate with each FRET pair
    %       (only needed for pixel-by-pixel spectral unmixing of FRET)
    ip.bkmd.exp.FRET       =   {};% struct('EKAR', {{'CFP','YFP'}}); % {[1,2]}
    ip.bkmd.exp.Light      =   'SOLA';           %Light source name ++
    ip.bkmd.exp.Exposure   =   [0, 0, 0];  %Exposure times (ms) ++
    ip.bkmd.exp.ExVolt     =   [0, 0, 0];     %Relative voltage (0-100) ++
    %   IF using the SPECTRAX light source, also define the following:
%     ip.bkmd.exp.ExLine     =   {[],[],[]};   %Excitation line
%     ip.bkmd.exp.ExWL       =   {[],[],[]};   %Wavelength
    %   ----- ----------------------------------------------- -----


%Show final image parameters (Comment out the following to avoid printing)
display(ip);
fprintf('\nip.indsz = \n\n');  display(ip.indsz); 
fprintf('\nip.bkmd = \n\n');   display(ip.bkmd);
%   ------ ------ END OF IMAGE DATA DEFINITION SECTION ------ ------


%% Define Operation Parameters
%   -------- EDIT THIS SECTION TO SETUP THE PROCESSING RUN --------
%Procedure scope - target ranges of the data to process
op.cind     = [];       %Indices of Channels to process (names or indices)
op.xypos    = [];       %Indices of XY positions to process
op.trng     = [];       %Start and End indices of Time to process
op.nW       = 1;        %Number of parallel workers to use
%   Optional procedures
op.objbias  = true;     %Correct for objective view bias
op.unmix    = false;  	%Linearly unmix color channel cross-talk
op.fixshift = false;  	%Correct any indicated frameshifts

%   Indicate which fields of bkmd should override any other MetaData source
%   (use when recorded MetaData are corrupted)
op.mdover   = {'exp', 'ExVolt'; 'exp', 'Exposure';'exp','Filter'; ...
    'exp','Channel'; 'exp','FPhore'};

%Procedure settings ***CRITICAL - Review carefully
%   Segmentation settings
op.seg.chan = [];       %Color channel on which to segment  !!
%   A second seg channel may be given (opposite nuc/cyt of first), to
%   filter cyto masks. Use op.seg.cyt to specify the FIRST seg channel
op.seg.cyt  = false; 	%Segment on cytoplasmic signal      !!
%       The following defaults are typical of MCF10As at 20x magnification
op.seg.maxD = 25;       %Maximum nuclear diameter (um)
op.seg.minD = 8;      	%Minimum nuclear diameter (um)
%       It is recommended to scale imaging for nuclei (segmentation
%       targets) to be greater than 10 pixels in diameter.
op.seg.maxEcc = 0.9;    %Maximum eccentricity [0-1] (opposite of circularity)
op.seg.Extent = [0.65, 0.85]; %Minimum fraction of bounding box filled 
%   Extent range is [0-1], value for a cirlce is PI/4
op.seg.maxSmooth = 0.1;     %Maximum fraction of convex hull without mask 
op.seg.sigthresh = []; %[Optional] Minimum intensity of 'good' signal
%   If you use sigthresh, pre-subtract any camera baseline
op.seg.hardsnr = false;  %Typically kept to FALSE.  TRUE makes the signal 
%   threshold 'hard', enforcing cutoff of any pixels below it.
op.seg.nrode = 0;   %[Optional] Approx. number of pixels to erode nuc mask
%   NOTE: op.seg.nrode is only to improved nearby cell separation - it is
%   undone in final masks.  Use op.msk.nrode to adjust final masks.

%    Masking settings
op.msk.rt = {};             %Pre-averaging channel ratios to take {{'',''}}
op.msk.storemasks = false;  %Save all segmentation masks
op.msk.saverawvals = true;  %Save all raw valcube entries (pre-tracking)
%   IF additional aggregation functions are desired, define here
%       For example:
%           op.msk.aggfun.name = 'var';
%           op.msk.aggfun.chan = {'CNAME1', 'CNAME2', ...};
%           op.msk.aggfun.loc = {'Nuc','Cyt'};
%           op.msk.aggfun.fun = @var;
op.msk.nrode = 0;   %[Optional] Approx. # pixels to erode final nuc mask
op.msk.cgap = 0;    %[Optional] Approx. # pixels to expand nuc-cyt mask gap
op.msk.cwidth = 0;  %[Optional] Approx. # pixels to expand cyto ring width
op.msk.nfilt = false; %[Optional] TRUE to filter cyto mask by nuc channel
%   Removes cyto mask where nuc signal > 5th percentile of nuc channel
%   IF significant varied cyto staining in nuc channel, do not use
op.msk.cfilt = false; %[Optional] TRUE to filter cyto mask by cyt channel
%   Removes cyto mask where cyto signal is at background (below sigthresh)

% Tracking settings (for utrack)
%   Typically, these do not require editing
op.trk.linkrad = 10;        % Radius (in um) to consider for frame-to-frame 
                            %   linking of cell movement
op.trk.minTrkLength = 5;    % Minimum number of consecutive frames required
                            %   for a track to be used in gap closing
op.trk.gaprad = 15;         % Radius (in um) for gap closing
op.trk.gapwin = 30;         % Time window (in minutes) for gap closing
op.trk.msRadMult = 1;       % Multiplier for search radius when finding 
                            %   merges and splites (Keep at 1 if tracking
                            %   splits is not critical. If increased, run
                            %   validation GUI to check flagged splits for
                            %   real/fake.

%   Display settings - select which options to show while running
op.disp.meta     = false;       %Final MetaData to be used
op.disp.samples  = false;       %Sparse samples of segmentation
op.disp.shifts   = false;       %Shifted images
op.disp.warnings = false;       %Procedural warning messages

%Show final operation parameters (uncomment to show)
% display(op);
% fprintf('\nop.seg = \n\n');  display(op.seg); 
% fprintf('\nop.disp = \n\n');  display(op.disp); 
%   ------ ------ END OF PROCESSING RUN SETUP SECTION ------ ------

%   Do NOT edit the following parameter transfers
op.cname    = ip.bkmd.exp.Channel;  %Copy channel names to op
op.msk.fret = ip.bkmd.exp.FRET;     %Copy FRET channel names to op


%% Run image analysis procedure
%Pre-Run validation
[ip, op] = celltracer_v2_partner('validate', ip, op);

%Optional segmentation check
%   Run the following commented function to view segmentation based on
%   current parameters.  Type help iman_segcheck for more details.  op
%   settings can then be adjusted and iman_segcheck rerun to test. Do not
%   reset pst to empty - it is needed to avoid reloading of the dataset.
% pst = [];
% pst = iman_segcheck(ip, op, 't',1, 'xy',1, 'pastinfo',pst);
%
%Optional tracking check
%   Run the following, after editing the time range ('t') and potentially
%   the range of frames ('frm') to plot ('splits' instructs the routine to
%   plot the frame before and after any track split).  Reset d.tF to empty
%   to rerun u-track if op.trk parameters have been changed.  Do not reset
%   d = [] - it is needed to avoid reloading of the dataset and rerunning
%   segmentation.
% d = [];       d.tF = [];
% d = iman_trkcheck(ip, op, d, 't',1:100, 'npth',10, 'frm','splits');
%
%  After running, can use the following to run segview
%   iman_segview(ip.sname,1,1,'printid',true,'idcolor',[0,0,0]);


%% Main procedure
[dao, GMD, dmx] = iman_celltracer(ip, op);



end





