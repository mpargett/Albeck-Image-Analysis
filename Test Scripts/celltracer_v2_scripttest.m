%CELLTRACER_V2_SCRIPTTEST
%   


function [ip, op] = celltracer_v2_scripttest()
%Initialize Celltracer Parameter Structures
[ip, op] = celltracer_v2_partner('initialize');


%% Define Imaging Parameters
%   --------- EDIT THIS SECTION TO DEFINE THE IMAGE DATA ---------
%Data locations
%   Target imaging data file name (include path)
ip.fname = 'L:\imageData\ProcessingValidation\2014-05-12-FRA1 C.9.nd2';  
%   Target save name for processed data files (include path)
ip.sname = 'celltracer_v2_test';  

%Data - basic statistics
%   Provide the expected size of imaging data
ip.indsz.xy = 10;   %Number of XY positions
ip.indsz.z  = 1;   %Number of Z positions
ip.indsz.t  = 6;   %Number of Time points
ip.indsz.c  = 3;   %Number of color Channels

%Data - detailed definitions
%   Provide the baseline image intensity value (from Camera)
ip.bval = 100;
%   Initialize the properly sized bkg structure
ip.bkg(ip.indsz.xy) = ip.bkg;   
[ip.bkg.fix]    = deal(false);  %Default fix is false
[ip.bkg.dyn]    = deal(true);   %Default dyn is true
[ip.bkg.bkonly] = deal(false);  %Default bkonly is false
%   Provide background regions or values for each XY position
ip.bkg(1).reg =  [50, 1; 230, 45];          
ip.bkg(2).reg =  [1260, 1; 1280, 130];      ip.bkg(2).dyn = false;
ip.bkg(3).reg =  [1000, 175, 275];          ip.bkg(3).fix = true;  
                                            ip.bkg(3).dyn = false;
ip.bkg(4).reg =  [1, 1; 50, 35];  
ip.bkg(5).reg =  [930, 700; 1070, 830]; 
ip.bkg(6).reg =  [10, 680; 105, 920]; 
ip.bkg(7).reg =  [10, 680; 105, 920];       ip.bkg(7).altxy = 8;
ip.bkg(8).reg =  [1183, 992; 1280, 1080];   ip.bkg(8).bkonly = true;
ip.bkg(9).reg =  [1245, 1; 1280, 75]; 
ip.bkg(10).reg = [1, 960; 80, 1080];

%   Provide expected (or measured) frame-shift events
ip.xyshift.frame = [];      %Time indices of first shifted frames
ip.xyshift.dx    = [];      %x magnitude of shift (omit for tracking)
ip.xyshift.dy    = [];      %y magnitude of shift (omit for tracking)
ip.xyshift.trackwell = [];  %XY index for use in tracking
ip.xyshift.trackchan = [];  %Channel to use in tracking
ip.xyshift.badtime   = [];  %Time points bad for tracking (optional)
ip.xyshift.shifttype = 'translation';  %Type of shift to perform

%   Provide backup imaging MetaData:
    %   Parameters of the objective
    %   ---------------------------------------------------------------
    ip.bkmd.obj.Desc       =   'Apo_20x';  %Description of Objective
    ip.bkmd.obj.Mag        =   20;         %Magnification
    ip.bkmd.obj.WkDist     =   1;          %Working Distance (mm)
    ip.bkmd.obj.RefIndex   =   1;          %Index of Refraction
    ip.bkmd.obj.NA         =   0.75;       %Numerical Aperture
    %   Parameters of the camera
    %   ---------------------------------------------------------------
    ip.bkmd.cam.Desc       =   'Zyla5';    %Descriptor of camera
    ip.bkmd.cam.PixSizeX   =   0.65;       %Pixel horizontal size (um)
    ip.bkmd.cam.PixSizeY   =   0.65;       %Pixel vertical size (um)
    ip.bkmd.cam.PixNumX    =   1280;       %Number of Pixels in X
    ip.bkmd.cam.PixNumY    =   1080;       %Number of Pixels in Y
    ip.bkmd.cam.BinSizeX   =   2;          %Number of Pixels per bin, X
    ip.bkmd.cam.BinSizeY   =   2;          %Number of Pixels per bin, Y
    ip.bkmd.cam.tsamp      =   5;          %Sampling time (min)
    %   Parameters of the exposures
    %   ---------------------------------------------------------------
    %   Declare here names for the color channels, in proper order.  
    %       All other entries in ip.bkmd.exp must be in the same order.
    ip.bkmd.exp.Channel    =   {'CFP', 'YFP', 'RFP'};  %Name(s) of Channel(s)
    %   Declare the Filters for each Channel
    ip.bkmd.exp.Filter     =   {'Filter_CFP', 'Filter_YFP', 'Filter_Cherry2'};
    %   Declare the target Fluorophores for each Channel
    ip.bkmd.exp.FPhore     =   {'mTurq2', 'YPet', 'mCherry'};
    %       -> See iman_naming for valid Filter and Fluorophore Names    
    %   Declare the Channel names associate with each FRET pair
    ip.bkmd.exp.FRET       =   struct('EKAR', {{'CFP','YFP'}}); % {[1,2]}
    ip.bkmd.exp.Light      =   'SOLA';           %Light source name
    ip.bkmd.exp.Exposure   =   [800, 500, 800];  %Exposure times (ms)
    ip.bkmd.exp.ExVolt     =   [15, 15, 30];     %Relative voltage (0-100)
    %   IF using the SPECTRAX light source, also define the following:
%     ip.bkmd.exp.ExLine     =   {[],[],[]};   %Excitation line
%     ip.bkmd.exp.ExWL       =   {[],[],[]};   %Wavelength

%Show final image parameters (Comment out the following to avoid printing)
display(ip);
fprintf('\nip.indsz = \n\n');  display(ip.indsz); 
fprintf('\nip.bkmd = \n\n');   display(ip.bkmd);
%   ------ ------ END OF IMAGE DATA DEFINITION SECTION ------ ------


%% Define Operation Parameters
%   -------- EDIT THIS SECTION TO SETUP THE PROCESSING RUN --------
%Procedure scope - target ranges of the data to process
op.cind     = 'YFP';%1:ip.indsz.c;    %Indices of Channels to process (names or indices)
op.xypos    = [1,2,3,7,8];     %Indices of XY positions to process
op.trng     = [2,6];       %Start and End indices of Time to process
op.nW       = 4;        %Number of parallel workers to use
op.unmix    = false;  	%Linearly unmix color channel cross-talk
op.fixshift = false;  	%Correct any indicated frameshifts
%   Indicate which fields of bkmd should override any other MetaData source
op.mdover   = { 'exp', 'Channel';   'exp', 'Exposure';  ...
                'exp', 'ExVolt';    'exp', 'Filter';    'exp', 'FPhore'};

%Procedure settings
%   Segmentation settings
op.seg.chan = 'YFP';       %Color channel on which to segment
op.seg.cyt  = true; 	%Segment on cytoplasmic signal
%       The following defaults are typical of MCF10As at 20x magnification
op.seg.maxD = 25;       %Maximum nuclear diameter (pixels)
op.seg.minD = 8;       %Minimum nuclear diameter (pixels)
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
%   NOTE: op.seg.nrode is only to improve nearby cell separation - it is
%   undone in final masks.  Use op.msk.nrode to adjust final masks.

op.msk.rt = {};%{{'CFP','YFP'}};  %Pre-averaging channel ratios to take
op.msk.storemasks = false;    %Save all segmentation masks
op.msk.saverawvals = true;    %Save all raw valcube entries (pre-tracking)
%   IF additional aggregation functions are desired, define here
%       For example:
%     op.msk.aggfun.name = 'var';
%     op.msk.aggfun.chan = {'YFP'};
%     op.msk.aggfun.loc = {'Nuc'};
%     op.msk.aggfun.fun = @var;
op.msk.nrode = 1;   %[Optional] Approx. # pixels to erode final nuc mask
op.msk.cgap = 0;    %[Optional] Approx. # pixels to expand nuc-cyt mask gap
op.msk.cwidth = 0;  %[Optional] Approx. # pixels to expand cyto ring width

% Tracking settings (for utrack)
op.trk.linkrad = 10;        % Radius (in um) to consider for frame-to-frame 
                            %   linking of cell movement
op.trk.gaprad = 15;         % Radius (in um) for gap closing
op.trk.nnwin = 10;          % Time window (in frames) to consider when 
                            %   evaluating Nearest Neighbors for density
op.trk.gapwin = 5;          % Time window (in frames) for gap closing
op.trk.minTrkLength = 2;    % Minimum number of consecutive frames required
                            %   for a track to be used in gap closing
op.trk.gapPenalty = 1.5;    % Added penalty in gap closing for longer gaps
                            %   n+1 frame gap is penalized by gapPenalty^n
op.trk.ampRat = [0.5, 1.5]; % Amplitude Ratio bounds to use for Merge/Split
                            %   [Min, Max] of AmpAfter / sum(AmpsBefore)
op.trk.ampgap = 5;          % Time Gap (in frames) around Merge/Split to
                            %   neglect when calculating Amp Ratios
op.trk.ampwin = 5;          % Time Window (in frames) around Merge/Split to
                            %   use for Amp Ratios (after skipping Gap)
op.trk.motionType = 1;      % Flag for motion model to use: 0 - Brownian, 
                            %   1 - Linear, 3 - Linear w/reverse (uTubes)
op.trk.bStdMult = 4;       % Standard Deviation multiplier to calculate 
                            %   Brownian motion search regions
op.trk.lStdMult = 4;       % Standard Deviation multiplier to calculate 
                            %   Linear motion search regions 
op.trk.msRadMult = 2;       %Multiplier for search radius when finding 
                            %   merges and splites (Keep at 1 if tracking
                            %   splits is not critical. If increased, run
                            %   validation GUI to check flagged splits for
                            %   real/fake.

%   Display settings - select which options to show while running
op.disp.meta     = true;       %Final MetaData to be used
op.disp.samples  = true;       %Sparse samples of segmentation
op.disp.shifts   = true;       %Shifted images
op.disp.warnings = false;%true;       %Procedural warning messages

%Show final operation parameters
display(op);
fprintf('\nop.seg = \n\n');  display(op.seg); 
fprintf('\nop.disp = \n\n');  display(op.disp); 
%   ------ ------ END OF PROCESSING RUN SETUP SECTION ------ ------
%   Do NOT edit the following parameter transfers
op.cname    = ip.bkmd.exp.Channel;      %Copy channel names to op
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
% pst = iman_segcheck(ip, op, 't',5, 'xy',1, 'pastinfo',pst);
%
% d = [];       d.tF = [];
% d = iman_trkcheck(ip, op, d, 't',1:100, 'npth',10, 'frm','splits');
%
%  After running, can use the following to run segview
%   iman_segview(ip.sname,1,1,'printid',true,'idcolor',[0,0,0]);
%
%   Run parameter optimization (Abandoned)
%   popt.seg = struct('maxEcc', [0.6,1]);
%   z = iman_spo(ip, op, 'seg', popt.seg, 'np', 64);

%Main procedure
[dao, GMD, dmx] = iman_celltracer(ip, op);



end





