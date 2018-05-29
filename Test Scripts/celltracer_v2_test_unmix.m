%CELLTRACER_V2_TEST_UNMIX
%   Tests the spectral unmixing function of the Celltracer image processing
%   pipeline, using an ND2 file for source images.


function [ip, op] = celltracer_v2_test_unmix()
%Initialize Celltracer Parameter Structures
[ip, op] = celltracer_v2_partner('initialize');


%% Define Imaging Parameters
%   --------- EDIT THIS SECTION TO DEFINE THE IMAGE DATA ---------
%Data locations
%   Target imaging data file name (include path)
ip.fname = ['L:\imageData\ProcessingValidation\',...
    '2016-03-15_cropMP.nd2'];    
%   Target save name for processed data files (include path)
ip.sname = 'celltracer_v2_test_unmix';  

%Data - basic statistics
%   Provide the expected size of imaging data
ip.indsz.xy = 1;   %Number of XY positions
ip.indsz.z  = 1;   %Number of Z positions
ip.indsz.t  = 31;   %Number of Time points
ip.indsz.c  = 5;   %Number of color Channels

%Data - detailed definitions
%   Provide the baseline image intensity value (from Camera)
ip.bval = 100;
%   Initialize the properly sized bkg structure
ip.bkg(ip.indsz.xy) = ip.bkg;   
[ip.bkg.fix]    = deal(false);  %Default fix is false
[ip.bkg.dyn]    = deal(true);   %Default dyn is true
[ip.bkg.bkonly] = deal(false);  %Default bkonly is false
%   Provide background regions or values for each XY position
ip.bkg(1).reg =  [616, 14; 747, 76];

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
    ip.bkmd.cam.tsamp      =   7.5;          %Sampling time (min)
    %   Parameters of the exposures
    %   ---------------------------------------------------------------
    %   Declare here names for the color channels, in proper order.  
    %       All other entries in ip.bkmd.exp must be in the same order.
    ip.bkmd.exp.Channel    =   {'BFP', 'CFP', 'YFP', 'RFP', 'CY5'};  %Name(s) of Channel(s)
    %   Declare the Filters for each Channel
    ip.bkmd.exp.Filter     =   {'Filter_DAPI', 'Filter_CFP', ...
        'Filter_YFP', 'Filter_Cherry2', 'Filter_Cy5'};
    %   Declare the target Fluorophores for each Channel
    ip.bkmd.exp.FPhore     =   {'mTagBFP2', 'mCer3', 'YPet', 'mCherry', 'mCherry'};
    %       -> See iman_naming for valid Filter and Fluorophore Names    
    %   Declare the Channel names associate with each FRET pair
    ip.bkmd.exp.FRET       =   []; % {[1,2]}
    ip.bkmd.exp.Light      =   'SOLA';           %Light source name
    ip.bkmd.exp.Exposure   =   [200, 100, 100, 200, 600];  %Exposure times (ms)
    ip.bkmd.exp.ExVolt     =   [20,  20,  17,  20,  100];     %Relative voltage (0-100)
    %   IF using the SPECTRAX light source, also define the following:
%     ip.bkmd.exp.ExLine     =   {[],[],[]};   %Excitation line
%     ip.bkmd.exp.ExWL       =   {[],[],[]};   %Wavelength

%Show final image parameters (Comment out the following to avoid printing)
% display(ip);
% fprintf('\nip.indsz = \n\n');  display(ip.indsz); 
% fprintf('\nip.bkmd = \n\n');   display(ip.bkmd);
%   ------ ------ END OF IMAGE DATA DEFINITION SECTION ------ ------


%% Define Operation Parameters
%   -------- EDIT THIS SECTION TO SETUP THE PROCESSING RUN --------
%Procedure scope - target ranges of the data to process
op.cind     = {'BFP','CFP','YFP','RFP'};    %Indices of Channels to process (names or indices)
op.xypos    = 1;     %Indices of XY positions to process
op.trng     = [];       %Start and End indices of Time to process
op.nW       = 4;        %Number of parallel workers to use
op.objbias  = true;     %Correct for objective view bias
op.unmix    = true;  	%Linearly unmix color channel cross-talk
op.fixshift = false;  	%Correct any indicated frameshifts
%   Indicate which fields of bkmd should override any other MetaData source
op.mdover   = { 'exp', 'Channel';   'exp', 'Exposure';  ...
                'exp', 'ExVolt';    'exp', 'Filter';    'exp', 'FPhore'};

%Procedure settings
%   Segmentation settings
op.seg.chan = 'CFP';       %Color channel on which to segment
op.seg.cyt  = false; 	%Segment on cytoplasmic signal
%       The following defaults are typical of MCF10As at 20x magnification
op.seg.maxD = 40;       %Maximum nuclear diameter (pixels)
op.seg.minD = 12;       %Minimum nuclear diameter (pixels)
op.seg.maxEcc = 0.9;    %Maximum eccentricity [0-1] (opposite of circularity)
op.seg.Extent = [0.65, 0.85]; %Minimum fraction of bounding box filled 
%   Extent range is [0-1], value for a cirlce is PI/4
op.seg.maxSmooth = 0.1;     %Maximum fraction of convex hull without mask 
%       It is recommended to scale imaging for nuclei (segmentation
%       targets) to be greater than 10 pixels in diameter.
op.seg.sigthresh = []; %[Optional] Minimum intensity of 'good' signal
%   If using sigthresh, pre-subtract camera baseline (default is 1.5*bkg)
op.seg.hardsnr = false;  %Typically kept to FALSE.  TRUE makes the signal 
%   threshold 'hard', enforcing cutoff of any pixels below it.
op.seg.nrode = 0;   %[Optional] Approx. number of pixels to erode nuc mask
%   NOTE: op.seg.nrode is only to improve nearby cell separation - it is
%   undone in final masks.  Use op.msk.nrode to adjust final masks.

%    Masking settings
op.msk.rt = {};
op.msk.storemasks = false;  	%Save all segmentation masks
op.msk.saverawvals = true;  %Save all raw valcube entries (pre-tracking)
op.msk.nrode = 1;   %[Optional] Approx. # pixels to erode final nuc mask
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
op.trk.minTrkLength = 3;    % Minimum number of consecutive frames required
                            %   for a track to be used in gap closing
op.trk.gaprad = 15;         % Radius (in um) for gap closing
op.trk.gapwin = 5;          % Time window (in frames) for gap closing
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

%Show final operation parameters
% display(op);
% fprintf('\nop.seg = \n\n');  display(op.seg); 
% fprintf('\nop.disp = \n\n');  display(op.disp); 
%   ------ ------ END OF PROCESSING RUN SETUP SECTION ------ ------
%   Do NOT edit the following parameter transfers
op.cname    = ip.bkmd.exp.Channel;      %Copy channel names to op
op.msk.fret = ip.bkmd.exp.FRET;     %Copy FRET channel names to op


%% Run image analysis procedure
%Pre-Run validation
[ip, op] = celltracer_v2_partner('validate', ip, op);

%Main procedure
[dao, GMD, dmx] = iman_celltracer(ip, op);


%Cytoplasmic filtering test
fprintf('\nBeginning test with mask filtering.\n');
op.seg.chan = {'CFP','YFP'};
op.msk.nfilt = true; op.msk.cfilt = true; 
[ip, op] = celltracer_v2_partner('validate', ip, op);
[dao, GMD, dmx] = iman_celltracer(ip, op);

end





