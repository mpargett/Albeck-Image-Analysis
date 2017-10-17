%CELLTRACER_V2_TEST_MPTIFF
%   Tests the basic functions of the the Celltracer image processing
%   pipeline, using a MultiPage Tiff for source images.


function [ip, op] = celltracer_v2_test_mptiff()
%Initialize Celltracer Parameter Structures
[ip, op] = celltracer_v2_partner('initialize');


%% Define Imaging Parameters
%   --------- EDIT THIS SECTION TO DEFINE THE IMAGE DATA ---------
%Data locations
%   Target imaging data file name (include path)
ip.fname = ['L:\imageData\ProcessingValidation\',...
    'MultiPage Tiff Samples\nd_crop.tif'];  
%   Target save name for processed data files (include path)
ip.sname = 'celltracer_v2_test_mptiff';  

%Data - basic statistics
%   Provide the expected size of imaging data
ip.indsz.xy = 2;   %Number of XY positions
ip.indsz.z  = 1;   %Number of Z positions
ip.indsz.t  = 10;   %Number of Time points
ip.indsz.c  = 2;   %Number of color Channels

%Data - detailed definitions
%   Provide the baseline image intensity value (from Camera)
ip.bval = 100;
%   Initialize the properly sized bkg structure
ip.bkg(ip.indsz.xy) = ip.bkg;   
[ip.bkg.fix]    = deal(false);  %Default fix is false
[ip.bkg.dyn]    = deal(true);   %Default dyn is true
[ip.bkg.bkonly] = deal(false);  %Default bkonly is false
%   Provide background regions or values for each XY position
ip.bkg(1).reg =  [268, 45; 435, 136]; 
ip.bkg(2).reg =  [40, 711; 110, 785];

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
    ip.bkmd.cam.tsamp      =   6;          %Sampling time (min)
    %   Parameters of the exposures
    %   ---------------------------------------------------------------
    %   Declare here names for the color channels, in proper order.  
    %       All other entries in ip.bkmd.exp must be in the same order.
    ip.bkmd.exp.Channel    =   {'CFP', 'YFP'};  %Name(s) of Channel(s)
    %   Declare the Filters for each Channel
    ip.bkmd.exp.Filter     =   {'Filter_CFP', 'Filter_YFP'};
    %   Declare the target Fluorophores for each Channel
    ip.bkmd.exp.FPhore     =   {'mCer3', 'YPet'};
    %       -> See iman_naming for valid Filter and Fluorophore Names    
    %   Declare the Channel names associate with each FRET pair
    ip.bkmd.exp.FRET       =   struct('EKAR', {{'CFP','YFP'}}); % {[1,2]}
    ip.bkmd.exp.Light      =   'SOLA';           %Light source name
    ip.bkmd.exp.Exposure   =   [400, 400];  %Exposure times (ms)
    ip.bkmd.exp.ExVolt     =   [50, 25];     %Relative voltage (0-100)
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
op.cind     = 1:ip.indsz.c;    %Indices of Channels to process (names or indices)
op.xypos    = 1:2;     %Indices of XY positions to process
op.trng     = [1];       %Start and End indices of Time to process
op.nW       = 10;        %Number of parallel workers to use
op.objbias  = true;     %Correct for objective view bias
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
op.seg.maxD = 40;       %Maximum nuclear diameter (pixels)
op.seg.minD = 12;       %Minimum nuclear diameter (pixels)
op.seg.maxEcc = 0.9;    %Maximum eccentricity [0-1] (opposite of circularity)
op.seg.Extent = [0.65, 0.85]; %Minimum fraction of bounding box filled 
%   Extent range is [0-1], value for a cirlce is PI/4

%       It is recommended to scale imaging for nuclei (segmentation
%       targets) to be greater than 10 pixels in diameter.
op.seg.sigthresh = []; %[Optional] Minimum intensity of 'good' signal
%   If you use sigthresh, pre-subtract any camera baseline
op.seg.hardsnr = false;  %Typically kept to FALSE.  TRUE makes the signal 
%   threshold 'hard', enforcing cutoff of any pixels below it.

%    Masking settings
op.msk.rt = {};
op.msk.storemasks = false;  	%Save all segmentation masks
op.msk.saverawvals = true;  %Save all raw valcube entries (pre-tracking)

% Tracking settings (for utrack)
op.trk.movrad = 25;     % Radius (in um) to consider for cell movement
op.trk.linkwin = 75;    % Time window (in min) to consider for tracking

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

% --- Main procedure, for different export varieties --- 
%   Image kept as Multi-Channel (channels as pages)
fprintf('\nBeginning for MultiPage TIFF with MultiChannel (as pages).\n');
[dao, GMD, dmx] = iman_celltracer(ip, op);
%   Mono image for each Channel
fprintf('\nBeginning for MultiPage TIFF with Mono Channels.\n');
ip.fname = ['L:imageData\ProcessingValidation\',...
    'MultiPage Tiff Samples\nd_crop_monoc1.tif'];  
[dao, GMD, dmx] = iman_celltracer(ip, op);
%   Channels Merged (to RGB) - Not supported or recommended
% ip.fname = ['L:\albeck\imageData\ProcessingValidation\',...
%     'MultiPage Tiff Samples\nd_crop_merge.tif'];



end





