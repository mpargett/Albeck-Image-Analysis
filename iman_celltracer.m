%IMAN_CELLTRACER
%   Main image analysis function for extracting single cell time-series
%   traces from sequential images.  The procedure extracts images from file
%   storage (ND2 files or TIFF images), segments each image to identify
%   cells, and generates masks to isolate pixels in both the nucleus and
%   cytoplasm of each identified cell.  The mean value of each mask is
%   stored as the intensity for that cell per color channel and time point.
%   Multiple frame locations (XY positions) can be processed in a single
%   run (serially).
%   
%   [dao, GMD, dmx] = iman_celltracer(p, op)
%
%   It is recommended to call this procedure using the celltracer script,
%   modified to suit the target data.  See celltracer_v2_template.m
%
%   The input parameter structure, p, should describe the data to be used.
%   It is best to fill out fields for all of the data available in a file
%   (or TIFF set). Fields in op (see below) can select subsets of the data
%   to process.  p must provide the following fields: 
%
%   fname - [str] Path and filename to image(s).
%   sname - [str] Base filename for saving processed data.
%   indsz - [struct] Gives expected sizes of image data dimensions.
%       xy - [num] Number of frame locations (XY positions).
%       z - [num]  Number of vertical location (Z positions).
%       t - [num]  Number of time points.
%       c - [num]  Number of color Channels.
%   bkg - [struct] Directs how to collect image background intensity.
%     Prepare an element [bkg(2) etc.] for each XY position.
%       reg - [num] 2x2 Array of pixel positions defining a box containing
%           background only, specified by upper left (ul) corner and lower
%           right (lr) corner [x_ul, y_ul; x_lr, y_lr]. Alternately, can
%           provide fixed intensity values as 1xN array.
%       fix - [log] TRUE if .reg contains fixed intensity values instead.
%       dyn - [log] TRUE if .reg is valid for all time points. If FALSE,
%           will only sample background at time 1.
%       altxy - [num] XY position of an alternate frame from which to
%           collect background.
%       bkonly - [log] TRUE if this XY is only for background (no cells).
%   xyshift - [struct] Directs where and how to correct for frame shifts.
%     Shifts may be specified by providing the size of displacement
%     (manually observed) or registered automatically.
%       frame - [num] Time value(s) of first point after each shift.
%       dx - [num] IF specifying shifts, distance in X (i.e. x2 - x1). Must
%           be a vector the same size as .frame. Leave empty if registering
%           shifts automatically. 
%       dy - [num] Y distances, as with dx.
%       trackwell - [num] IF registering shifts, XY position to use. Leave
%           empty to register each XY position separately (recommended).
%       trackchan - [num,str] IF registering shifts, color Channel to use.
%           May be a channel index or name.
%       badtime - [num] Vector of any time points with empty frames.
%       shifttype - [string] Specify the type of shift to consider
%           ('translation', 'rigid').
%   bkmd - [struct] Defines the MetaData for the images, as a backup to
%     that extracted from storage.
%       obj - [struct] Parameters of the objective
%           Desc - [str]     Descriptive string
%           Mag - [num]      Magnification value
%           WkDist - [num]   Working distance (mm)
%           RefIndex - [num] Index of Refraction
%           NA - [num]       Numerical Aperture
%       cam - [struct] Parameters of the camera
%           Desc - [str]     Descriptive string (model)
%           PixSizeX - [num] Horizontal size per pixel (um)
%           PixSizeY - [num] Vertical size per pixel (um)
%           PixNumX - [num]  Horizontal number of pixels
%           PixNumY - [num]  Vertical number of pixels
%           BinSizeX - [num] Number of pixels per bin (horiz)
%           BinSizeY - [num] Number of pixels per bin (vert)
%           tsamp - [num]    Time between samples (min)
%       exp - [struct] Parameters of the exposure
%           Channel - [cell] Name(s) of each color Channel
%           Filter - [cell]  Name(s) of Filter used per Channel
%           FPhore - [cell]  Name(s) of Fluorophore in each Channel
%           FRET - [struct]  Give names of Channels involved in FRET (donor
%               first, e.g.: .FRET.EKAR = {'CFP', 'YFP'};
%           Light - [str]    Descriptive name of ligtht source
%           Exposure - [num] Exposure time, per Channel (ms)
%           ExVolt - [num]   Relative voltage, per Channel (0-100)
%           ExLine - [cell]  IF a multiline source, line indices used.
%           ExWL - [cell]    IF a multiline source, line WaveLengths.
%   bval - [num] Baseline intensity value from camera used.
%
%   The operation parameter structure, op, must provide:
%
%   cname - [cell]      Channel names used (typ. copied frop p.bkmd)
%   cind - [num,cell]   Channel indices to use (by index or names)
%   xypos - [num]       XY positions to use
%   trng - [num]        Time range to use [start, end]
%   nW - [num]          Number of parallel workers to use
%   objbias - [log]     TRUE to correct objective view bias
%   unmix - [log]       TRUE to run unmixing
%   fixshift - [log]    TRUE to run frameshift correction
%   mdover - [cell]     List of bkmd fields to overwrite onto final
%       MetaData, Nx2 cell with 1st field in 1st column etc.  For example,
%       op.mdover = {'exp', 'Filter'; 'cam', 'tsamp'};
%   seg - [struct]      Defines segmentation procedure options.
%       chan - [num,str]    Index/Name of Channel to use for segmentation.
%       cyt - [log]         TRUE if need cytoplasmic segmentation.
%       maxD - [num]        Max expected diameter of cell nuclei.
%       minD - [num]        Min expected diameter of cell nuclei.
%       minF - [num]        Min expected circularity of nuclei.
%   msk - [struct]      Defines masking procedure options.
%       fret - [struct]      FRET Channels (typ. copied from p.bkmd)
%       freti - [num]        FRET indices (typ. autofilled in process)
%       rt - [cell]          Channel ratios to calculate pre-averaging, as
%           either indices or names.
%       storemasks - [log]   TRUE to store nuclear/cytoplasm masks.
%       storerawvals - [log] TRUE to store all unaligned segmented data.
%   disp - [struct]     Defines what to display during processing.
%       meta - [log]        TRUE for MetaData as extracted and finalized.
%       samples - [log]     TRUE for sample images mid-process.
%       shifts - [log]      TRUE for examples of shifts and correction.
%       warnings - [log]    TRUE to print warnings.
%   pth - [struct]      Defines paths (filled by iman_assertjavapaths)
%       jc - [cell]         Java paths (e.g. for Bioformats).
%       ml - [cell]         MATLAB paths.
%
%   Version 2.0

%   M. Pargett, Albeck Lab, UC Davis, MCB.  2016


%FIXME still needs to accommodate the Mixing Calibration Table
%   In iman_unmix.

function [dao, GMD, dmx] = iman_celltracer(p, op, varargin)
%Version checking provision
if strcmpi(p,'version')
    dao = 'v2.0'; 
    %Return structure detailing compatible versions of required code
    GMD = struct('imageaccess',1, 'getframe',1, 'getmeta',1, ...
        'cellid',2, 'cellmask',2, 'xyshift',2, 'utrack_call',1, ...
        'trackcoords',1, 'aligntrackdata',1, 'unmix',2, 'naming',1, ...
        'refine',1);
    return; 
end

%Get MATLAB version for syntax selection
vv = version;  vv = str2double(regexp(vv, '^\d*\.\d*', 'match'));

%% INPUT PARSING AND CHECKING  -------------------------------------------
%IF warnings are prohibited, disable
if ~op.disp.warnings
    warning('off', 'MATLAB:Java:DuplicateClass');
    warning('off', 'IMAN:Celltrace');
end

%Evaluate scale of processing
nxy = numel(op.xypos);      %Number of XY positions to be run
nc = numel(op.cind);        %Number of channels to be used
sz = 1;                     %Z-STACKS NOT DEFINED FOR THIS PROCEDURE
%Check number of time points requested
if      isempty(op.trng);  allt = true;    op.trng = 1;
else    allt = false;   nT = op.trng(end) - op.trng(1) + 1;
end

%Check on background information provided
nbk   = numel(p.bkg);
if nbk > nxy; p.bkg = p.bkg(1:nxy);
    warning('IMAN:Celltrace', ['Extra background information found: ',...
        'numel(bkg) > numel(xypos). Removing last ',num2str(nbk - nxy),...
        ' element(s) from bkg.']);
elseif nbk < nxy; error(['Fewer backgrounds provided than XY positions.',...
        ' The background structure (bkg) must contain an element for ',...
        'each XY position to be used.']);
end

%Establish processing order (to accommodate alernate background refs)
    %Set indices to re-order run when some XYs depend on others' background
use_altXY = ~cellfun('isempty',{p.bkg.altxy});  xyrev(op.xypos) = 1:nxy;
altXYi = NaN(1,nxy); altXYi(use_altXY) = xyrev([p.bkg(use_altXY).altxy]);
%Check validity of re-ordering structure
    bkvalid = use_altXY;
    bkvalid(~use_altXY) = arrayfun(@(x)( x.fix && numel(x.reg) == nc )...
        || ( ~x.fix && isequal(size(x.reg),[2,2]) ),   p.bkg(~use_altXY));
    unaltXYi = unique(altXYi(use_altXY));    %Check for mislabeled alts
    bkvalid( unaltXYi(use_altXY(unaltXYi)) ) = false;
    if ~all(bkvalid);
        error(['An XY indicated as an alternate background is not',...
            ' itself properly defined.  Check background structure (bkg)',...     
            ' for XY(s) ', num2str(op.xypos(~bkvalid)),'.']);
    end

%Establish parallel processing environment (syntax depends on version)
if op.nW > 1 && (  vv <  8.2 && matlabpool('size') == 0  || ...
                vv >= 8.2 && isempty(gcp('nocreate'))    ) %#ok<*DPOOL>
    %Start the local parallel pool
    if vv < 8.2;  matlabpool('local', op.nW);  
    else parpool('local', op.nW, 'IdleTimeout', Inf);  end
    %Ensure availability of proper java classes from BioFormats
    if iscell(op.pth.jc); runtext = ['javaclasspath({''',op.pth.jc{1}];
        for s = 2:numel(op.pth.jc); runtext = [runtext,''',''',...
                op.pth.jc{s}]; end;  runtext = [runtext,'''});']; %#ok<AGROW>
    else runtext = ['javaclasspath ',op.pth.jc,';'];
    end;     pctRunOnAll(runtext);
end
runid = ceil(rand*1e4);                     %Assign random run ID number


%% DATA PREPARATION  -----------------------------------------------------

%Establish Data Access (e.g. BioFormats Reader for ND2 files)
msg = sprintf(['Loading Data Access Object (DAO). ',...
    'Start time is: ', datestr(now, 'mm/dd HH:MM PM')]);    disp(msg);
dao = iman_imageaccess(p.fname, p.indsz);   isTS = dao.info.ThreadSafe; 
if ~isTS; warning('off', 'MATLAB:Java:ConvertFromOpaque'); end
msg = sprintf(['DAO loaded. ', datestr(now, 'mm/dd HH:MM PM')]);
disp(msg); msg = [];   %Clear msg to prevent partial erasing on next update

%Collect image MetaData (and assert completeness)
GMD = iman_getmeta(dao, p.bkmd, op.mdover, op.disp.meta);


%% PREPROCESSING PREPARATION  --------------------------------------------
%Prepare linear unmixing (channel cross-talk correction) if applicable
if op.unmix; [~, dmx, mxs, op.msk.freti] = iman_unmix([],GMD);  %#ok<ASGLU>
else dmx = [];  mxs = []; end    %#ok<NASGU>

%Save MetaData and XTalk information
%   Verify target save directory
sd = regexp(p.sname, '^(.*)[\\/][^\\/]*$', 'tokens');  
if~isempty(sd); sd = sd{1}{1}; if ~exist(sd,'dir'); mkdir(sd); end; end
%   Save metadata etc.
save([p.sname, '_Global.mat'], 'GMD', 'dmx', 'mxs', 'p', 'op');
%   Short circuit if no XY positions specified (i.e. run was for meta-data)
%       Or if abort condition was satisfied in the DAO
if nxy == 0 || dao.abort;  dao.abort = false; return; end

%Establish data output sequence (to be stored for reference)
vcorder = iman_cellmask([GMD.cam.PixNumY, GMD.cam.PixNumX, nc], ...
    [], op, []); %#ok<NASGU>
    
%Evaluate geometric objective-based intensity correction (per metadata)
im = zeros(GMD.cam.PixNumY, GMD.cam.PixNumX, nc);
[~, e_inv] = iman_refine(im, GMD, p.bval, op.objbias);

    % - If displaying samples, setup parameters - 
    if op.disp.samples||op.disp.shifts; dsz = subf_sampdisp('setup',GMD); end


%% ----- Sequentially process each XY position desired -----
bki = cell(nxy,1);  hi = 1:nxy;     %Initialize background variables
t_tot = tic;  t_perxy = [];  emsg = [];  errid = 0;  nrun = 0; ctrans = [];
for h = [hi(~use_altXY), hi(use_altXY)] %First run XYs with own background
    t_xy = tic;         nrun = nrun + 1;
    try  %Protect other xy runs from failures
        %Ensure memory is not occupied from previous run
        clear valcube masks movieInfo tracksFinal coord
        
        %Get number of time points for this XY
        %   Allows restriction of time to a specified range by trng
        if allt
            if isfield(dao.info.imax,'t'); nT = dao.info.imax.t;
            elseif isfield(dao.info, 'mpdim'); nT = dao.info.mpdim.t;
            end; 	op.trng(2) = nT;
        end
        
        %Get non-dynamic background intensities, if applicable
        if ~use_altXY(h)
            if p.bkg(h).fix;   bki{h} = p.bkg(h).reg - p.bval;  %IF fixed
            elseif ~p.bkg(h).dyn   %IF static and region provided for image 1
                %Load first image for this XY
                im = zeros(GMD.cam.PixNumY, GMD.cam.PixNumX, nc);
                for sc = 1:nc;
                    im(:,:,sc) = iman_getframe(dao, ...
                        [op.trng(1), op.cind(sc), op.xypos(h), sz]);
                end
                %Refine image
                im = iman_refine(im, GMD, p.bval, e_inv);
                %Get background valued from first image region
                bki{h} = mean(mean( ...
                    im( p.bkg(h).reg(3):p.bkg(h).reg(4), ...
                    p.bkg(h).reg(1):p.bkg(h).reg(2), : ), 1 ), 2);
            end
        elseif ~p.bkg(h).dyn || p.bkg(h).fix; bki{h} = bki{altXYi(h)};
        end
        
        %Determine slicing structure for parallel processing
        ipw = ceil( nT./op.nW );               %Approx. images per Worker
        nper = ipw*ones(1,op.nW);  dif = op.nW*ipw - nT - 1; %Fix overestimate
        nper(end-dif:end) = nper(end-dif:end) - 1;  %Minus 1 excess per Worker
        
        %% Segment and mask images in parallel
            % - Display status message -
            t_el = toc(t_tot)/60; t_rem = t_perxy*(nxy-h+1) - toc(t_xy)/60;
            [msg, emsg] = subf_update_timer(t_el, t_rem, msg, emsg);
            msg = sprintf(['Segmenting cells for XY ', ...
                '%d (%d/%d).'], op.xypos(h), nrun, nxy);    disp(msg);
        
        	%Pre-allocate loop variables as needed
            movieInfo = cell(1,op.nW);  bkis = cell(op.nW,1);
            valcube = cell(1,op.nW);    masks = cell(1,op.nW);  
            dbk = masks;  sbk = dbk;  %Backups for display/debug
        parfor ps = 1:op.nW   %Parallel Loop 
            %Load Data Access on each Worker (may need to reinitialize)
            if isTS;    daop = dao;          
            else daop = iman_imageaccess(p.fname, p.indsz);  %#ok<PFBNS>
            end 

            %Pre-allocate and run serial sub-loop
            im = zeros(GMD.cam.PixNumY, GMD.cam.PixNumX, nc); %#ok<PFBNS>
            for s = 1:nper(ps)
                %Get current time index (as start + already completed)
                st = op.trng(1)-1 + sum( nper(1:(ps-1)) ) + s;  %#ok<PFBNS> 
                %Load Image
                for sc = 1:nc
                    im(:,:,sc) = iman_getframe(daop, ...
                        [st, op.cind(sc), op.xypos(h), sz]);
                end
                
                %Refine Image (storing examples of first image)
                    if ps == 1 && s == 1; dbk{ps}.orim = im; end
                im = iman_refine(im, GMD, p.bval, e_inv);
                    
                %Get dynamic background image values, if applicable
                if p.bkg(h).dyn && ~p.bkg(h).fix
                    if use_altXY(h); bkit = bki{altXYi(h)}(st-op.trng(1)+1,:); %#ok<PFBNS>
                    else    bkit = mean(mean( ...
                            im( p.bkg(h).reg(3):p.bkg(h).reg(4), ...
                            p.bkg(h).reg(1):p.bkg(h).reg(2), : ), 1 ), 2);
                        bkis{ps}(s,:) = bkit;    %Store background used
                        %IF for background only, skip cell ID
                        if p.bkg(h).bkonly; continue; end
                    end
                else bkit = bki{h};
                end
                
                %Remove background from image
                im = bsxfun( @minus, im, reshape(bkit, 1, 1, numel(bkit)) );
                %Perform linear channel unmixing (optional)
                if op.unmix; im = iman_unmix(im,[],[],dmx,op.msk.freti); end
                
                %Perform Cell ID by segmentation
                [movieInfo{ps}{1,s}, nmask] = iman_cellid(im, op, bkit);
                
                %Perform Cell signal measurement, by masking
                [valcube{ps}{s}, masks{ps}{s}] = iman_cellmask(im, ...
                    movieInfo{ps}{1,s}, op, nmask);
                
                    % - Backup for display/debugging, when needed - 
                    if ~isempty(dbk{ps}); dbk{ps}.mz = movieInfo{ps}{s};
                        dbk{ps}.vc = valcube{ps}{s};  dbk{ps}.rfim = im;
                        dbk{ps}.msk = masks{ps}{s};
                    end
                    if op.disp.shifts && ~isempty(p.xyshift) && ...
                            any(st == [p.xyshift.frame-1, p.xyshift.frame]) 
                        sx = length(sbk{ps}) + 1; sbk{ps}(sx).frame = st;  
                        sbk{ps}(sx).im = single(im);
                    end
            end
            %   IF not storing masks, empty the variable
            if ~op.msk.storemasks;  masks{ps} = {}; end
        end     %END PARALLEL LOOP
        
        if p.bkg(h).dyn    %IF using dynamic background, store values
            if use_altXY(h); bki{h} = bki{altXYi(h)};
            else bki{h} = cat(1, bkis{:}); end
        end
        %IF this well is for background only, skip further evaluation
        if p.bkg(h).bkonly; continue; end
        
            % - Display sample images, if applicable - 
            if op.disp.samples; subf_sampdisp('cellid', runid, nc, dbk, ...
                    GMD, op, dsz); end
        
        %Collect parallelized outputs
        movieInfo = cat(2, movieInfo{:});     %Cat parallel layer
        movieInfo = cat(2, movieInfo{:});     %Cat serial layer
        
        valcube = cat(2, valcube{:});  %valcube = cat(2, valcube{:});
        masks = cat(2, masks{:});      masks = cat(2, masks{:}); %#ok<NASGU>
        
        
        %% Perform Cell Tracking
            % - Update timing and display -
            t_el = toc(t_tot)/60; t_rem = t_perxy*(nxy-h+1) - toc(t_xy)/60;
            [msg, emsg] = subf_update_timer(t_el, t_rem, msg, emsg);
            msg = sprintf(['Building Tracks for XY ', ...
                '%d (%d/%d).'], op.xypos(h), nrun, nxy);    disp(msg);
            
        %Correct for frameshifts in images, IF applicable
        if  op.fixshift
                % - Update display if fixing frameshifts
                msg2 = sprintf('\b - Fixing frame shift events.'); disp(msg2);
            if isempty(p.xyshift.trackwell) || isempty(ctrans)
                p.xyshift.curwell = op.xypos(h);
                [movieInfo, ctrans] = iman_xyshift(movieInfo, p.xyshift, 1, dao);
            else
                movieInfo = iman_xyshift(movieInfo, ctrans, 1, dao);
            end
                % - Update timing and display -
                t_el = toc(t_tot)/60; t_rem = t_perxy*(nxy-h+1)-toc(t_xy)/60;
                [~, emsg] = subf_update_timer(t_el,t_rem,[msg,msg2(2:end)],emsg);
                disp(msg);
        end
        
        %Prepare for uTrack - movieInfo fields must include a 2nd column of
        %   zeros, for Std Deviation of position and amplitude 
        for s = 1:numel(movieInfo)
            movieInfo(s) = structfun(@(x)[x, zeros(size(x,1),1)], ...
                movieInfo(s), 'UniformOutput', false);
        end
        %Call uTrack (uses trackCloseGapsKalmanSparse)
        tracksFinal = iman_utrack_call(movieInfo, GMD.cam);
        
        %Perform track cleanup to deliver xC and yC matrices
        %   Provides sliced matrices for contouring (ensuring integers)
        coord = round(iman_trackcoords(tracksFinal));
        %   If only one time point, coord will be empty.  Fill.
        if numel(movieInfo) == 1;
            coord = cat(3,movieInfo.xCoord(:,1), movieInfo.yCoord(:,1));
        end
        %If coords is not full length (no tracks are full), pad with NaN
        coord(:,(end+1):nT,:) = NaN;
        
        %Reverse the XY shift, if used
        if  op.fixshift;  coord = iman_xyshift(coord, ctrans, -1);  end
        
        %Align data (from masked regions) to tracked coordinates
        if op.msk.saverawvals; vcraw = valcube; else vcraw = []; end %#ok<NASGU>
        valcube = iman_aligntrackdata(coord, valcube); %#ok<NASGU>
        
            % - Display the adjustments for frame shift, as applicable -
            if op.disp.shifts && ~isempty(p.xyshift) && ~isempty(p.xyshift.frame)
                subf_shift_display('shift', movieInfo, coord, p.xyshift, ...
                    sbk, op, ctrans, dsz, runid);
            end

        
        %% Finalize and Save
        %Save Processed Data to desired location ('0' as placeholder)
        dh = '0'; xystr = ['_xy', dh(op.xypos(h)<10), num2str(op.xypos(h))];
        bkg = bki{h}; %#ok<NASGU>
        save([p.sname, xystr, '.mat'], 'valcube', 'vcorder', 'vcraw', ...
            'bkg', 'masks', 'ctrans');
        
        %Display status message
        t_el = toc(t_tot)/60;  t_rem = t_perxy*(nxy-h);
        [msg, emsg] = subf_update_timer(t_el, t_rem, msg, emsg);
        msg = sprintf(['Processing completed for XY ', ...
            '%d (%d/%d). - %s'], op.xypos(h), nrun, nxy, datestr(now));
        disp(msg);  msg = []; emsg = []; %Clears message to keep this note.
        
    catch ME  %OPERATIONS TO PERFORM IN CASE OF ERROR
        errid = errid + 1;  %Set the new error ID number
        %Display an error message
        dmsg = sprintf(repmat('\b', 1, numel(msg)+numel(emsg)+2)); 
        msg = sprintf(['Processing FAILED on XY ', ...
            '%d (%d/%d).'], op.xypos(h), h, nxy);  disp(dmsg);  disp(msg);
        disp(ME.message);  emsg = []; msg = [];
        %Produce an Error Report and proceed to the next XY position
        MElog{errid}(1:2) = { msg, [ME.identifier, ' - ', ME.message] }; %#ok<AGROW>
        for se = 1:length(ME.stack)
            MElog{errid}(2 + se) =   {sprintf('%s  - Line %d', ...
                ME.stack(se).file, ME.stack(se).line) };  
        end
        %Append a structure with relevant info to the end of MElog
        MElog{errid}{end} = struct('h',h); %#ok<AGROW>
        continue
    end
    
        % - Update timing counter -
        t_el = toc(t_tot)/60;   t_perxy = t_el/nrun;
end


%% End of procedure wrap-up
%Close the parallel enviroment and return (syntax depends on version)
if vv < 8.2;    if matlabpool('size') > 0;    matlabpool close;   end
else            delete(gcp('nocreate'));
end

%Error recording: Generate logfile for any errors encountered
if errid > 0
    %Create an error log file
    ct = now;    erfname = ['CellTrace_ErrorLog_',datestr(ct,30),'.txt'];
    erfid = fopen(erfname, 'w+');  %Open log file for writing
    fprintf(erfid, '%s\n%s\n\n', ['Error Log File for iman_celltrace_main,',...
        ' run ', datestr(ct, 31), ' for file: '], p.fname);
    %Write to file all error information in series
    for se = 1:errid
        fprintf(erfid, '\nError %d\n', se);         %New title line
        fprintf(erfid, '\nIn XY Position %d\n', ...
            op.xypos(MElog{se}{end}.h));          %XY Pos
        for sse = 1 : (length(MElog{se}) - 1)    %Print error messages
            fprintf(erfid, '%s\n', MElog{se}{sse});
        end
    end
    fclose(erfid);  %Close error file to finalize
end

end




%% SUBFUNCTIONS

% --- Timer display update function ---
function [msg, emsg] = subf_update_timer(t_el, t_rem, msg, emsg)
fprintf(repmat('\b', 1, numel(msg)  + ~isempty(msg) + ...
    numel(emsg) + ~isempty(emsg)  ));
msg = [];
%IF t_rem not yet defined, skip time estimation
if isempty(t_rem);      return;     end
emsg = sprintf(['Time elapsed: %dh %dm.  Estimate time remaining: ',...
    '%dh %dm.'], floor(t_el/60),  round(mod(t_el,60)), ...
    floor(t_rem/60), round(mod(t_rem,60))  );
disp(emsg);
end


% --- Sample display management ---
function out = subf_sampdisp(in, varargin)
%   Called in different contexts to reduce main code bloat
switch in
    case 'setup'
        %Use: dsz = subf_sampdisp('setup', GMD);
        GMD = varargin{1};
        out.thet = 0:(2*pi/119):(2*pi);  %For displaying ~nuclear location
        %Set sample display parameters  (preserves aspect ratio)
        out.x = min(640,GMD.cam.PixNumX); out.y = min(540,GMD.cam.PixNumY);
        out.c = min([out.x./GMD.cam.PixNumX, out.y./GMD.cam.PixNumY]);
        out.x = GMD.cam.PixNumX*out.c;    out.y = GMD.cam.PixNumY*out.c;
        
    case 'cellid'
        %Use: subf_sampdisp('cellid', runid, nChan, dbk, GMD, op, dsz);
        [runid, nChan, dbk, GMD, op, dsz] = deal(varargin{:});
        dbk = [dbk{:}];  %Reduce backup image cell array
        
        %SHOW REFINED IMAGES
        figure(runid); clf reset; %Figure with unique ID for this run
        set(runid, 'Position', [80, 50, 700, 920]);
        for ds = 1:nChan
            subplot(nChan,2,1+(ds-1)*2);
            set(gca,'Position', get(gca,'Position') + ...
                [-0.08, -0.02*ds, 0.1, 0.06] );
            imshow(dbk.orim(:,:,ds), [],...
                'Border', 'tight', 'InitialMagnification', 'fit');
            lh = ylabel(GMD.exp.Channel{ds});
            lhp = get(lh,'Position'); ax = axis;
            set(lh, 'Position', [-0.025*ax(4), lhp(2:end)]);
            if ds == 1; title('Original'); end
            subplot(nChan,2,2+(ds-1)*2);
            set(gca,'Position', get(gca,'Position') + ...
                [-0.05, -0.02*ds, 0.1, 0.06] );
            imshow(imresize(dbk.rfim(:,:,ds), [dsz.y,dsz.x]/2), [], ...
                'Border', 'tight', 'InitialMagnification', 'fit');
            if ds == 1; title('Refined'); end
        end
        %SHOW CELL ID RESULTS
        figure(runid+1); clf reset; %Figure with unique ID for this run
        set(runid+1, 'Position', [80, 50, 710, 600]);
        %Show Image
        imshow( imresize(dbk.rfim(:,:,op.seg.chan), [dsz.y,dsz.x]), [], ...
            'Border', 'tight', 'InitialMagnification', 'fit');  hold on;
        %Plot centers
        scatter(dbk.mz.xCoord.*dsz.c, ...
            dbk.mz.yCoord.*dsz.c, 75, 'r+');
        %Plot circles representing ~nuclear area
        for ds = 1:length(dbk.mz.xCoord)
            ampy = sqrt( dbk.mz.are(ds).*dsz.c.*dsz.c ./ pi);
            plot(dbk.mz.xCoord(ds).*dsz.c + ampy.*cos(dsz.thet), ...
                dbk.mz.yCoord(ds).*dsz.c + ampy.*sin(dsz.thet), 'b:');
        end
        %Plot reference cell areas
        plot( (op.seg.maxD/2 + 0.5*op.seg.maxD.*cos(dsz.thet)).*dsz.c, ...
            (op.seg.maxD/2 + 0.5*op.seg.maxD.*sin(dsz.thet)).*dsz.c, 'r:' );
        plot( (op.seg.maxD/2 + 0.5.*op.seg.minD.*cos(dsz.thet)).*dsz.c, ...
            (op.seg.maxD/2 + 0.5.*op.seg.minD.*sin(dsz.thet)).*dsz.c, 'r:' );
        drawnow;
        
        figure(runid+2); clf reset; %Figure with unique ID for this run
        set(runid+2, 'Position', [80, 50, 710, 600]);
        
        %Get mask boundaries from stored mask images
        nmim = bwboundaries(bwunpack(dbk.msk.nuc, GMD.cam.PixNumY));
        cmim = bwboundaries(bwunpack(dbk.msk.cyt, GMD.cam.PixNumY));
        imshow( dbk.rfim(:,:,op.seg.chan), [], ...
            'Border', 'tight', 'InitialMagnification', 'fit'); hold on;
        for ds = 1:numel(nmim)
            plot(nmim{ds}(:,2),nmim{ds}(:,1), 'g-');
        end
        for ds = 1:numel(cmim)
            plot(cmim{ds}(:,2),cmim{ds}(:,1), 'b-');
        end
        drawnow;
        
    case 'shift' % --- Sample image display for shift comparisons ---
        %Use: [] = subf_shift_display('shift', movieInfo, coord, xyshift, ...
        %         sbk, op, ctrans, dsz, runid);
        [mz, coord, xyshift, sbk, op, ctrans, dsz, runid] = deal(varargin{:});
        sbk = [sbk{:}];  %Cat shift display images
        
        sfts = [sbk.frame];      tpts = sfts - ctrans(1).stime + 1;
        for ds = 1:numel(xyshift.frame)
            ds1 = xyshift.frame(ds)-1 == sfts;
            ds2 = xyshift.frame(ds) == sfts;
            figure(runid+10+ds); clf reset; %Figure with unique ID for this run
            set(runid+10+ds, 'Position', [80, 50, 1400, 600]);
            subplot('position', [0,0,0.475,0.95]); %Start with pre-shift on left
            
            %Show pre-shift image
            imshow( imresize(sbk(ds1).im(:,:,op.seg.chan), [dsz.y,dsz.x]), [], ...
                'Border', 'tight', 'InitialMagnification', 'fit');  hold on;
            %Plot pre-shift centers (post-process)
            scatter(coord(:,tpts(ds1),1).*dsz.c, ...
                coord(:,tpts(ds1),2).*dsz.c, 75, 'r+');
            %Plot post-shift centers (post-process)
            scatter(coord(:,tpts(ds2),1).*dsz.c, ...
                coord(:,tpts(ds2),2).*dsz.c, 75, 'bx');
            title(['Pre-shift image, frame ', num2str(sbk(ds1).frame),...
                '. Crosses: post, Circles: pre']);
            
            %Show also on post-shift image
            subplot('position', [0.525,0,0.475,0.95]); %Show post-shift on right
            %Show post-shift image, aligned
            imshow( imresize(iman_xyshift(sbk(ds2).im(:,:,op.seg.chan), ...
                ctrans, 1, [], sbk(ds2).frame), [dsz.y,dsz.x]), [], ...
                'Border', 'tight', 'InitialMagnification', 'fit');  hold on;
            
            %Plot pre-shift circles representing nuclear area (pre-process)
            scatter(mz(tpts(ds1)).xCoord(:,1).*dsz.c, ...
                mz(tpts(ds1)).yCoord(:,1).*dsz.c, 75, 'r+');
            %Plot post-shift circles (pre-process)
            for ds3 = 1:length(mz(tpts(ds2)).xCoord)
                ampy = sqrt( mz(tpts(ds2)).are(ds3).*dsz.c.*dsz.c ./ pi);
                plot(mz(tpts(ds2)).xCoord(ds3).*dsz.c + ampy.*cos(dsz.thet), ...
                    mz(tpts(ds2)).yCoord(ds3).*dsz.c + ampy.*sin(dsz.thet), 'b:');
            end
            title('Post-shift image, aligned, with ');
        end
        drawnow;
end

end
