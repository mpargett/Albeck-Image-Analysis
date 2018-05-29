% IMAN_CURATE_LINEAGE 
%   With iman_curate_lineage.fig, opens a graphical user interface (GUI)
%   window to curate the lineage flags (splits, merges) identified in the
%   tracking portion of the CellTracer image processing.
%
%   iman_curate_lineage(filename, ...)
%       uses the file and path specified in filename to identify and open
%       both processed data and image data for an experiment.  filename
%       must be a processed data MAT file, and must be in the same folder
%       as the _Global processed data file and any other files from the
%       same experiment to be curated.  If no file or a bad filename is
%       provided, a dialog box will open to browse to the desired file.
%       Curation will begin for the file specified, and may be continued to
%       other XYs in numerical order.  filename may be the _Global file, in
%       which case, curation will begin with XY_01.
%
%   Curation consists of viewing up to 5 cropped frames, the last of which
%   is annotated with the centroids of the two joined cells (blue and
%   yellow circles), with a magenta line joining them.  Flag the Event
%   (either a track split, or merge) by hitting the Validate button,
%   the Reject button, or the Bad Mask button.  
%
%   Additional parameters may be provided as Name, Value pairs:
%   cid     - The Curator's ID, your name or initials.  If this parameter
%               is unset, a dialog box will open for you to provide it.
%               The ID is printed in a note, saved in the processed data.
%   imth    - The threshold percentile of the image to clip when
%               displaying. Images will be scaled from 0 to the imth
%               percentile of their intensity.  Default is 99.
%   recurate -Flag to redo curation of previously curated files.  Default
%               is FALSE.

function varargout = iman_curate_lineage(varargin)
%MATLAB code for iman_curate_lineage.fig
%      IMAN_CURATE_LINEAGE, by itself, creates a new IMAN_CURATE_LINEAGE or raises the existing
%      singleton*.
%
%      H = IMAN_CURATE_LINEAGE returns the handle to a new IMAN_CURATE_LINEAGE or the handle to
%      the existing singleton*.
%
%      IMAN_CURATE_LINEAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAN_CURATE_LINEAGE.M with the given input arguments.
%
%      IMAN_CURATE_LINEAGE('Property','Value',...) creates a new IMAN_CURATE_LINEAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iman_curate_lineage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iman_curate_lineage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iman_curate_lineage

% Last Modified by GUIDE v2.5 16-May-2018 09:31:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iman_curate_lineage_OpeningFcn, ...
                   'gui_OutputFcn',  @iman_curate_lineage_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
lf = localfunctions;     %Get list of local functions
if nargin && ischar(varargin{1}) && ...
        any(strcmpi(varargin{1},cellfun(@char,lf,'Un',0)))
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



function iman_curate_lineage_OpeningFcn(hObject, eventdata, h, varargin)
%% --- Executes just before iman_curate_lineage is made visible.
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iman_curate_lineage (see VARARGIN)

%Default parameters
h.p.imth = 99;          %Image percentile threshold, for display
h.p.cid = 'Anonymous';  %Curator ID (Name)
h.p.recurate = false;   %Flag to re-curate if files already done

%Get input file from front of input arguments
if ~isempty(varargin)
    infile = varargin{1};  varargin = varargin(2:end);
else infile = [];
end

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   h.p.(lower(varargin{s})) = varargin{s+1};   end

%Get Curator ID, if not provided
if strcmpi(h.p.cid, 'anonymous')
    inputdlg('Enter your name or initials.', 'Curator ID', 1, {h.p.cid});
end

%Parse data file name
%   Collect input file name, or open a dialog box for it
if  exist(infile,'file');
    fp = regexpi(infile,'^(?<dir>.*[\\/])(?<name>[^\\/]*)$', 'names');
else  [fp.name,fp.dir] = uigetfile();
end
%   Processes file name
fp = parse_target_file(fp);
%   Store file path info
h.fp = fp;

%Get XY information
%   Load global data
gf = load([fp.dir,fp.name,'Global.mat']);
%   Load Data Access Object for imagery
h.xy.dao = iman_imageaccess(gf.p.fname);
%   Set Channel and Z indices
h.xy.chan = gf.op.seg.chan;     h.xy.z = gf.p.indsz.z;

%Parse XY list to curate
if nin > 1 && isnumeric(varargin{2});    xys = unique(varargin{2});
    [badxy, wi] = setdiff(xys, fp.xy);
    if ~isempty(badxy); xys = xys(~wi); %Remove xys not in the data
        warning('IMAN:Curate:BadXY',['XYs ', num2str(badxy), ...
            ' are not present in the datafile(s) with base name ', ...
            [fp.dir,fp.name], '. They have been removed from the ', ...
            'curation set.']);
    end
else xys = find(strcmpi(regexp(fp.apdx,'\d*$','match'), fp.xystr));
    xys = fp.xy(xys:end); %Take xys from starting file to end
end
%Output the xy list to curate, and current position in the list
h.xy.list = xys;  h.xy.n = numel(xys);    h.xy.i = 1;

%   Load new XY data and image object
[h.pd, h.xy] = get_new_xy(h.fp, h.xy, h.p.recurate);
if isempty(h.pd)    %IF no data returned, No new XYs
    %Set message text
    str = ['NO new XYs to Curate. ',...
        'Press ''Finish'' to close the curation GUI. Run with ',...
            '''recurate'' set to TRUE to review already curated XYs.'];
    set(h.progress_text, 'String', str);    %Show message
else
    %Get Event list for the current XY
    h.ev = get_new_eventlist(h.pd);
    % Update handles structure
    update_image(h);
end

% Choose default command line output for iman_curate_lineage
set( findall( hObject, '-property', 'Units' ), 'Units', 'Normalized' )
h.output = hObject;     guidata(hObject, h);    

% UIWAIT makes iman_curate_lineage wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function varargout = iman_curate_lineage_OutputFcn(hObject, eventdata, handles)
%% --- Outputs from this function are returned to the command line.
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%MP No output, but save revised data with a note that it is post curation

% Get default command line output from handles structure
varargout{1} = handles.output;


function validate_button_Callback(hObject, eventdata, h)
%% --- Executes on button press in validate_button.
% hObject    handle to validate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Do nothing with this Event (keep it)
%Iterate to next Event
h.ev.i = h.ev.i + 1;

guidata(hObject,h);   update_image(h);


function reject_button_Callback(hObject, eventdata, h)
%% --- Executes on button press in reject_button.
% hObject    handle to reject_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Flag the current Event for removal
h.ev.rej(h.ev.i) = 1;
%Iterate to next Event
h.ev.i = h.ev.i + 1;

guidata(hObject,h);     update_image(h);


function bad_mask_button_Callback(hObject, eventdata, h)
%% --- Executes on button press in bad_mask_button.
% hObject    handle to bad_mask_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if 1st, 2nd, or both masks are bad
bm = questdlg('Which track is a bad mask?', 'Bad Mask', ...
    'Blue', 'Yellow', 'Both', 'Both');

%Flag this event to indicate remove of whole tracks
%   2: 1st track (Mother) to be removed
%   3: 2nd track (Daughter) to be removed
%   4: Both tracks to be removed
switch bm
    case 'Blue';        h.ev.rej(h.ev.i) = 2;
    case 'Yellow';      h.ev.rej(h.ev.i) = 3;
    case 'Both';        h.ev.rej(h.ev.i) = 4;
end

%Iterate to next Event
h.ev.i = h.ev.i + 1;

guidata(hObject,h);     update_image(h);


function back_button_Callback(hObject, eventdata, h)
%% --- Executes on button press in back_button.
% hObject    handle to back_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%IF this is not the first Event, revert to the previous event
if h.ev.i > 1;       h.ev.i = h.ev.i - 1;
    guidata(hObject,h);     update_image(h);
end


function finish_button_Callback(hObject, eventdata, h)
%% --- Executes on button press in finish_button.
% hObject    handle to finish_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check to save unsaved data (i.e. if 'Next XY' hasn't been pushed)?
if ~isfield(h.pd, 'Curation')
    %Process rejection flags, modify data, and save
    process_curation(h);
end
%Close GUI
closereq;


function next_xy_button_Callback(hObject, eventdata, h)
%% --- Executes on button press in next_xy_button.
% hObject    handle to next_xy_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%IF the end of Events in this XY has not been reached, double-check
if h.ev.i < h.ev.n
    x = questdlg(['You have not curated all Events in this XY.  ',...
        'Do you really want to go to the next XY?'], ...
        'Skip remaining events?', 'Yes', 'No', 'No');
    if strcmpi(x,'no');   return;   end
end

%Process rejection flags, modify data, and save
process_curation(h);

%IF another XY is to be curated
if h.xy.i < h.xy.n;      
    %Iterate XY index
    h.xy.i = h.xy.i + 1;
    %Load new Data
    [h.pd, h.xy] = get_new_xy(h.fp, h.xy, h.p.recurate);
    if isempty(h.pd)    %IF no data returned, No new XYs
        %Set message text
        str = ['NO more new XYs to Curate. ',...
            'Press ''Finish'' to close the curation GUI. Run with ',...
            '''recurate'' set to TRUE to review already curated XYs.'];
        set(h.progress_text, 'String', str);    %Show message
    else
        %Get Event list for the current XY
        h.ev = get_new_eventlist(h.pd);
    end
else %IF not, display end of curation message
    %Set message text
    str = [ 'END of XYs in this dataset. ',...
        'Press ''Finish'' to close the curation GUI.'];
    %Show message
    set(h.progress_text, 'String', str);
end


function [fp] = parse_target_file(fp)
%% Parse input filename
%Parse filename into base name and appendix
tmp = regexpi(fp.name, '^(?<name>.*_)(?<apdx>Global|xy\d*)\.mat$', 'names');
%   Allocate name and appendix
fp.name = tmp.name; fp.apdx = tmp.apdx;
%Get available files
d = dir(fp.dir);    %Get all files in directory, then restrict to matching
d = d(~cellfun(@isempty, regexpi({d.name}, [fp.name,'xy'])));
%Get matching xy flags (reducing tokens cell array)
tmp = regexpi({d.name}, '_xy(\d*).mat$', 'tokens');     tmp = [tmp{:}]; 
%   Store xy strings and numeric values
fp.xy = cellfun(@str2double, tmp);      fp.xystr = [tmp{:}];
%   Sort xy maps
[fp.xy, xyi] = sort(fp.xy, 'ascend'); fp.xystr = fp.xystr(xyi);
%Setup xy string mapping function
fp.xyfun = @(x)fp.xystr{fp.xy == x};


function [hpd, hxy] = get_new_xy(hfp, hxy, hrc)
%% Open new XY
hpd = [];   %Initialize in case of no new XYs
%Allow iteration through list of XYs, if already curated
while hxy.i <= hxy.n && isempty(hpd)
    %Point to the current XY in the list
    fld = [hfp.dir,hfp.name,'xy',hfp.xyfun(hxy.list(hxy.i)),'.mat'];
    %Check if current XY has been curated already, unless re-curating
    if ~hrc && isfield(matfile(fld), 'Curation')
        hxy.i = hxy.i + 1;  %Iterate to next XY
    else
        hpd = load(fld);    %Load data file
    end
end


function hev = get_new_eventlist(hpd)
%% Assemble list of Events for a given XY dataset
%Get coordinate channel indices
xi = find(strcmpi(hpd.vcorder, 'XCoord'));   %Get X channel
yi = find(strcmpi(hpd.vcorder, 'YCoord'));   %Get Y channel
%Initialize merge/split assembly
ms = ~isnan(hpd.linfo);  nc = size(hpd.valcube,1);
%   Get list of start and end times for all tracks
dm = diff([true(nc,1),isnan(hpd.valcube(:,:,xi)),true(nc,1)],1,2);
se = nan(nc,2);
for s = find( any(dm ~= 0, 2) )'
    [~,se(s,1)] = find(dm(s,:) < 0,1,'first');
    [~,se(s,2)] = find(dm(s,:) > 0,1,'last');
end
%Merge and Split start/end vectors, Events includes both
ev = [find(ms), hpd.linfo(ms), se(ms)];
ev = ev(~any(isnan(ev),2),:);

%Assemble Event list with central coordinates
%   Event list: LInfo # (Daughter track), Mother track, Frame
nev = size(ev,1);   hev.list = ev;  %Store event list
for s = 1:nev
    %Get coordinates of Event involved tracks
    hev.coords{s} = squeeze(hpd.valcube(ev(s,1:2), :, [xi,yi]));
end
%Set Event size, counter, and rejection index
hev.n = size(hev.list,1);    hev.i = 1;
hev.rej = zeros(hev.n,1);


function update_image(h)
%% Update the display image
if h.ev.i <= h.ev.n  %IF not finished with Merge/Splits for this XY
    %Assemble and show frame
    update_event_frames(h);
    %Update display text
    str = ['Event ', num2str(h.ev.i),' of ',num2str(h.ev.n)];
    set(h.progress_text, 'String', str);  %Show display text
else    %IF finished, show a blank frame, and prompt to load next 
    imshow(zeros(101,525),'Parent',h.axes1);  %Show empty frame
    if h.xy.i <= numel(h.xy.list)
        %Define display text for end of XY
        str = [ 'END of Events in this XY. ',...
                'Press ''Next XY'' to curate additional wells, or ',... 
                'press ''Finish'' to close the curation GUI.'];
    else
        %Define display text for end of curation
        str = [ 'END of Events in this XY. ',...
                'END of XYs in this dataset. ',...
                'Press ''Finish'' to close the curation GUI.'];
    end
    %Show display text
    set(h.progress_text, 'String', str);
end


function [] = update_event_frames(h)
%% Assemble display frame for current Event
%Get scale of data
[ntrk, ntime, ~] = size(h.pd.valcube);
%Number of previous frames and frames after the Event
pfrm = 3; afrm = 1;
%Event indices (assert dependent track index from linfo index)
ei = h.ev.list(h.ev.i,:);   ei(1) = rem(ei(1)-1, ntrk) + 1;

%Initialize frame to be shown
frm = zeros(101,525);
%Get centroid coordinates
cc = reshape(h.ev.coords{h.ev.i}(:,ei(3),:), 2,2);
%Get average center coordinates
acc = round(nanmean(cc,1))';
%Get frame range around center coordinates
xy = num2cell(bsxfun(@plus,acc, -50:50),2);
%Adjust for frame edge clipping
%   Get clipping bounds
cb = cell(2,1);
cb{1} = xy{1} > 0 & xy{1} <= h.xy.dao.read.getSizeX;
cb{2} = xy{2} > 0 & xy{2} <= h.xy.dao.read.getSizeY;
cb = cellfun(@(x)find(x,1,'first') + [0, nnz(x)-1], cb, 'Un', 0);
%   Clip the image frame range (xy)
xy = cellfun(@(x,c)x(c(1):c(2)), xy, cb, 'Un', 0);
%   Assign target frame regions
tfrx = cb{1}(1):cb{1}(2);   tfry = cb{2}(1):cb{2}(2);

%Assemble frame, one time point at a time
t1 = max(-afrm, ei(3)-ntime); t2 = min(pfrm, ei(3)-1);
ctc = cell(t2-t1,1);
for s = t1:t2 
    %Get full frame for this time
    im = iman_getframe(h.xy.dao, ...
        [ei(3)-s, h.xy.chan, h.fp.xy(h.xy.i), h.xy.z]);
    %Place cropped region in frame to be shown
    frm(tfry, (pfrm-s)*106 + tfrx) = im(xy{2}, xy{1});
    %Get centroids for the frames to be shown (cropped and tiled)
    ctc{s-t1+1} = ([106*(pfrm-s); 0] - [xy{1}(1);xy{2}(1)] + ...
                   [cb{1}(1);cb{2}(1)] - 1)*ones(1,2) ...
        + reshape(h.ev.coords{h.ev.i}(:,ei(3)-s,:), 2,2)';
end

%Show image
mxi = prctile(frm(frm>0), h.p.imth);
imshow(frm, [0, mxi], 'Parent', h.axes1);  hold(h.axes1,'on');
%Overplot centroids in different colors
sty = {1, 'ko', 4; 1, 'bo', 3; 2, 'ko', 4; 2, 'yo', 3}; %Color styles
for s = 1:numel(ctc)    %PER frame
    for sc = 1:size(sty,1)  %PER style
        plot(h.axes1, ctc{s}(1,sty{sc,1}), ctc{s}(2,sty{sc,1}), ...
            sty{sc,2}, 'MarkerSize', 10, 'LineWidth', sty{sc,3});
    end
end
%Overplot the Merge/Split line on Event frame
plot(h.axes1, ctc{1-t1}(1,:), ctc{1-t1}(2,:), 'k-', 'LineWidth', 3);
plot(h.axes1, ctc{1-t1}(1,:), ctc{1-t1}(2,:), '-', ...
    'LineWidth', 2, 'Color', [1,0.4,0.8]);


function [] = process_curation(h)
%% Process curation flags and save
%Parse flags
%   Get scope
nt = size(h.pd.linfo,1);    %Number of tracks (original) 
dti = rem(h.ev.list(:,1)-1, nt)+1; 	%Dependent track index
iti = h.ev.list(:,2);               %Independent track index
%   Get track removal index
trem = [    dti(h.ev.rej == 2);     %Dependent track bad
            iti(h.ev.rej == 3);     %Independent track bad
            dti(h.ev.rej == 4);     %Both bad
            iti(h.ev.rej == 4);     ];

%Process rejections from linfo
h.pd.linfo(h.ev.list(h.ev.rej == 1,1)) = NaN;

%Process removal of whole tracks
%   Get index of tracks to keep
tki = true(size(h.pd.linfo,1),1);  tki(trem) = false;
%   Remove tracks from data
h.pd.valcube = h.pd.valcube(tki,:,:);
%   Remove tracks from linfo
h.pd.linfo = h.pd.linfo(tki,:);
%   Get index mapping for linfo adjustment
imap = nan(size(tki));  imap(tki) = 1:nnz(tki);
%   Adjust remaining linfo references
for sl = unique(h.pd.linfo(~isnan(h.pd.linfo)))'
	h.pd.linfo( h.pd.linfo == sl ) = imap(sl);
end

%Set curation note
h.pd.Curation = ['This datafile was curated for track splitting ',...
    'and merging on ', datestr(now, 31), ', by ', h.p.cid];
%   Temporary structure for saving
tmpsave = h.pd; %#ok<NASGU>
%Save data as original filename
save([h.fp.dir, h.fp.name, 'xy', h.fp.xyfun(h.xy.list(h.xy.i))], ...
    '-struct', 'tmpsave');
