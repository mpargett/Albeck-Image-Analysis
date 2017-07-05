function varargout = iq_showspectraGUI(varargin)
% IQ_SHOWSPECTRAGUI MATLAB code for iq_showspectraGUI.fig
%      IQ_SHOWSPECTRAGUI, by itself, creates a new IQ_SHOWSPECTRAGUI or raises the existing
%      singleton*.
%
%      H = IQ_SHOWSPECTRAGUI returns the handle to a new IQ_SHOWSPECTRAGUI or the handle to
%      the existing singleton*.
%
%      IQ_SHOWSPECTRAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IQ_SHOWSPECTRAGUI.M with the given input arguments.
%
%      IQ_SHOWSPECTRAGUI('Property','Value',...) creates a new IQ_SHOWSPECTRAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iq_showspectraGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iq_showspectraGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iq_showspectraGUI

% Last Modified by GUIDE v2.5 13-Jan-2016 14:05:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iq_showspectraGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @iq_showspectraGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before iq_showspectraGUI is made visible.
function iq_showspectraGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iq_showspectraGUI (see VARARGIN)

% Choose default command line output for iq_showspectraGUI
handles.output = hObject;

%MP ------------------------------------------------------------------- MP
%Initialize and gather data, define indices and names.  Create user
%   data.
handles.Spectra = iq_getspectralpar();
%Define the Color Table for Wavelength-related plotting
handles.CTable = GetColorTable(handles.Spectra.WaveLength);

%Initialize lists and associated indices (all hidden)
set(handles.Hidden_FPhores, 'String', handles.Spectra.FPhoreNames(:,1));
Hud.Index = 1:size(handles.Spectra.FPhoreNames,1);
set(handles.Hidden_FPhores, 'UserData', Hud);
set(handles.Hidden_Filters, 'String', handles.Spectra.FilterNames(:,1));
Hud.Index = 1:size(handles.Spectra.FilterNames,1);
set(handles.Hidden_Filters, 'UserData', Hud);
Dud.Index = [];
set(handles.Disp_FPhores, 'UserData', Dud);
set(handles.Disp_Filters, 'UserData', Dud);

%Initialize Plot Axes (Labels)
xlabel(handles.EmAxes, 'WaveLength (nm)');
ylabel(handles.ExAxes, 'Excitation');
ylabel(handles.EmAxes, 'Emission');
%MP ------------------------------------------------------------------- MP

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iq_showspectraGUI wait for user response (see UIRESUME)
% uiwait(handles.MainFigure);


% --- Outputs from this function are returned to the command line.
function varargout = iq_showspectraGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Hidden_FPhores.
function Hidden_FPhores_Callback(hObject, eventdata, handles)
% hObject    handle to Hidden_FPhores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%MP ------------------------------------------------------------------- MP
ud = get(hObject,'UserData');
handles.Spectra.HidFPcurI = ud.Index(get(hObject,'Value'));
guidata(handles.MainFigure, handles);
%MP ------------------------------------------------------------------- MP

% Hints: contents = cellstr(get(hObject,'String')) returns Hidden_FPhores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Hidden_FPhores


% --- Executes during object creation, after setting all properties.
function Hidden_FPhores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hidden_FPhores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%MP ------------------------------------------------------------------- MP
% fpn = handles.Spectra.FPhoreNames(handles.Indices.FP.Hid);
% set(hObject,'String', fpn);  %Set names to list
%MP ------------------------------------------------------------------- MP

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Hidden_Filters.
function Hidden_Filters_Callback(hObject, eventdata, handles)
% hObject    handle to Hidden_Filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: contents = cellstr(get(hObject,'String')) returns Hidden_Filters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Hidden_Filters


% --- Executes during object creation, after setting all properties.
function Hidden_Filters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hidden_Filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PB_FilterShow.
function PB_FilterShow_Callback(hObject, eventdata, handles)
% hObject    handle to PB_FilterShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%MP ------------------------------------------------------------------- MP
%Get UserData, including Hidden List and Selected Value
Hud = get(handles.Hidden_Filters,'UserData');
Dud = get(handles.Disp_Filters,'UserData');
curI = get(handles.Hidden_Filters,'Value');
if ~isempty(Hud.Index)
    %Append to Displayed List
    Dud.Index = unique([Dud.Index, Hud.Index(curI)]);
    set(handles.Disp_Filters,'UserData',Dud);
    set(handles.Disp_Filters,'String', ...
        handles.Spectra.FilterNames(Dud.Index,1));
    %Remove from Hidden List
    Hud.Index(curI) = [];
    set(handles.Hidden_Filters,'String', ...
        handles.Spectra.FilterNames(Hud.Index,1));
    set(handles.Hidden_Filters,'UserData',Hud);
    set(handles.Hidden_Filters,'Value', 1);
    %Display on Axes
    UpdatePlot(handles)
end
%MP ------------------------------------------------------------------- MP



% --- Executes on button press in PB_FPhoreShow.
function PB_FPhoreShow_Callback(hObject, eventdata, handles)
% hObject    handle to PB_FPhoreShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%MP ------------------------------------------------------------------- MP
%Get UserData, including Hidden List and Selected Value
Hud = get(handles.Hidden_FPhores,'UserData');
Dud = get(handles.Disp_FPhores,'UserData');
curI = get(handles.Hidden_FPhores,'Value');
if ~isempty(Hud.Index)
    %Append to Displayed List
    Dud.Index = unique([Dud.Index, Hud.Index(curI)]);
    set(handles.Disp_FPhores,'UserData',Dud);
    set(handles.Disp_FPhores,'String', ...
        handles.Spectra.FPhoreNames(Dud.Index,1));
    %Remove from Hidden List
    Hud.Index(curI) = [];
    set(handles.Hidden_FPhores,'String', ...
        handles.Spectra.FPhoreNames(Hud.Index,1));
    set(handles.Hidden_FPhores,'UserData',Hud);
    set(handles.Hidden_FPhores,'Value', 1);
    %Display on Axes
    UpdatePlot(handles)
end
%MP ------------------------------------------------------------------- MP



% --- Executes on button press in PB_FPhoreHide.
function PB_FPhoreHide_Callback(hObject, eventdata, handles)
% hObject    handle to PB_FPhoreHide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%MP ------------------------------------------------------------------- MP
%Get UserData, including Hidden List and Selected Value
Hud = get(handles.Hidden_FPhores,'UserData');
Dud = get(handles.Disp_FPhores,'UserData');
curI = get(handles.Disp_FPhores,'Value');
if ~isempty(Dud.Index)
    %Append to Hidden List
    Hud.Index = unique([Hud.Index, Dud.Index(curI)]);
    set(handles.Hidden_FPhores,'UserData',Hud);
    set(handles.Hidden_FPhores,'String', ...
        handles.Spectra.FPhoreNames(Hud.Index,1));
    %Remove from Displayed List
    Dud.Index(curI) = [];
    set(handles.Disp_FPhores,'String', ...
        handles.Spectra.FPhoreNames(Dud.Index,1));
    set(handles.Disp_FPhores,'UserData',Dud);
    set(handles.Disp_FPhores,'Value', 1);
    %Display on Axes
    UpdatePlot(handles)
end
%MP ------------------------------------------------------------------- MP

% --- Executes on button press in PB_FilterHide.
function PB_FilterHide_Callback(hObject, eventdata, handles)
% hObject    handle to PB_FilterHide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%MP ------------------------------------------------------------------- MP
%Get UserData, including Hidden List and Selected Value
Hud = get(handles.Hidden_Filters,'UserData');
Dud = get(handles.Disp_Filters,'UserData');
curI = get(handles.Disp_Filters,'Value');
if ~isempty(Dud.Index)
    %Append to Hidden List
    Hud.Index = unique([Hud.Index, Dud.Index(curI)]);
    set(handles.Hidden_Filters,'UserData',Hud);
    set(handles.Hidden_Filters,'String', ...
        handles.Spectra.FilterNames(Hud.Index,1));
    %Remove from Displayed List
    Dud.Index(curI) = [];
    set(handles.Disp_Filters,'String', ...
        handles.Spectra.FilterNames(Dud.Index,1));
    set(handles.Disp_Filters,'UserData',Dud);
    set(handles.Disp_Filters,'Value', 1);
    %Display on Axes
    UpdatePlot(handles)
end
%MP ------------------------------------------------------------------- MP


% --- Executes on selection change in Disp_FPhores.
function Disp_FPhores_Callback(hObject, eventdata, handles)
% hObject    handle to Disp_FPhores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Update display on Axes
UpdatePlot(handles)

% Hints: contents = cellstr(get(hObject,'String')) returns Disp_FPhores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Disp_FPhores


% --- Executes during object creation, after setting all properties.
function Disp_FPhores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Disp_FPhores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Disp_Filters.
function Disp_Filters_Callback(hObject, eventdata, handles)
% hObject    handle to Disp_Filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Update display on Axes
UpdatePlot(handles)

% Hints: contents = cellstr(get(hObject,'String')) returns Disp_Filters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Disp_Filters


% --- Executes during object creation, after setting all properties.
function Disp_Filters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Disp_Filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on Disp_FPhores and none of its controls.
function Disp_FPhores_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Disp_FPhores (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


%MP: Plotting function, called after each update
function UpdatePlot(handles)

%Axes names
axn = {'ExAxes', 'EmAxes'};  
%Data fieldnames
dn = {'ab','ex';'em','em'};  %Axes x Elements
%Element names
en = {'Disp_FPhores', 'Disp_Filters'};  lstyles = {'-','--'};

for sa = 1:2  %FOR each Axis (Excitation/Emission)
    cla(handles.(axn{sa})); hold(handles.(axn{sa}), 'on');
    for ss = 1:2  %FOR each Element (FPhore/Filter)
        snames = get(handles.(en{ss}), 'String');    
        nshown = numel(snames);
        %Get currently selected element
        ci = get(handles.(en{ss}), 'Value');
        for s = 1:nshown
            k = s==ci;
            [~,clri] = max(handles.Spectra.(snames{s}).(dn{sa,ss}));
            wl = handles.Spectra.WaveLength(clri);
            lcr = interp1(300:800,handles.CTable,wl,'nearest')';
            plot(handles.(axn{sa}), handles.Spectra.WaveLength, ...
                handles.Spectra.(snames{s}).(dn{sa,ss}), ...
                lstyles{ss}, 'LineWidth', 2 + k*2, 'Color', lcr);
            if k
                plot(handles.(axn{sa}), handles.Spectra.WaveLength, ...
                    handles.Spectra.(snames{s}).(dn{sa,ss}), ...
                    '-', 'LineWidth', 1.5, 'Color', 'w');
            end
        end
        
    end

end

function ct = GetColorTable(wl)
wl_x = [300, 380, 440, 490, 510, 580, 645, 800];
r_x =  [  0,   1,   0,   0,   0,    1,  1,   0];
g_x =  [  0,   0,   0,   1,   1,    1,  0,   0];
b_x =  [  0,   1,   1,   1,   0,    0,  0,   0];

ct = interp1(wl_x, [r_x',g_x',b_x'], wl, 'linear', 'extrap');
