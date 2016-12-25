function varargout = GROUP_GUI(varargin)
% GROUP_GUI MATLAB code for GROUP_GUI.fig
%      GROUP_GUI, by itself, creates a new GROUP_GUI or raises the existing
%      singleton*.
%
%      H = GROUP_GUI returns the handle to a new GROUP_GUI or the handle to
%      the existing singleton*.
%
%      GROUP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GROUP_GUI.M with the given input arguments.
%
%      GROUP_GUI('Property','Value',...) creates a new GROUP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GROUP_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GROUP_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% // 
% //                                  _oo8oo_
% //                                 o8888888o
% //                                 88" . "88
% //                                 (| -_- |)
% //                                 0\  =  /0
% //                               ___/'==='\___
% //                             .' \\|     |// '.
% //                            / \\|||  :  |||// \
% //                           / _||||| -:- |||||_ \
% //                          |   | \\\  -  /// |   |
% //                          | \_|  ''\---/''  |_/ |
% //                          \  .-\__  '-'  __/-.  /
% //                        ___'. .'  /--.--\  '. .'___
% //                     ."" '<  '.___\_<|>_/___.'  >' "".
% //                    | | :  `- \`.:`\ _ /`:.`/ -`  : | |
% //                    \  \ `-.   \_ __\ /__ _/   .-` /  /
% //                =====`-.____`.___ \_____/ ___.`____.-`=====
% //                                  `=---=`
% // 
% //                               God Bless You!!
% //
% //                          
% //

% Edit the above text to modify the response to help GROUP_GUI

% Last Modified by GUIDE v2.5 16-Apr-2015 00:22:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GROUP_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GROUP_GUI_OutputFcn, ...
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


% --- Executes just before GROUP_GUI is made visible.
function GROUP_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GROUP_GUI (see VARARGIN)

% Choose default command line output for GROUP_GUI
handles.output = hObject;

set(hObject,'toolbar','figure');
set(gcf,'HandleVis','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GROUP_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GROUP_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ReadXls.
function ReadXls_Callback(hObject, eventdata, handles)
% hObject    handle to ReadXls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ismember('control',get(gcf,'currentModifier'))
    winopen('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm');
    return;
end

set(handles.num_entries,'string',0); drawnow;

handles.XlsData = ReadXls('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',str2num(get(handles.sheetN,'string')),str2num(get(handles.headerN,'string')));

% Clear mask cache
handles.N = size(handles.XlsData.num,1);
handles.mask = true(handles.N,1);

set(handles.items,'string',fieldnames(handles.XlsData.header));
set(handles.num_entries,'string',num2str(handles.N));

guidata(hObject,handles);


% --- Executes on button press in Microstimulation.
function Microstimulation_Callback(hObject, eventdata, handles)
% hObject    handle to Microstimulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ismember('control',get(gcf,'currentModifier'))
    edit Group_MicroStim;
    return;
end

data = handles.XlsData;
data.num = data.num(handles.mask,:);
data.txt = data.txt(handles.mask,:);
data.raw = data.raw(handles.mask,:);

Group_MicroStim(data);


% --- Executes on button press in Psychometric.
function Psychometric_Callback(hObject, eventdata, handles)
% hObject    handle to Psychometric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ismember('control',get(gcf,'currentModifier'))
    edit Group_Psychometric;
    return;
end

data = handles.XlsData;
data.num = data.num(handles.mask,:);
data.txt = data.txt(handles.mask,:);
data.raw = data.raw(handles.mask,:);

if str2num(get(handles.sheetN,'string'))==2
    Group_Psychometric(data);
else
    Group_Psychometric_Polo(data);
end

% --- Executes on button press in LIP_HD.
function LIP_HD_Callback(hObject, eventdata, handles)
% hObject    handle to LIP_HD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ismember('control',get(gcf,'currentModifier'))
    edit Group_HD;
    return;
end

data = handles.XlsData;
data.num = data.num(handles.mask,:);
data.txt = data.txt(handles.mask,:);
data.raw = data.raw(handles.mask,:);

Group_HD(data);


% --- Executes on selection change in items.
function items_Callback(hObject, eventdata, handles)
% hObject    handle to items (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns items contents as cell array
%        contents{get(hObject,'Value')} returns selected item from items
select_item = get(hObject,'value');
if sum(~isnan(handles.XlsData.num(:,select_item)))>0
    unique_item_values = unique(handles.XlsData.num(:,select_item));
else
    unique_item_values = unique(handles.XlsData.txt(:,select_item));
end
set(handles.values,'value',[]);
set(handles.values,'string',unique_item_values);


% --- Executes during object creation, after setting all properties.
function items_CreateFcn(hObject, eventdata, handles)
% hObject    handle to items (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in values.
function values_Callback(hObject, eventdata, handles)
% hObject    handle to values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns values contents as cell array
%        contents{get(hObject,'Value')} returns selected item from values


% --- Executes during object creation, after setting all properties.
function values_CreateFcn(hObject, eventdata, handles)
% hObject    handle to values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Mask_add.
function Mask_add_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
items = get(handles.items,{'value','string'});
item_select = items{2}{items{1}};

mask_temp = false(handles.N,1);

values = get(handles.values,{'value','string'});

for i = 1:length(values{1})
    if getfield(handles.XlsData.type,item_select) == 1 % Number type
        mask_temp = mask_temp | (handles.XlsData.num(:,items{1}) == str2num(values{2}(values{1}(i),:)));
    else
        mask_temp = mask_temp | (strcmp(handles.XlsData.txt(:,items{1}),values{2}{values{1}(i),:}));
    end
end

handles.mask = handles.mask & mask_temp;
set(handles.num_entries,'string',num2str(sum(handles.mask)));

guidata(hObject,handles);


% --- Executes on button press in Mask_clear.
function Mask_clear_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask = true(handles.N,1);
set(handles.num_entries,'string',num2str(sum(handles.mask)));

guidata(hObject,handles);

% --- Executes on selection change in info.
function info_Callback(hObject, eventdata, handles)
% hObject    handle to info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns info contents as cell array
%        contents{get(hObject,'Value')} returns selected item from info


% --- Executes during object creation, after setting all properties.
function info_CreateFcn(hObject, eventdata, handles)
% hObject    handle to info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sheetN_Callback(hObject, eventdata, handles)
% hObject    handle to sheetN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sheetN as text
%        str2double(get(hObject,'String')) returns contents of sheetN as a double


% --- Executes during object creation, after setting all properties.
function sheetN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sheetN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function headerN_Callback(hObject, eventdata, handles)
% hObject    handle to headerN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of headerN as text
%        str2double(get(hObject,'String')) returns contents of headerN as a double


% --- Executes during object creation, after setting all properties.
function headerN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to headerN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text1.
function text1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guide GROUP_GUI;
edit GROUP_GUI;


% --- Executes on selection change in monkey.
function monkey_Callback(hObject, eventdata, handles)
% hObject    handle to monkey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns monkey contents as cell array
%        contents{get(hObject,'Value')} returns selected item from monkey


% --- Executes during object creation, after setting all properties.
function monkey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to monkey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_figs = findall(0,'type','figure');
h_figs(h_figs == gcbf) = [];
close(h_figs);
