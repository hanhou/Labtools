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

% Last Modified by GUIDE v2.5 12-Sep-2016 17:36:21

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

% set(hObject,'toolbar','figure');
set(gcf,'HandleVis','off');
set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'off');
set(findall(handles.uipanel2, '-property', 'Enable'), 'Enable', 'off');

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

clear global group_result;

if ismember('control',get(gcf,'currentModifier'))
    winopen('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm');
    return;
end

set(handles.num_entries,'string',0); drawnow;

handles.XlsData = ReadXls('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',str2num(get(handles.sheetN,'string')),str2num(get(handles.headerN,'string')));
handles.N = size(handles.XlsData.num,1);

set(handles.num_entries,'string',num2str(handles.N));

set(findall(handles.uipanel2, '-property', 'Enable'), 'Enable', 'on');
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
set(handles.close_figs,'enable','on');
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

set([handles.level1 handles.level2 handles.fig_selected_txt],'string',[]);
set([handles.num_all_units],'string',[]);
handles.level1_selected = [];
handles.level2_selected = [];

data = handles.XlsData;
set(handles.close_figs,'enable','on');

set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findall(handles.cell_counter, '-property', 'Enable'), 'Enable', 'off');


% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425
function_handles = Group_Psychometric(data,get(handles.tolerance,'value'));
handles.function_handles = function_handles;
guidata(hObject,handles);
load_function_handles(hObject,handles); % Update possible function_handles


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

set([handles.level1 handles.level2 handles.fig_selected_txt],'string',[]);
% set([handles.num_all_units handles.num_su handles.num_t_cells],'string','0');
set([handles.num_all_units],'string',[]);

handles.level1_selected = [];
handles.level2_selected = [];

set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findall(handles.cell_counter, '-property', 'Enable'), 'Enable', 'on');

% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425
function_handles = Group_HD(data);

handles.function_handles = function_handles;
guidata(hObject,handles);
load_function_handles(hObject, handles);


function load_function_handles(hObject, handles)
% Clear previous values
set([handles.level1 handles.level2],'value',[]);

% Update possible function_handles
handles.level1_n = size(handles.function_handles,1)-1; % The last one reserved for functions that is accessible to GROUP_GUI 
                                                       % but not to manual selection. (cell_selection undapte) HH20150723
handles.level1_txt = cell(handles.level1_n,1);

handles.fig_all = [];
for i = 1:handles.level1_n
    handles.level1_txt{i} = handles.function_handles{i,1};
    level2_len = size(handles.function_handles{i,2},1);
    handles.fig_all = [handles.fig_all; repmat(i,level2_len,1) (1:level2_len)'];
end
set(handles.level1,'string',handles.level1_txt);

% Clear fig_selected_txt cache
handles.fig_selected = [];
guidata(hObject,handles);

cc = get(gcbf,'color');
set(gcbf,'color','y');
pause(0.1);
set(gcbf,'color',cc);
drawnow;

% --- Executes on selection change in level1.
function level1_Callback(hObject, eventdata, handles)
% hObject    handle to level1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns level1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from level1

if ~isfield(handles,'function_handles'); return; end

select_item = get(hObject,'value');

handles.level1_selected = select_item;

set(handles.level2,'value',[]);
handles.level2_selected = [];
handles.level2_txt = [];

if length(select_item) == 1
    handles.level2_n = size(handles.function_handles{select_item,2},1);
    for i = 1:handles.level2_n
        handles.level2_txt{i} = handles.function_handles{select_item,2}{i,1};
    end
end

set(handles.level2,'string',handles.level2_txt);
set(handles.debug,'enable','off');

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function level1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to level1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in level2.
function level2_Callback(hObject, eventdata, handles)
% hObject    handle to level2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns level2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from level2

handles.level2_selected = get(hObject,'value');
if length(handles.level2_selected) == 1
    set(handles.debug,'enable','on');
else
    set(handles.debug,'enable','off');
end

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function level2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to level2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fig_add.
function fig_add_Callback(hObject, eventdata, handles)
% hObject    handle to fig_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_add = [];

if ~isempty(handles.level2_selected) % Add selected level 2 entries
    fig_add = [repmat(handles.level1_selected,length(handles.level2_selected),1) handles.level2_selected'];    
else % Add all level 2 entries for each level 1 entry
    for j = 1:length(handles.level1_selected)
        level2_len = size(handles.function_handles{handles.level1_selected(j),2},1);
        fig_add = [fig_add ; repmat(handles.level1_selected(j),level2_len,1) (1:level2_len)'];
    end
end

handles.fig_selected = [handles.fig_selected; fig_add];
handles.fig_selected = sortrows(munique(handles.fig_selected));

guidata(hObject,handles);
update_fig_selected(handles);


function update_fig_selected(handles)
txt = [];
for i = 1:size(handles.fig_selected,1)
    l1 = handles.fig_selected(i,1);
    l2 = handles.fig_selected(i,2);
    txt{i} = sprintf('%s | %s\n',handles.function_handles{l1,1},handles.function_handles{l1,2}{l2,1});
end
set(handles.fig_selected_txt,'string',txt);
set(handles.fig_selected_txt,'value',[]);


% --- Executes on button press in fig_remove_all.
function fig_remove_all_Callback(hObject, eventdata, handles)
% hObject    handle to fig_remove_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fig_selected = [];
guidata(hObject,handles);
update_fig_selected(handles);


% --- Executes on selection change in fig_selected_txt.
function fig_selected_txt_Callback(hObject, eventdata, handles)
% hObject    handle to fig_selected_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fig_selected_txt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fig_selected_txt


% --- Executes during object creation, after setting all properties.
function fig_selected_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fig_selected_txt (see GCBO)
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


% --- Executes on button press in close_figs.
function close_figs_Callback(hObject, eventdata, handles)
% hObject    handle to close_figs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_figs = findall(0,'type','figure');
h_figs(h_figs == gcbf) = [];
close(h_figs);


% --- Executes on button press in fig_add_all.
function fig_add_all_Callback(hObject, eventdata, handles)
% hObject    handle to fig_add_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fig_selected = handles.fig_all;
guidata(hObject,handles);
update_fig_selected(handles);

% --- Executes on button press in fig_reverse.
function fig_reverse_Callback(hObject, eventdata, handles)
% hObject    handle to fig_reverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
all = handles.fig_all;
sel = handles.fig_selected;

if isempty(sel)
    handles.fig_selected = all;
else
    [~,reverse_ind] = setdiff(all , sel ,'rows');
    handles.fig_selected = all(reverse_ind,:);
end
guidata(hObject,handles);
update_fig_selected(handles);

% --- Executes on button press in fig_remove.
function fig_remove_Callback(hObject, eventdata, handles)
% hObject    handle to fig_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
to_remove = get(handles.fig_selected_txt,'value');
handles.fig_selected(to_remove,:) = [];
guidata(hObject,handles);
update_fig_selected(handles);

% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.fig_selected)
    if isempty(handles.level2_selected)
        to_do = handles.fig_all(ismember(handles.fig_all(:,1) ,handles.level1_selected),:);
    else
        to_do = [repmat(handles.level1_selected,length(handles.level2_selected),1) handles.level2_selected'];
    end
%     handles.fig_selected = to_do;
%     guidata(hObject,handles);
%     update_fig_selected(handles);
%    
else
    to_do = handles.fig_selected;
end

% Do analysis
cell_num = str2num(get(handles.num_all_units,'String'));

if isempty(cell_num) || cell_num(end,1) > 0 % More than one cell included
    for i = 1:size(to_do,1)
        feval(handles.function_handles{to_do(i,1),2}{to_do(i,2),2},false);  % Not debugging
    end
else
    beep; pause(0.1); beep;
end

figure(gcbf);


% --- Executes during object creation, after setting all properties.
function go_CreateFcn(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


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


% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
% hObject    handle to debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(handles.level2_selected) == 1
    feval(handles.function_handles{handles.level1_selected,2}{handles.level2_selected,2},true);  % Debug mode
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    dbquit('all');
catch
end


% --- Executes on button press in Polo_data.
function Polo_data_Callback(hObject, eventdata, handles)
% hObject    handle to Polo_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t_cell_selection_num = get(handles.t_criterion,'value'); 
feval(handles.function_handles{end,2}{1},t_cell_selection_num);  % Update cell_selection

% Hint: get(hObject,'Value') returns toggle state of Polo_data


% --- Executes on button press in Messi_data.
function Messi_data_Callback(hObject, eventdata, handles)
% hObject    handle to Messi_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t_cell_selection_num = get(handles.t_criterion,'value'); 
feval(handles.function_handles{end,2}{1},t_cell_selection_num);  % Update cell_selection

% Hint: get(hObject,'Value') returns toggle state of Messi_data


% --- Executes on button press in HD_dt.
function HD_dt_Callback(hObject, eventdata, handles)
% hObject    handle to HD_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to LIP_HD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ismember('control',get(gcf,'currentModifier'))
    edit Group_HD_dt;
    return;
end

data = handles.XlsData;

set([handles.level1 handles.level2 handles.fig_selected_txt],'string',[]);
% set([handles.num_all_units handles.num_su handles.num_t_cells],'string','0');
set([handles.num_all_units],'string',[]);

handles.level1_selected = [];
handles.level2_selected = [];

set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findall(handles.cell_counter, '-property', 'Enable'), 'Enable', 'on');

% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425
function_handles = Group_HD_dt(data);

handles.function_handles = function_handles;
guidata(hObject,handles);
load_function_handles(hObject, handles);


% --- Executes on selection change in t_criterion.
function t_criterion_Callback(hObject, eventdata, handles)  % HH20160912
% hObject    handle to t_criterion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t_cell_selection_num = get(handles.t_criterion,'value'); 
feval(handles.function_handles{end,2}{1},t_cell_selection_num);  % Update cell_selection

% Hints: contents = cellstr(get(hObject,'String')) returns t_criterion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t_criterion


% --- Executes during object creation, after setting all properties.
function t_criterion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_criterion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
