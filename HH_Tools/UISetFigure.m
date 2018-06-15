function varargout = UISetFigure(varargin)
% UISETFIGURE MATLAB code for UISetFigure.fig
%      UISETFIGURE, by itself, creates a new UISETFIGURE or raises the existing
%      singleton*.
%
%      H = UISETFIGURE returns the handle to a new UISETFIGURE or the handle to
%      the existing singleton*.
%
%      UISETFIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UISETFIGURE.M with the given input arguments.
%
%      UISETFIGURE('Property','Value',...) creates a new UISETFIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UISetFigure_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UISetFigure_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UISetFigure

% Last Modified by GUIDE v2.5 14-Jun-2018 20:41:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UISetFigure_OpeningFcn, ...
                   'gui_OutputFcn',  @UISetFigure_OutputFcn, ...
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


% --- Executes just before UISetFigure is made visible.
function UISetFigure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UISetFigure (see VARARGIN)

% Choose default command line output for UISetFigure
handles.output = hObject;

% UIWAIT makes UISetFigure wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% set(gcf,'HandleVisibility','off')
handles.showLabel = 0;
handles = UpdateAxes([],[],handles);

% Update handles structure
guidata(hObject, handles);


function handles = UpdateAxes(~,~, handles)
% Show all axes of current figure
allAxes = findobj(gcf,'type','axes');

for aa = 1:length(allAxes)
    thisAxis = allAxes(aa);
    handles.Axes.String{aa} = sprintf('Axes %g', aa);
    if handles.showLabel
        handles.AxesLabel(aa) = text(thisAxis, min(xlim(thisAxis)),max(ylim(thisAxis)),num2str(aa),'FontSize',30,'Color','m');
    end
end

if ~handles.showLabel && isfield(handles,'AxesLabel')
    delete(handles.AxesLabel);
end

handles.allAxes = allAxes;

% --- Outputs from this function are returned to the command line.
function varargout = UISetFigure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function FontSize_Callback(hObject, eventdata, handles)
% hObject    handle to FontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FontSize as text
%        str2double(get(hObject,'String')) returns contents of FontSize as a double


% --- Executes during object creation, after setting all properties.
function FontSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetFigure.
function SetFigure_Callback(hObject, eventdata, handles)
% hObject    handle to SetFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(gcbo,'HandleVisibility','off')

fontSize = str2double(handles.FontSize.String);
SetFigure(fontSize);


% --- Executes on selection change in Axes.
function Axes_Callback(hObject, eventdata, handles)
% hObject    handle to Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Axes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Axes


function Axes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ToggleLabel.
function ToggleLabel_Callback(hObject, eventdata, handles)
% hObject    handle to ToggleLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showLabel = ~ handles.showLabel;

handles = UpdateAxes([],[],handles);

tmp = gcf;
handles.CurrentFigure.String = sprintf('Current: %g',tmp.Number);

% Update handles structure
guidata(hObject, handles);


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to FontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FontSize as text
%        str2double(get(hObject,'String')) returns contents of FontSize as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetFigure.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to SetFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AxesEqualSize.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to AxesEqualSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in Axes.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Axes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Axes


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ToggleLabel.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to ToggleLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function BatchCommand_Callback(hObject, eventdata, handles)
% hObject    handle to BatchCommand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BatchCommand as text
%        str2double(get(hObject,'String')) returns contents of BatchCommand as a double


% --- Executes during object creation, after setting all properties.
function BatchCommand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BatchCommand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BatchRun.
function BatchRun_Callback(hObject, eventdata, handles)
% hObject    handle to BatchRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
command = handles.BatchCommand.String;
command = [command repmat(',',size(command,1),1)];
command = command';
commandMultiLine = command(:)';  % Enable multiplelines

% ax = findall(gcf,'type','axes');
axisTodo = handles.allAxes(handles.Axes.Value);

for aa = 1:length(axisTodo)
    axes(axisTodo(aa)); 
    eval(commandMultiLine);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text5.
function text5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guide UISetFigure
edit UISetFigure


% --- Executes on button press in printEPS.
function printEPS_Callback(hObject, eventdata, handles)
% hObject    handle to printEPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% To eps file
set(gcf,'paperpositionmode','auto');
% print('-depsc','-noui','-painters','Z:\matlab_print_eps.eps');
printeps(gcf,'Z:\matlab_print_eps');

% Fix the problem of rendering pathces after MATLAB2014
if ~isempty(findobj(gcf,'type','patch','-or','type','contour'))
%    epsclean('Z:\matlab_print_eps.eps','combineAreas',true...
%             ,'removeBoxes',true,'groupSoft',true);
%     epsclean('Z:\matlab_print_eps.eps','Z:\matlab_print_eps_clean.eps');
%     winopen('Z:\matlab_print_eps_clean.eps');
else
    %epsclean('Z:\matlab_print_eps.eps',...
    %   'removeBoxes',false,'groupSoft',false);
end

winopen('Z:\matlab_print_eps.eps');

% To clipboard
% print('-dmeta','-r600',gcf);


% --- Executes on button press in copyFig.
function copyFig_Callback(hObject, eventdata, handles)
% hObject    handle to copyFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
copyobj(gcf,0);


% --- Executes on button press in AxesEqualSize.
function AxesEqualSize_Callback(hObject, eventdata, handles)
% hObject    handle to AxesEqualSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(gcbo,'HandleVisibility','off')

% axes = findall(gcf,'type','axes');

axisTodo = handles.allAxes(handles.Axes.Value);

if length(axisTodo) >= 1
    pos = reshape([axisTodo.Position],4,[]);
    
    if ~isempty(gco) && strcmp(get(gco,'Type'),'axes')
        gg = gco;
        aver_size = gg.Position([3 4]);
    else
        aver_size = mean(pos([3 4],:),2)';
    end
    
    for i = 1:length(axisTodo)
        handles.LastDo.cachePos(i,:) = axisTodo(i).Position;
        axisTodo(i).Position([3 4]) = aver_size;
    end
end
handles.LastDo.axisTodo = axisTodo;

guidata(hObject,handles);

% --- Executes on button press in AlignUp.
function AlignUp_Callback(hObject, eventdata, handles)
% hObject    handle to AlignUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axisTodo = handles.allAxes(handles.Axes.Value);

if length(axisTodo) >= 1
    pos = reshape([axisTodo.Position],4,[]);
    
    if ~isempty(gco) && strcmp(get(gco,'Type'),'axes')
        gg = gco;
        commonUp = sum(gg.Position([2 4]));
    else
        commonUp = mean(sum(pos([2 4],:)),2);
    end
    
    for i = 1:length(axisTodo)
        handles.LastDo.cachePos(i,:) = axisTodo(i).Position;
        axisTodo(i).Position(2) = commonUp - axisTodo(i).Position(4);
    end
end
handles.LastDo.axisTodo = axisTodo;

guidata(hObject,handles);

% --- Executes on button press in AlignDown.
function AlignDown_Callback(hObject, eventdata, handles)
% hObject    handle to AlignDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axisTodo = handles.allAxes(handles.Axes.Value);

if length(axisTodo) >= 1
    pos = reshape([axisTodo.Position],4,[]);
    
    if ~isempty(gco) && strcmp(get(gco,'Type'),'axes')
        gg = gco;
        commonDown = gg.Position(2);
    else
        commonDown = mean(pos(2,:),2);
    end
    
    for i = 1:length(axisTodo)
        handles.LastDo.cachePos(i,:) = axisTodo(i).Position;
        axisTodo(i).Position(2) = commonDown;
    end
end
handles.LastDo.axisTodo = axisTodo;

guidata(hObject,handles);


% --- Executes on button press in AlignLeft.
function AlignLeft_Callback(hObject, eventdata, handles)
% hObject    handle to AlignLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axisTodo = handles.allAxes(handles.Axes.Value);

if length(axisTodo) >= 1
    pos = reshape([axisTodo.Position],4,[]);
    
    if ~isempty(gco) && strcmp(get(gco,'Type'),'axes')
        gg = gco;
        commonDown = gg.Position(1);
    else
        commonDown = mean(pos(1,:),2);
    end
    
    for i = 1:length(axisTodo)
        handles.LastDo.cachePos(i,:) = axisTodo(i).Position;
        axisTodo(i).Position(1) = commonDown;
    end
end
handles.LastDo.axisTodo = axisTodo;

guidata(hObject,handles);

% --- Executes on button press in AlignRight.
function AlignRight_Callback(hObject, eventdata, handles)
% hObject    handle to AlignRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axisTodo = handles.allAxes(handles.Axes.Value);

if length(axisTodo) >= 1
    pos = reshape([axisTodo.Position],4,[]);
    
    if ~isempty(gco) && strcmp(get(gco,'Type'),'axes')
        gg = gco;
        commonRight = sum(gg.Position([1 3]));
    else
        commonRight = mean(sum(pos([1 3],:)),2);
    end
    
    for i = 1:length(axisTodo)
        handles.LastDo.cachePos(i,:) = axisTodo(i).Position;
        axisTodo(i).Position(1) = commonRight - axisTodo(i).Position(3);
    end
end
handles.LastDo.axisTodo = axisTodo;

guidata(hObject,handles);

% --- Executes on button press in EqualHeight.
function EqualHeight_Callback(hObject, eventdata, handles)
% hObject    handle to EqualHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axisTodo = handles.allAxes(handles.Axes.Value);

if length(axisTodo) >= 1
    pos = reshape([axisTodo.Position],4,[]);
    
    if ~isempty(gco) && strcmp(get(gco,'Type'),'axes')
        gg = gco;
        commonHeight = sum(gg.Position(4));
    else
        commonHeight = mean(pos(4,:),2);
    end
    
    for i = 1:length(axisTodo)
        handles.LastDo.cachePos(i,:) = axisTodo(i).Position;
        axisTodo(i).Position(4) = commonHeight;
    end
end
handles.LastDo.axisTodo = axisTodo;

guidata(hObject,handles);

% --- Executes on button press in EqualWidth.
function EqualWidth_Callback(hObject, eventdata, handles)
% hObject    handle to EqualWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axisTodo = handles.allAxes(handles.Axes.Value);

if length(axisTodo) >= 1
    pos = reshape([axisTodo.Position],4,[]);
    
    if ~isempty(gco) && strcmp(get(gco,'Type'),'axes')
        gg = gco;
        commonWidth = sum(gg.Position(3));
    else
        commonWidth = mean(pos(3,:),2);
    end
    
    for i = 1:length(axisTodo)
        handles.LastDo.cachePos(i,:) = axisTodo(i).Position;
        axisTodo(i).Position(3) = commonWidth;
    end
end
handles.LastDo.axisTodo = axisTodo;

guidata(hObject,handles);


% --- Executes on button press in RollBack.
function RollBack_Callback(hObject, eventdata, handles)
% hObject    handle to RollBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'LastDo') && ~isempty(handles.LastDo)
    for aa = 1:length(handles.LastDo.axisTodo)
       handles.LastDo.axisTodo(aa).Position = handles.LastDo.cachePos(aa,:);
    end    
end
handles.LastDo = [];
guidata(hObject,handles);
