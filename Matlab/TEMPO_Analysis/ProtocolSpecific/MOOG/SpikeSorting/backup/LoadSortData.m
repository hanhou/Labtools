function varargout = LoadSortData(varargin)
% LOADSORTDATA M-file for LoadSortData.fig
% Last Modified by GUIDE v2.5 26-Feb-2006 12:17:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoadSortData_OpeningFcn, ...
                   'gui_OutputFcn',  @LoadSortData_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin==0
    fig=openfig(mfilename,'reuse');
   %Generate a structure of handles to pass to callbacks, and store it. 
    handles=guihandles(fig);    
    handles.isFileOpen=0;
    clear global spsData2;
    guidata(fig,handles);   
    if nargout>0
        varargout{1}=fig;
    end
elseif nargin && ischar(varargin{1})
    %gui_State.gui_Callback = str2func(varargin{1});
    try 
        if(nargout)
            [varargout{1:nargout}]=feval(varargin{:});
        else
            feval(varargin{:});
        end
    catch
        disp(lasterr);
    end        
end

% if nargin && ischar(varargin{1})
%     gui_State.gui_Callback = str2func(varargin{1});
%     
% end

% if nargout
%     [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
% else
%     gui_mainfcn(gui_State, varargin{:});
% end
% End initialization code - DO NOT EDIT

% --- Executes just before LoadSortData is made visible.
function LoadSortData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoadSortData (see VARARGIN)

% Choose default command line output for LoadSortData
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LoadSortData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LoadSortData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
if (handles.isFileOpen==1)
    fclose(handles.spike2file);
    disp('Spike2 data file closed');
end
delete(gcf);



function edit_FileName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FileName as text
%        str2double(get(hObject,'String')) returns contents of edit_FileName as a double


% --- Executes during object creation, after setting all properties.
function edit_FileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseButton.
function BrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    dfltDir = 'Z:\Data\MOOG\Jedi\Analysis';% Default data directory
    dfltSuffix = '*.smr';      % Default file type.
    cd(dfltDir);

    % Locate the batch file and check if anything was chosen.
    [dataFileName, dataPath] = uigetfile(dfltSuffix, 'Choose Data File');
    if (dataFileName == 0)
        return;
    end
    set(handles.edit_FileName, 'String', [dataPath, dataFileName]);
    handles.dataFileName=dataFileName;
    handles.dataPath=dataPath;
    
    set(handles.OpenPlotButton,'Enable','on');
    
    guidata(hObject, handles);



% --- Executes on button press in ExportButton.
function ExportButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Set the sorted spikes output file
global spsData2;
% slashIndex=findstr(handles.dataPath,'\')
% OutFileName=[handles.dataPath
% handles.dataFileName(1:length(handles.dataFileName)-3),'.mat'];
% eval(['save ', OutFileName, ' spsData2']);


slashIndex = findstr(handles.dataPath, '\');
monkeyName = handles.dataPath(slashIndex(3)+1:slashIndex(4)-1);
OutFileName=['Z:\Data\Moog\', monkeyName, '\Analysis\SortedSpikes2\', handles.dataFileName(1:length(handles.dataFileName) - 3), 'mat']; 
eval(['save ', OutFileName, ' spsData2']);
guidata(hObject,handles);



% --- Executes on button press in OpenPlotButton.
function OpenPlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to OpenPlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get the file to open
fileName=get(handles.edit_FileName,'string')

%Check to see if a file was given
if(isempty(fileName))
    errordlg('You must choose a file to open!','File Error');
    return;
end

%Close any old open files
if (handles.isFileOpen==1)
    fclose(handles.spike2file);
    disp('Old file closed');
end

%Open up the Spike2 data file
handles.spike2file=fopen(fileName);
handles.isFileOpen=1;
disp('File openend successfully!');

%Load ADCMarker (WaveMark) data CHAN3
global CHAN3;
[CHAN3, CHAN3header]=SONGetADCMarkerChannel(handles.spike2file,3);
plot(CHAN3.adc');

%Get the Markers of the sorted neurons
Markers=double(CHAN3.markers(:,1));
[NeuronID,SpikeNumber]=munique(Markers)

%reject some neuron if the firing rate is too low
MaxSpikeNumber=max(SpikeNumber);
k=0;
for i=1:size(NeuronID,1)
    if NeuronID(i)>0
        if SpikeNumber(i)>0.1*MaxSpikeNumber
            k=k+1;
            SelectNeuronID(k)=NeuronID(i);        
        end
    end
end

%find the timings of the SelectNeuronID
for i=1:length(SelectNeuronID)
    Index=find(Markers(:,1)==SelectNeuronID(i));
    Neuron(i).SpikeTiming=CHAN3.timings(Index);
    %spsData2(i).SpikeTiming=CHAN3.timings(Index);
    clear Index;
end

%Convert to tempo format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Event channel 32 into memory
global CHAN32;
[CHAN32 CHAN32Header]=SONGetMarkerChannel(handles.spike2file,32);

index = [];
for (i = 1:length(CHAN32.timings))
    % If we find a stimuls start code, then we look to see if we find
    % a success code before we get to another start code.
    if(real(CHAN32.markers(i,1)) == 4)
        j = i + 1;
        while (j <= length(CHAN32.timings) & real(CHAN32.markers(j,1)) ~= 4)
            if (real(CHAN32.markers(j,1)) == 12)
                index = [index, i];
                break;
            else
                j = j + 1;
            end
        end
        i = j;
    end
end

global spsData2;
% Stuff the spsData2 struct array with information about each successfull trial.
for k=1:length(SelectNeuronID)
    PreEventBuffer=1;
    PostEventBuffer=4;
    spsData2(k).sampleRate = 25000;
    spsData2(k).prebuffer = round(PreEventBuffer * 25000);
    spsData2(k).postbuffer = round(PostEventBuffer * 25000);
    
    for (i = 1:length(index))   
        % This is basically descriptive data about the spike area analyzed.   
        spsData2(k).spikeInfo(i).startCodeTime = CHAN32.timings(index(i))*spsData2(k).sampleRate;     
        spsData2(k).spikeInfo(i).startTime = spsData2(k).spikeInfo(i).startCodeTime - spsData2(k).prebuffer + 1;    
        spsData2(k).spikeInfo(i).endTime = spsData2(k).spikeInfo(i).startCodeTime + spsData2(k).postbuffer;
        
        % Store all the event codes for later reference. 
        binCount = (spsData2(k).postbuffer + spsData2(k).prebuffer) /spsData2(k).sampleRate * 1000;
        slotsperbin = (spsData2(k).postbuffer + spsData2(k).prebuffer) / binCount;
        spsData2(k).spikeInfo(i).eventCodes = zeros(1, binCount);
        mrstart = spsData2(k).spikeInfo(i).startTime;
        mrend = spsData2(k).spikeInfo(i).endTime;
        mrsuckass = spsData2(k).sampleRate*[CHAN32.timings];      
        mrstupid = find(mrsuckass >= mrstart & mrsuckass <= mrend);
        
        a = [CHAN32.timings(mrstupid)]*spsData2(k).sampleRate;
        CHAN32markers=real(CHAN32.markers(:,1));
        spsData2(k).spikeInfo(i).eventCodes([ceil((a - spsData2(k).spikeInfo(i).startTime + 1) / 25)]) =CHAN32markers(mrstupid);
        
        %Find spike times of each trial
        clear mrsuckass;mrsuckass=Neuron(k).SpikeTiming'*spsData2(k).sampleRate;
        clear mrstupid;mrstupid=find(mrsuckass >= mrstart & mrsuckass <= mrend);
        spsData2(k).spikeInfo(i).SpikeTimes=mrsuckass(mrstupid);    
    end % End for (i = 1:length(index))    
end

guidata(hObject,handles);

% --- Executes on button press in ClosePlotButton.
function ClosePlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClosePlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Close the .smr file
fclose(handles.spike2file);
disp('Spike2 data file closed');
%Clear the axes and the channel data
cla;
clear global CHAN1;
clear global CHAN2;
clear global CHAN32;
clear global spsData2;
disp('Axes cleared and channel data wiped!')
handles.isFileOpen=0;

set(handles.edit_FileName, 'String', '');
set(handles.OpenPlotButton,'Enable','off');
guidata(hObject,handles);

