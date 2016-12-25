function varargout = LoadSortData_BATCH(varargin)
% LOADSORTDATA_BATCH M-file for LoadSortData_BATCH.fig
% Last Modified by GUIDE v2.5 24-Feb-2006 15:11:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoadSortData_BATCH_OpeningFcn, ...
                   'gui_OutputFcn',  @LoadSortData_BATCH_OutputFcn, ...
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


% --- Executes just before LoadSortData_BATCH is made visible.
function LoadSortData_BATCH_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoadSortData_BATCH (see VARARGIN)

% Choose default command line output for LoadSortData_BATCH
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LoadSortData_BATCH wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LoadSortData_BATCH_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_BatchFileName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_BatchFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_BatchFileName as text
%        str2double(get(hObject,'String')) returns contents of edit_BatchFileName as a double


% --- Executes during object creation, after setting all properties.
function edit_BatchFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BatchFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
    dfltDir = 'Z:\Users\Aihua\SpikeSorting\Batch';% Default data directory
    dfltSuffix = '*.m';      % Default file type.
    cd(dfltDir);

    % Locate the batch file and check if anything was chosen.
    [batchFileName, batchPath] = uigetfile(dfltSuffix, 'Choose Batch File');
    if (batchFileName == 0)
        return;
    end
    set(handles.edit_BatchFileName, 'String', [batchPath, batchFileName]);

    % Use textread to grab the path and filename of all files to be processed.
    [handles.filePath, handles.fileName] = textread([batchPath, batchFileName], '%s%s%*[^\n]', 'commentstyle', 'matlab');
    set(handles.Files_listbox,'String',[handles.fileName]);
    
    % Store the monkey's name.
    handles.monkeyName = batchPath(13:length(batchPath) - 1);
    guidata(hObject, handles);

% --- Executes on button press in StartButton.
function StartButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %Make sure these variables are clear
    clear global CHAN1;
    clear global CHAN2;
    clear global CHAN32;
    clear global spsData2;
    
    %Clear the FileList
    set(handles.FileList,'Value',1);
    set(handles.FileList,'String','');
    
    %buff=fopen('Z:\Users\Aihua\SpikeSorting\LoadSortData_BATCH.log','w');
    
    for (i=1:length(handles.fileName))
        fileList=get(handles.FileList,'String');
        listLength=length(fileList);
        handles.smrFileName=[handles.fileName{i}(1:length(handles.fileName{i})-4),'.smr'];
        handles.dataFileName=handles.smrFileName;
        
        %Add the file being processed to the list
        fileList{listLength+1}=['Analyzing --',handles.smrFilePath, handles.smrFileName,'...........'];
        set(handles.FileList,'String',fileList);
        set(handles.FileList,'Value',listLength+1);
        
        %Call all the processing functions
        disp('Loading spike2 sorted data');
        [handles,fixable]=LoadSpikeData(handles,i);
        
        %Go to the next file if there was a problem opeing the file.
        if(handles.spike2file<0 | fixable==0)
            %Show that the file was bad or corrupt
            if(fixable==0)
                fileList{ListLength+1}=[fileList{listLength+1},'Corrupt File'];
            else
                fileList{ListLength+1}=[fileList{listLength+1},'Bad File'];
            end
            %Updata the file list
            set(handles.FileList,'String',fileList);
            set(handles.FileList,'Value',listLength+1);
            
            %write out string to log file to keep record of results
            fprintf(buff, '%s \n', fileList{listLength + 1});
            continue;
        end

        %Close the .smr file
        fid()
        
        %Export spsData2
        global sps2Data2
        eval(['save ', handles.smrFilePath, '/SortedSpikes/', handles.smrFileName(1:length(handles.smrFileName) - 3), 'mat', ' spsData wdata']);

        % Clear all data channels.
        clear global CHAN1;
        clear global CHAN32;
        clear global spsData;
        
        % Show that the file is done being processed.
        fileList{listLength + 1} = [fileList{listLength + 1}, 'Finished'];
        
        % Update the file list.
        set(handles.FileList, 'String', fileList);
        set(handles.FileList, 'Value', listLength + 1);
        
        % write out string to log file to keep record of results
        fprintf(buff, '%s \n', fileList{listLength + 1});        
    end
    fclose(buff);
    
    % Indicate that everything is finished.
    fl = get(handles.FileList, 'String');
    fl{length(fl) + 1} = 'Done';
    set(handles.FileList, 'String', fl);
    set(handles.FileList, 'Value', length(fl));
    
    guidata(hObject, handles);    
    

 
% % ----------------------------------------------------------------------------
% %   Loads spike data from a file and processes it.
% % ----------------------------------------------------------------------------
% function [handles, fixable] = LoadSpikeData(handles, dumdex)
%     fixable = 1;
%     
%     % Look at the corresponding .htb file and load the pre/post buffer times.
%     handles.htdFile = [handles.filePath{dumdex}, handles.fileName{dumdex}];
%     fid = htbOpen(handles.htdFile);     % Open the htb file.
%     ndbs = htbCount(fid);    % Find out the # of databases it holds.
%     % Find the Events database and get the times.
%     for (i = 1:ndbs)
%         hd = htbGetHd(fid, i);   % Get the database header.
%         
%         if (strcmp(hd.title, 'Events'))
%             hertz = hd.speed_units / hd.speed;    % see p366 in TEMPO v9 manual
%    	        binwidth = (hd.skip + 1) / hertz;
% 	        epoch_start = hd.offset * binwidth;
%    	        epoch_period = hd.period * binwidth;
%             handles.PreEventBuffer = epoch_start;
%             handles.PostEventBuffer = -epoch_start + epoch_period;
%             break;
%         end
%     end
%     htbClose(fid);
%     
%     % Open up the Spike2 data file.
%     handles.spike2file = cedfunction('SonOpenOldFile', [handles.smrFilePath, handles.smrFileName], 1);
%     if (handles.spike2file >= 0)
%         disp('File opened successfully');
%     else
%         disp('Could not open smr file.');
%         return;
%     end
%     
%     % Load ADC channel 1 into memory.
%     global CHAN1;
%     handles.maxTime1 = cedfunction('SonChanMaxTime', handles.spike2file, 0);
%     handles.chand = cedfunction('SonChanDivide', handles.spike2file, 0);
%     numpts = cedfunction('SonGetNumADCPts', handles.spike2file, 0, 10000, 0, handles.maxTime1, handles.chand);
%     [CHAN1, num] = cedfunction('SonGetAllADCData', handles.spike2file, 0, numpts, 0, handles.maxTime1, handles.chand);
%     
%     % Load ADC Channel 2 into memory (Window Discriminator).
%     global CHAN2;
%     handles.maxTime2 = cedfunction('SonChanMaxTime', handles.spike2file, 1);
%     if (handles.maxTime2 > 0)
%         [num, CHAN2] = cedfunction('SonGetEventData', handles.spike2file, 1, 1000000, 0, handles.maxTime2);
%         handles.wData = 1;
%     else
%         disp('No Channel 2');
%         handles.maxTime2 = 0;
%         handles.wData = 0;
%     end
%     
%     % Load Event channel 32 into memory.
%     global CHAN32;
%     handles.maxTime32 = cedfunction('SonChanMaxTime', handles.spike2file, 31);
%     [num, CHAN32] = cedfunction('SonGetMarkData', handles.spike2file, 31, 1000000, 0, handles.maxTime32);
%     
%     % Determine max and min range values for the voltages.
%     maxVoltage = abs(double(max(CHAN1)) / (2^15 - 1) * 5);
%     minVoltage = abs(double(min(CHAN1)) / (2^15 - 1) * 5);
%     handles.rangeMax = round(max([maxVoltage, minVoltage])) + .5;
%     
%     % Close the open spike file.
%     cedfunction('SonCloseFile', handles.spike2file);
%     
%     % We call this function to fix any errors in CHAN32.
%     disp('Checking Data Integrity...');
%     [fixable, errsFound] = FixSMRData(handles);
%     
%     % Calculate the single unit spontaneous level from the window
%     % discriminator data.
%     if handles.maxTime2 > 0
%         handles.WinDiscrimSpontRate = FindWinDiscrimSpontRate(handles);
%     else
%         handles.WinDiscrimSpontRate = 0;
%     end
%     
%     disp(['SU spont rate: ', num2str(handles.WinDiscrimSpontRate)]);
%     
%     return;





% --- Executes on selection change in Files_listbox.
function Files_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Files_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Files_listbox

% --- Executes during object creation, after setting all properties.
function Files_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


