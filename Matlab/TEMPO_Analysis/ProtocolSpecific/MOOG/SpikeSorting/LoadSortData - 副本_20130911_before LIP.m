function varargout = LoadSortData(varargin)
% LOADSORTDATA M-file for LoadSortData.fig
% Last Modified by GUIDE v2.5 10-Jul-2006 18:57:15

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
    handles.spikeFigure = []; % HH20130905 Store the spike figures.
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
    handles.isFileOpen = 0;
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


function varargout=ChannelEdit_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to ChannelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ChannelEdit as text
%        str2double(get(hObject,'String')) returns contents of ChannelEdit as a double

% --- Executes during object creation, after setting all properties.
function ChannelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelEdit (see GCBO)
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

%    dfltDir = 'Z:\Data\MOOG\Jedi\Analysis';% Default data directory
%    dfltDir = 'Z:\Data\MOOG\'; % Default data directory
dfltDir = 'Z:\Data\MOOG\Hetao\raw';% Default data directory for HH20130826
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

dataFileName = get(handles.edit_FileName,'string');
slashIndex = findstr(dataFileName,'\');
if ~isempty(slashIndex)
    dataFileName = dataFileName(slashIndex(end)+1:end-4);
end

dataFileName=dataFileName(dataFileName~=' ');
OutFileName=['Z:\Data\Moog\', monkeyName, '\Analysis\SortedSpikes2\', dataFileName, '.mat'];
if exist(OutFileName,'file')
    delete(OutFileName);
    disp('Overwriting...');
end
eval(['save ', OutFileName, ' spsData2']);
guidata(hObject,handles);
disp(['File Exported: ' OutFileName]);

% Window vibration
currPos = get(gcf,'Position');
for ii = 1:6
    set(gcf,'Position',[currPos(1) + 2*mod(ii,2) currPos(2) currPos(3) currPos(4)]);
    drawnow;
    pause(0.02);
end



% --- Executes on button press in OpenPlotButton.
function OpenPlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to OpenPlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get the file to open
%try

fileName=get(handles.edit_FileName,'string')


%Check to see if a file was given
if(isempty(fileName))
    errordlg('You must choose a file to open!','File Error');
    return;
end

% Automatically add default path if you input cell No. by hand.  HH20130826
if isempty(findstr(fileName,'\'))    % If the directory has been ignored...
    fileName = horzcat(fileName(fileName~=' '),'.smr');
    handles.dataPath = 'Z:\Data\MOOG\Hetao\raw\';
    handles.dataFileName = fileName;
    fileName = horzcat(handles.dataPath,fileName) 
end

% %Close any old open files
% if (handles.isFileOpen==1)
%     fclose(handles.spike2file);
%     disp('Old file closed');
% end

% Call closePlotButton directly. HH20130905
ClosePlotButton_Callback(hObject, eventdata, handles);

%Open up the Spike2 data file
handles.spike2file=fopen(fileName);
handles.isFileOpen=1;
disp('File openend successfully!');

%Load ADCMarker (WaveMark) data CHAN3
global CHAN3;
%ChannelWaveMark=4;
ChannelWaveMark=str2num(get(handles.ChannelEdit,'String'));
% %handles.ChannelWaveMark=get(hObject,'String');
% set(handles.ChannelEdit, 'String', handles.ChannelWaveMark);
% guidata(hObject,handles);

[CHAN3, CHAN3header]=SONGetADCMarkerChannel(handles.spike2file,ChannelWaveMark);
%[CHAN3, CHAN3header]=SONGetADCMarkerChannel(handles.spike2file,3);
%Get the Markers of the sorted neurons
Markers=double(CHAN3.markers(:,1));
[NeuronID,SpikeNumber]=munique(Markers)

[pc, zscores, pcvars] = princomp(double(CHAN3.adc));
cumsum(pcvars./sum(pcvars) * 100);
gscatter(zscores(:,1),zscores(:,2),Markers,'krgbcym','.******');
xlabel('PC1');ylabel('PC2');
title('Principal Component Scatter Plot with Colored Clusters');
% percent_explained = 100*pcvars/sum(pcvars);
% pareto(percent_explained)
% xlabel('Principal Component')
% ylabel('Variance Explained (%)')
%clear pc zscores pcvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the cluster figure;
%scatter plot
% figure;subplot(2,1,1);
% percent_explained = 100*pcvars/sum(pcvars);
% pareto(percent_explained)
% xlabel('Principal Component')
% ylabel('Variance Explained (%)')

% Color='krgbcym';
% az=[0 0 180];el=[90 0 90] ;
% for j=1:3
%     subplot(2,3,j);
%     for i=1:length(NeuronID)
%         clear index;index=find(Markers==NeuronID(i));
%         scatter3(zscores(index,1),zscores(index,2),zscores(index,3),[Color(i),'.']);hold on;
%     end
%     xlabel('PC1');ylabel('PC2');zlabel('PC3');
%     view(az(j),el(j));
%    % title('Principal Component Scatter Plot with Colored Clusters');
% end
% figure;subplot(4,1,1);plot(zscores(1,:));ylabel('PC1');
% subplot(4,1,2);plot(zscores(2,:));ylabel('PC2');
% subplot(4,1,3);plot(zscores(3,:));ylabel('PC3');
% clear pc zscores pcvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reject some neuron if the firing rate is too low
MaxSpikeNumber=max(SpikeNumber);
k=0;
for i=1:size(NeuronID,1)
    if NeuronID(i)>0
        if SpikeNumber(i)>0.01*MaxSpikeNumber
            k=k+1;
            SelectNeuronID(k)=NeuronID(i);
        end
    end
end

%find the timings of the SelectNeuronID
for k=1:length(SelectNeuronID)
    Index=find(Markers(:,1)==SelectNeuronID(k));
    Neuron(k).SpikeTiming=CHAN3.timings(Index);
    %spsData2(k).SpikeTiming=CHAN3.timings(Index);
    clear Index;
end

%Convert to tempo format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Event channel 32 into memory
global CHAN32;
[CHAN32 CHAN32Header]=SONGetMarkerChannel(handles.spike2file,32);

index = [];   % Store the indices of marker 04 of success trials
index_end = [];  % Store the indices of marker 05 of success trials. HH20130508

for i = 1:length(CHAN32.timings)
    % If we find a stimuls start code, then we look to see if we find
    % a success code before we get to another start code.
    if(real(CHAN32.markers(i,1)) == 4)   % Visual begin code 04
        j = i + 1;
        while (j <= length(CHAN32.timings) && real(CHAN32.markers(j,1)) ~= 4)
            if (real(CHAN32.markers(j,1)) == 12)   % Success code 12
                % Now we know it's a success trial, and we save the indices of 04 and 05
                index = [index, i];
                index_end = [index_end tempEndIndex];
                break;
            elseif (real(CHAN32.markers(j,1)) == 5)  % Visual end code 05
                tempEndIndex = j;  % Just store it here temporarily
                j = j + 1;  % Move on to the next marker
            else
                j = j + 1;
            end
        end
        i = j;  % Let's move on to the next trial
    end
end

spikeNumberInGUI = zeros(1:length(SelectNeuronID));
global spsData2;
disp('SpikeNum in GUI, SNR, meanWaveform:');

% Stuff the spsData2 struct array with information about each successfull trial.
for k=1:length(SelectNeuronID)
    
    PreEventBuffer=1; % s
    PostEventBuffer=4; % s
    spsData2(k).sampleRate = 1/CHAN3header.sampleinterval;   % Make it soft-coded . HH20130508
    spsData2(k).prebuffer = round(PreEventBuffer *  spsData2(k).sampleRate);
    spsData2(k).postbuffer = round(PostEventBuffer *  spsData2(k).sampleRate);
    
    for i = 1:length(index) % Each trial
        % This is basically descriptive data about the spike area analyzed.
        spsData2(k).spikeInfo(i).startCodeTime = CHAN32.timings(index(i))*spsData2(k).sampleRate;  % 04
        spsData2(k).spikeInfo(i).endCodeTime = CHAN32.timings(index_end(i))*spsData2(k).sampleRate;   % 05
        spsData2(k).spikeInfo(i).startTime = spsData2(k).spikeInfo(i).startCodeTime - spsData2(k).prebuffer + 1;
        spsData2(k).spikeInfo(i).endTime = spsData2(k).spikeInfo(i).startCodeTime + spsData2(k).postbuffer;
        
        % Store all the event codes for later reference.
        binCount = (spsData2(k).postbuffer + spsData2(k).prebuffer) /spsData2(k).sampleRate * 1000;
        slotsperbin = (spsData2(k).postbuffer + spsData2(k).prebuffer) / binCount;
        spsData2(k).spikeInfo(i).eventCodes = zeros(1, round(binCount));
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
        
        % Spike number that should be shown in TEMPO_GUI (often more than the real number).  HH20130508
        spikeNumberInGUI(k)=spikeNumberInGUI(k)+length(spsData2(k).spikeInfo(i).SpikeTimes);
        
    end % End for (i = 1:length(index))
    
    %% Rater Plot and PSTH.   HH20130508
    len = 0;
    h_f = figure(20+k);     
    set(h_f,'position',[200 200 800 700])
    handles.spikeFigure = [handles.spikeFigure h_f];
    
    set(gcf, 'DefaultAxesXTickMode', 'auto','color','white');
%    subplot(2,2,1);
    axes('Position',[0.1 0.35 0.4 0.55]);
    box on;
    
    PSThist = zeros(1,100);
    
    for i = 1:length(spsData2(k).spikeInfo)
        hold on;
        if ~isempty(spsData2(k).spikeInfo(i).SpikeTimes)  % If there is more than 1 spike in this trial
            spikeTiming = (spsData2(k).spikeInfo(i).SpikeTimes-spsData2(k).spikeInfo(i).startCodeTime)/spsData2(k).sampleRate;
            plot(spikeTiming,i,'b.','MarkerSize',5);
            n = histc(spikeTiming,linspace(-PreEventBuffer,PostEventBuffer,100));
            PSThist = PSThist + n;
        end
        % Plot the marker 05 (Vis Off)
        endThisTrial = (spsData2(k).spikeInfo(i).endCodeTime-spsData2(k).spikeInfo(i).startCodeTime)/spsData2(k).sampleRate;
        plot([endThisTrial endThisTrial],[i-0.3 i+0.3],'k','LineWidth',2);
     
        len = len + length(spsData2(k).spikeInfo(i).SpikeTimes);
    end
    
    PSThist = PSThist/length(spsData2(k).spikeInfo)/((PostEventBuffer+PreEventBuffer)/100);
    
    plot([0 0],[0 i],'k--','LineWidth',2);
    title([get(handles.edit_FileName,'string') ' SU' num2str(SelectNeuronID(k)) ', N = ' ...
        num2str(SpikeNumber(NeuronID==SelectNeuronID(k))) '(' num2str(spikeNumberInGUI(k)) ')' ]);
    % text(0,i+10,'04-VisualOn');
    xlim([-PreEventBuffer,PostEventBuffer]);
    ylim([0 i]);
    set(gca,'xtick',[]);

    axes('Position',[0.1 0.12 0.4 0.2]);
    h_bar = bar(linspace(-PreEventBuffer,PostEventBuffer,100), PSThist ,'style','histc');
    set(h_bar,'FaceColor','b','EdgeColor','b');
    hold on; plot([endThisTrial endThisTrial],[0 max(PSThist)*1.1],'k--','LineWidth',2);
    plot([0 0],[0 max(PSThist)*1.1],'k--','LineWidth',2);
    xlim([-PreEventBuffer,PostEventBuffer]);
    ylim([0 max(PSThist)*1.1]);
    xlabel('Time to visual onset (s)');
    ylabel('Firing rate (Hz)');
    
    %% Spike Waveform Plot.   HH20130508
    subplot(2,2,2);
    waveForm = double(CHAN3.adc((CHAN3.markers(:,1)==SelectNeuronID(k)),:))';  % All the kth NEURON
    % waveForm = waveForm - repmat(mean(waveForm(1:20,:)),size(waveForm,1),1);
    meanWaveForm = mean(waveForm,2);
    stdWaveForm = std(waveForm,0,2);
    
    % Signal-to-Noise Ratio
    SNR(k) = range(meanWaveForm)/mean(stdWaveForm);
    
    hold on; plot((1:size(waveForm,1))*1000/spsData2(k).sampleRate,waveForm,'color',[0.5 0.5 1]);
    plot((1:size(meanWaveForm,1))*1000/spsData2(k).sampleRate,meanWaveForm,'k','LineWidth',3);
    plot((1:size(stdWaveForm,1))*1000/spsData2(k).sampleRate,[meanWaveForm - stdWaveForm meanWaveForm + stdWaveForm],'k');
    set(gca,'ytick',[]);
    xlim([0 size(meanWaveForm,1)*1000/spsData2(k).sampleRate]);
    xlabel('(ms)');
    title(['SNR = ' num2str(SNR(k))]);
    box on;
    
    %% ISI Plot. HH20130508
    subplot(2,2,4);
    timing = CHAN3.timings((CHAN3.markers(:,1)==SelectNeuronID(k)),:);
    hist(diff(timing)*1000,0:0.1:100);
    xlim([0 10]);
    xlabel('Inter-spike-interval (ms)');

    % Output to command window for Excel. HH20130829
    fprintf('SU%g, %g, %g, %s\n',SelectNeuronID(k), spikeNumberInGUI(k), SNR(k), num2str((meanWaveForm'-min(meanWaveForm))/range(meanWaveForm)));
end

spikeNumberInGUI'

% catch lasterr
%     rethrow(lasterr);
%     keyboard
% end

%
guidata(hObject,handles);

%% --- Executes on button press in ClosePlotButton.
function ClosePlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClosePlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Close the .smr file
if handles.isFileOpen==1
    try fclose(handles.spike2file); end
    disp('Spike2 data file closed');
    handles.isFileOpen=0;
end

%Clear the axes and the channel data
cla;
clear global CHAN1;
clear global CHAN2;
clear global CHAN32;
clear global spsData2;
disp('Axes cleared and channel data wiped!')


% set(handles.edit_FileName, 'String', '');
% set(handles.ChannelEdit, 'String', '');
% set(handles.OpenPlotButton,'Enable','off');

if ~isempty(handles.spikeFigure)
    for i = 1:length(handles.spikeFigure)
        try close(handles.spikeFigure(i)); end
    end
    handles.spikeFigure=[];
end

guidata(hObject,handles);



%{

TRIAL_START_CD		= 1;	//trial started                                     01
FP_ON_CD			= 2;	//Fixation Point on                                     02
IN_FIX_WIN_CD		= 3;	//Entered fixation window                             03
VSTIM_ON_CD 		= 4;	//visual stimulus on                                  04
VSTIM_OFF_CD 		= 5;	//visual stimulus off                                 05
TARGS_ON_CD			= 6;	//targets on                                          06
SACCADE_BEGIN_CD	= 7;	//saccade has begun; monkey left fixation window    07
IN_T1_WIN_CD		= 8;	//monkey is in the T1 window                          08
IN_T2_WIN_CD		= 9;	//monkey is in the T2 window                          09
BROKE_FIX_CD		= 10;	//monkey broke fixation                               0A
BROKE_VERG_CD		= 11;	//monkey broke vergence                               0B
SUCCESS_CD			= 12;	//trial was successful                                0C
REWARD_CD			= 13;	//reward is delivered                                   0D
PUNISH_CD			= 14;	//beep is delivered                                     0E
TRIAL_END_CD		= 15;	//trial ended                                         0F
MICROSTIM_ON_CD 	= 16;	//microstim turned on                               10
MICROSTIM_OFF_CD 	= 17;	//microstim turned off                              11
MICROSTIM2_ON_HH = 18; //microstim2 turned on. HH120523                     12
MICROSTIMboth_ON_HH = 19; //both microstims turned on. HH120523             13

%}


