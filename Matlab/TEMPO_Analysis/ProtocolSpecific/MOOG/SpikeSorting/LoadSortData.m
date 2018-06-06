function varargout = LoadSortData(varargin)
% LOADSORTDATA M-file for LoadSortData.fig
% Last Modified by GUIDE v2.5 18-Mar-2015 11:02:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LoadSortDataz_OpeningFcn, ...
    'gui_OutputFcn',  @LoadSortData_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

if nargin==0 || nargin == 1  %HH20140830
    fig = openfig(mfilename,'reuse');
    set(fig,'toolbar','figure','name','=== Load Spike2 Data ===');
%     set(fig,'HandleVisibility','off');
    %Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    handles.isFileOpen=0;
    handles.spikeFigure = []; % HH20130905 Store the spike figures.
    
    if nargin==1
        set(handles.edit_FileName,'string',varargin{1});
    end
    
    clear global spsData2;
    
    guidata(fig,handles);
    if nargout>0
        varargout{1}=fig;
    end
    
elseif nargin && ischar(varargin{1})
    %gui_State.gui_Callback = str2func(varargin{1});
    %     try
    if(nargout)
        [varargout{1:nargout}]=feval(varargin{:});
    else
        feval(varargin{:});
    end
    %     catch error
    %         disp(lasterr);
%         beep;
%         keyboard;
%     end
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

fileName = GetFileName(hObject, eventdata, handles);

% The below line is important, because at the end of GetFileName(), we just updated the GUIDATA.
% If we do not have this line, the 'handles' here will fail to catch up with the changes made in 
% GetFilename(), e.g., handles.dataPath. As a result, when we do guidata()
% again at the end of this function, the changes we made in GetFileName()
% would be overwritten and missing.  HH20141005
handles = guidata(hObject); 

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
handles.toClip = [];
disp('File openend successfully!');

%Load ADCMarker (WaveMark) data CHAN3
global CHAN3;
%ChannelWaveMark=4;
ChannelWaveMark = str2num(get(handles.ChannelEdit,'String'));
% %handles.ChannelWaveMark=get(hObject,'String');
% set(handles.ChannelEdit, 'String', handles.ChannelWaveMark);
% guidata(hObject,handles);

[CHAN3, CHAN3header]=SONGetADCMarkerChannel(handles.spike2file,ChannelWaveMark);
%[CHAN3, CHAN3header]=SONGetADCMarkerChannel(handles.spike2file,3);
%Get the Markers of the sorted neurons
Markers=double(CHAN3.markers(:,1));
[NeuronID,SpikeNumber]=munique(Markers)

% [pc, zscores, pcvars] = princomp(double(CHAN3.adc));
% cumsum(pcvars./sum(pcvars) * 100);

% ==== Get MU from Chan21 if needed. HH20150209 ====
allChan = SONChanList(handles.spike2file);
Chan21 = 21;
Chan20 = 20;

if_explicit_chan21 = get(handles.chan21,'value') == 1 && ~isempty(find([allChan.number] == Chan21,1)); % If needed & exist
if_explicit_chan20 = get(handles.chan21,'value') == 1 && ~isempty(find([allChan.number] == Chan20,1)); % If needed & exist

if if_explicit_chan21
    MU_chan21 = SONGetEventChannel(handles.spike2file,Chan21);
end

if if_explicit_chan20
    MU_chan20 = SONGetEventChannel(handles.spike2file,Chan20);
end

% axes(handles.axes1);
% gscatter(zscores(:,1),zscores(:,2),Markers,'krgbcym','.******');
% xlabel('PC1');ylabel('PC2');
% title('Principal Component Scatter Plot with Colored Clusters');

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
SelectNeuronID = [];
for i=1:size(NeuronID,1)
    if NeuronID(i)>0
%         if SpikeNumber(i) > 0.01*MaxSpikeNumber
        if SpikeNumber(i) > 0.0 *MaxSpikeNumber   % HH20141115. Stick to CED
            k=k+1;
            SelectNeuronID(k) = NeuronID(i);
        end
    end
end

% Added two types of units:    @HH20150206  @HH20150410
% 1. If "Generate MU21 & 20 Explicitly" checked and Chan21/Chan20 exist, they are generated explicitly from Spike2 
% 2. Else, MU21 and MU20 will be generated from Chan5: MU20 = Marker 0; MU21 = Marker 0 + All SUs (all cross threshold)

MU20ID = 20;
SelectNeuronID = [SelectNeuronID MU20ID MU20ID+1]; 

%find the timings of the SelectNeuronID
for k=1:length(SelectNeuronID)
    
    switch SelectNeuronID(k)
        case MU20ID % Marker 0
           Index = find(Markers(:,1)==0);
        case MU20ID+1 % All cross threshold
           Index = find(Markers(:,1)>=0);
        otherwise % Real SUs
           Index = find(Markers(:,1)==SelectNeuronID(k));
    end
    
    if SelectNeuronID(k) == MU20ID && if_explicit_chan20  % If needed, we fetch MU from Chan20
        Neuron(k).SpikeTiming = MU_chan20;  % HH20150410
    elseif SelectNeuronID(k) == MU20ID + 1 && if_explicit_chan21  % If needed, we fetch MU from Chan20
        Neuron(k).SpikeTiming = MU_chan21;  % HH20150210
    else % Real SUs
        Neuron(k).SpikeTiming = CHAN3.timings(Index);
    end
    
    %spsData2(k).SpikeTiming=CHAN3.timings(Index);
    clear Index;
end

%Convert to tempo format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Event channel 32 into memory
global CHAN32;
[CHAN32 CHAN32Header]=SONGetMarkerChannel(handles.spike2file,32);

index = nan(1,sum(CHAN32.markers(:,1)==12));   % Store the indices of marker 04 of success trials
index_end = index;  % Store the indices of marker 05 of success trials. HH20130508
trial_condition = index; % Store trial condition number. HH20150409
good_trial_counter = 0;

for i = 1:length(CHAN32.timings)
    % If we find a stimuls start code, then we look to see if we find
    % a success code before we get to another start code. (i.e., only 'good trials')
    if(real(CHAN32.markers(i,1)) == 4)   % Visual begin code 04
        j = i + 1;
        while (j <= length(CHAN32.timings) && real(CHAN32.markers(j,1)) ~= 4)
            if (real(CHAN32.markers(j,1)) == 12)   % Success code 12
                % Now we know it's a success trial, and we save the indices of 04 and 05
                good_trial_counter = good_trial_counter + 1;
                index(good_trial_counter) = i;
                index_end(good_trial_counter) = tempEndIndex;
                
                % Get the condition number
                % Condition number starts from 48 (ASCII = '0') HH20150409
                trial_condition(good_trial_counter) = real(CHAN32.markers(find(real(CHAN32.markers(j:end,1))>=48,1)+j-1,1)) - 48;

                break;  % Break the whild loop. Let's move on to the next trial
                
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
for k= length(SelectNeuronID) :-1: 1
    
    PreEventBuffer=1; % s
    PostEventBuffer=4; % s
    
    spsData2(k).UnitId = SelectNeuronID(k);  % Make it clear. 10 = Marker 0; 11 = All cross threshold @HH20150207
    spsData2(k).sampleRate = 1/CHAN3header.sampleinterval;   % Make it soft-coded . HH20130508
    spsData2(k).prebuffer = round(PreEventBuffer *  spsData2(k).sampleRate);
    spsData2(k).postbuffer = round(PostEventBuffer *  spsData2(k).sampleRate);
    
    for i = 1:length(index) % Each trial
        % This is basically descriptive data about the spike area analyzed.
        spsData2(k).spikeInfo(i).startCodeTime = CHAN32.timings(index(i))*spsData2(k).sampleRate;  % 04
        spsData2(k).spikeInfo(i).endCodeTime = CHAN32.timings(index_end(i))*spsData2(k).sampleRate;   % 05
        spsData2(k).spikeInfo(i).startTime = spsData2(k).spikeInfo(i).startCodeTime - spsData2(k).prebuffer + 1;
        spsData2(k).spikeInfo(i).endTime = spsData2(k).spikeInfo(i).startCodeTime + spsData2(k).postbuffer;
        spsData2(k).spikeInfo(i).trialCondition = trial_condition(i);
        
        % Store all the event codes for later reference.
        binCount = (spsData2(k).postbuffer + spsData2(k).prebuffer) /spsData2(k).sampleRate * 1000;
        slotsperbin = (spsData2(k).postbuffer + spsData2(k).prebuffer) / binCount;
        spsData2(k).spikeInfo(i).eventCodes = zeros(1, round(binCount));
        mrstart = spsData2(k).spikeInfo(i).startTime;
        mrend = spsData2(k).spikeInfo(i).endTime;
        mrsuckass = spsData2(k).sampleRate*[CHAN32.timings];  % Who chose this ugly names? HH
        mrstupid = find(mrsuckass >= mrstart & mrsuckass <= mrend);
        
        a = [CHAN32.timings(mrstupid)]*spsData2(k).sampleRate;
        CHAN32markers=real(CHAN32.markers(:,1));
        spsData2(k).spikeInfo(i).eventCodes([ceil((a - spsData2(k).spikeInfo(i).startTime + 1) / 25)]) =CHAN32markers(mrstupid);
        
        %Find spike times of each trial
        clear mrsuckass; mrsuckass=Neuron(k).SpikeTiming'*spsData2(k).sampleRate;
        clear mrstupid; mrstupid=find(mrsuckass >= mrstart & mrsuckass <= mrend);
        spsData2(k).spikeInfo(i).SpikeTimes=mrsuckass(mrstupid);
        
        % Spike number that should be shown in TEMPO_GUI (often larger than the real number).  HH20130508
        spikeNumberInGUI(k)=spikeNumberInGUI(k)+length(spsData2(k).spikeInfo(i).SpikeTimes);

    end % End for (i = 1:length(index))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spike2 file m5c30r2.smr lost the first trial, so here I faked the 1st
    % trial with the trial of the same condition 
    % (Stim type = 1, heading = -12; outcome = correct, 63# in logfile, 54# success trial in the original smr file). HH20140531
    
    fakeSpikeData = {   % FileName     # Success trial in original Spike2
        'Z:\Data\MOOG\Polo\raw\m5c30r2.smr', 54;
        'Z:\Data\MOOG\Polo\raw\m5c77r1.smr', 5;
        'Z:\Data\MOOG\Polo\raw\m5c74r2.smr', 27;
        'Z:\Data\MOOG\Polo\raw\m5c70r2.smr', 12;
        'Z:\Data\MOOG\Polo\raw\m5c104r1.smr', 5;      
        'Z:\Data\MOOG\Polo\raw\m5c100r1.smr', 5;
        'Z:\Data\MOOG\Polo\raw\m5c80r2.smr', 6;
        'Z:\Data\MOOG\Polo\raw\m5c69r2.smr', 5;
        'Z:\Data\MOOG\Polo\raw\m5c26r1.smr', 6;
        'Z:\Data\MOOG\Polo\raw\m5c188r1.smr', 9;
        'Z:\Data\MOOG\Hetao\raw\m2c135r1.smr',8;
        'Z:\Data\MOOG\Polo\raw\m5c248r3.smr', 2 ;
        'Z:\Data\MOOG\Polo\raw\m5c295r2.smr', 9;
        'Z:\Data\MOOG\Messi\raw\m10c198r1.smr', 1;
        'Z:\Data\MOOG\Messi\raw\m10c230r2.smr', -2;
        };
    
    for ff = 1:size(fakeSpikeData,1)
        if strcmp(fileName,fakeSpikeData{ff,1})
            if fakeSpikeData{ff,2} > 0  % Fake the first trial
                temp = spsData2(k).spikeInfo(fakeSpikeData{ff,2});
                spsData2(k).spikeInfo(end+1) = spsData2(k).spikeInfo(end);
                spsData2(k).spikeInfo(2:end)=spsData2(k).spikeInfo(1:end-1);
                spsData2(k).spikeInfo(1) = temp;
                disp('Trial faked... please see the code to make sure... Type ''dbcont'' to continue');
                keyboard;
            else  % Fake the last trial
                spsData2(k).spikeInfo(end+1) = spsData2(k).spikeInfo(-fakeSpikeData{ff,2});
                disp('Trial faked... please see the code to make sure... Type ''dbcont'' to continue');
                keyboard;
            end
        end
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Rater Plot and PSTH.   HH20130508
    len = 0;
    h_f = figure(20+k);     clf;
    set(h_f,'position',[100 100 800 700])
    handles.spikeFigure = [handles.spikeFigure h_f];
    
    set(gcf, 'DefaultAxesXTickMode', 'auto','color','white');
%    subplot(2,2,1);
    axes('Position',[0.1 0.35 0.4 0.55]);
    box on;
    
    PSThist = zeros(1,100);
    
    % For raster plotting
    y = 1:length(spsData2(k).spikeInfo);
    t = -PreEventBuffer * 1000 : 1: PostEventBuffer * 1000; % each bin for raster plot = 1 ms
    [tt,yy] = meshgrid(t,y);
    spike_in_bin = zeros(length(y),length(t));
    
    for i = 1:length(spsData2(k).spikeInfo)
        hold on;
        
        if ~isempty(spsData2(k).spikeInfo(i).SpikeTimes)  % If there is more than 1 spike in this trial
            spikeTiming = (spsData2(k).spikeInfo(i).SpikeTimes-spsData2(k).spikeInfo(i).startCodeTime)/spsData2(k).sampleRate;
            
            spike_in_bin (i,fix(spikeTiming * 1000) + PreEventBuffer * 1000) = 1;
%             plot(spikeTiming,i,'b.','MarkerSize',5);

            n = histc(spikeTiming,linspace(-PreEventBuffer,PostEventBuffer,100));
            PSThist = PSThist + n;
        end
        
        % Plot the marker 05 (Vis Off)
        endThisTrial = (spsData2(k).spikeInfo(i).endCodeTime-spsData2(k).spikeInfo(i).startCodeTime)/spsData2(k).sampleRate;
        plot([endThisTrial endThisTrial],[i-0.3 i+0.3],'k','LineWidth',2);
     
        len = len + length(spsData2(k).spikeInfo(i).SpikeTimes);
    end
    
    % Use meshgrid to plot spike raster
    plot(tt(logical(spike_in_bin))/1000,yy(logical(spike_in_bin)),'.b');
    
    PSThist = PSThist/length(spsData2(k).spikeInfo)/((PostEventBuffer+PreEventBuffer)/100);
    
    plot([0 0],[0 i],'k--','LineWidth',2);
    title([get(handles.edit_FileName,'string') ' SU' num2str(SelectNeuronID(k)) ', N = ' ...
        num2str(length(Neuron(k).SpikeTiming)) '(GUI:' num2str(spikeNumberInGUI(k)) ')' ]);
    % text(0,i+10,'04-VisualOn');
    xlim([-PreEventBuffer,PostEventBuffer]);
    ylim([0 max(1,i)]);
    set(gca,'xtick',[]);

    axes('Position',[0.1 0.12 0.4 0.2]);
    h_bar = bar(linspace(-PreEventBuffer,PostEventBuffer,100), PSThist ,'histc');
    set(h_bar,'FaceColor','b','EdgeColor','b');
    hold on; plot([endThisTrial endThisTrial],[0 max(PSThist)*1.1],'k--','LineWidth',2);
    plot([0 0],[0 max(PSThist)*1.1],'k--','LineWidth',2);
    xlim([-PreEventBuffer,PostEventBuffer]);
    ylim([0 max(1,max(PSThist)*1.1)]);
    xlabel('Time to visual onset (s)');
    ylabel('Firing rate (Hz)');
    
    %% Spike Waveform Plot.   HH20130508
    colors = {'b','g','c','r','m'};
    
    if SelectNeuronID(k) < MU20ID  % Only plot waveform and ISI for real SUs
        subplot(2,2,2);
        waveForm = double(CHAN3.adc((CHAN3.markers(:,1)==SelectNeuronID(k)),:))';  % All the kth NEURON
        % waveForm = waveForm - repmat(mean(waveForm(1:20,:)),size(waveForm,1),1);
        meanWaveForm = mean(waveForm,2);
        stdWaveForm = std(waveForm,0,2);
        
        % Save spike waveform. HH20160920
        spsData2(k).meanWaveForm = meanWaveForm;
        spsData2(k).stdWaveForm = stdWaveForm;

        % Signal-to-Noise Ratio
        SNR(k) = range(meanWaveForm)/mean(stdWaveForm);

        hold on; plot((1:size(waveForm,1))*1000/spsData2(k).sampleRate,waveForm(:,1:ceil(size(waveForm,2)/1000):end),'color',[0.5 0.5 1]);
        plot((1:size(meanWaveForm,1))*1000/spsData2(k).sampleRate,meanWaveForm,'k','LineWidth',3);
        plot((1:size(stdWaveForm,1))*1000/spsData2(k).sampleRate,[meanWaveForm - stdWaveForm meanWaveForm + stdWaveForm],'k');
        set(gca,'ytick',[]);
        xlim([0 size(meanWaveForm,1)*1000/spsData2(k).sampleRate]);
        xlabel('(ms)');
        title(['SNR = ' num2str(SNR(k))]);
        box on;

        % ISI Plot. HH20130508
        axes('Position', [0.5703    0.1100    0.3347    0.1743]);
        timing = CHAN3.timings((CHAN3.markers(:,1)==SelectNeuronID(k)),:);
        hist(diff(timing)*1000,0:0.1:100);
        xlim([0 30]);
        xlabel('Inter-spike-interval (ms)');
        
        % Overdraw waveform
        axes(handles.axes1);
        hold on;
        hh = shadedErrorBar((1:size(meanWaveForm,1))*1000/spsData2(k).sampleRate,meanWaveForm,stdWaveForm,colors{SelectNeuronID(k)},0.6);
        set(hh.mainLine,'LineW',2);
        axis tight;
        
        figure(h_f);
    end
    
    %% Spike count over time
    axes('Position',[0.56625 0.334285714285714 0.345 0.171428571428571]);
    plot(1:size(spike_in_bin,1),sum(spike_in_bin,2)'); axis tight; box off;
    hold on;
    hh = axes('Position',[0.56625 0.504285714285714 0.345 0.01]);
    set(plot(CHAN32.timings(index),0),'visible','off'); axis tight;
    set(hh,'XAxisLocation','top','color','none','yticklabel','','box','off');

    % Output to command window for Excel. HH20130829
    
    if SelectNeuronID(k) < MU20ID  % Real SUs
        handles.toClip{k} = sprintf('SU%g\t %g\t %g\t %s\n',SelectNeuronID(k), spikeNumberInGUI(k), SNR(k), num2str((meanWaveForm'-min(meanWaveForm))/range(meanWaveForm)));
    else
        handles.toClip{k} = sprintf('MU%g\t %g\t %g\t %g\n',SelectNeuronID(k), spikeNumberInGUI(k), nan, nan);
    end
        
    fprintf(handles.toClip{k});
    
    if k == 1  % Copy the first unit to clipboard
        clipboard('copy',handles.toClip{k});
    end
    
    figure(gcbf);

end


% spikeNumberInGUI'

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
cla(handles.axes1);
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




% --- Executes on button press in OpenSpike2.
function OpenSpike2_Callback(hObject, eventdata, handles)
% hObject    handle to OpenSpike2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileName = GetFileName(hObject, eventdata, handles);
winopen(fileName);



function fileName = GetFileName(hObject, eventdata, handles)

fileName= get(handles.edit_FileName,'string');

%Check to see if a file was given
if(isempty(fileName))
    errordlg('You must choose a file to open!','File Error');
    return;
end

% Automatically add default path if you input cell No. by hand.  HH20130826
if isempty(findstr(fileName,'\'))    % If the directory has been ignored...
    fileName = horzcat(fileName(fileName~=' '),'.smr');
    
    monkeyN = str2num(fileName(strfind(fileName,'m')+1:strfind(fileName,'c')-1));
    switch monkeyN
        case 2
        handles.dataPath = 'Z:\Data\MOOG\Hetao\raw\';   
        case 5
        handles.dataPath = 'Z:\Data\MOOG\Polo\raw\';   
        case 10
        handles.dataPath = 'Z:\Data\MOOG\Messi\raw\';   

    end
    
    handles.dataFileName = fileName;
    fileName = horzcat(handles.dataPath,fileName) 
end

guidata(hObject, handles);   % This is important!!




% --- Executes on button press in chan21.
function chan21_Callback(hObject, eventdata, handles)
% hObject    handle to chan21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chan21


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
set(gcbf,'color','y');
pause(0.05);
set(gcbf,'color',[0.8 0.8 0.8]);
drawnow;

switch eventdata.Key
    case 'f2'
        OpenSpike2_Callback(hObject, eventdata, handles);
    case 'f3'
        OpenPlotButton_Callback(hObject, eventdata, handles);
    case 'f4'
        ExportButton_Callback(hObject, eventdata, handles);
        figure(1000);
end

if ~isempty(eventdata.Modifier) &&  ~isempty(eventdata.Character)
    k = str2double(eventdata.Key)
    try
        clipboard('copy',handles.toClip{k});
        beep;
    catch
    end
end




