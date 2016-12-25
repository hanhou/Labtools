function LoadSortedSpikes_cah()
%Load Spike2 sorted spikes and save in matlab format to run tempo_gui
%AHC 02-22-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Opening a file
dfltSuffix = '*.smr';
[dataFileName, dataPath] = uigetfile(dfltSuffix, 'Choose Data File');
fid=fopen([dataPath dataFileName]);

%Read ADCMarker(WaveMark)data CHAN3
[CHAN3,header]=SONGetADCMarkerChannel(fid,3);
%plot(data.adc(5:10,:)')

%Get the Markers of the sorted neurons
Markers=double(CHAN3.markers(:,1));
[NeuronID, SpikeNumber] = munique(Markers);

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
% fileName='m7c73.smr';
% fid=fopen(fileName)
[CHAN32 CHAN32Header]=SONGetMarkerChannel(fid,32);

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

% Stuff the spsData struct array with information about each successfull trial.
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

%Set the sorted spikes output file
slashIndex = findstr(dataPath, '\');
monkeyName = dataPath(slashIndex(3)+1:slashIndex(4)-1);
OutFileName=['Z:\Data\Moog\', monkeyName, '\Analysis\SortedSpikes2\', dataFileName(1:length(dataFileName) - 3), 'mat']; 
eval(['save ', OutFileName, ' spsData2']);

%close the smr file
fclose(fid)
disp('Spike2 data file closed');