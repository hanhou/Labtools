function [good_data, outputChan] = PackData(good_data,fName)

TEMPO_Defs;

%Load sorted spikes data
load(fName);

%If spsData2 is empty just return
if(isempty(spsData2))
    beep;
    disp('*** CAUTION ***: Empty Spike2 file!');
    outputChan = -1;
    return;
end

% ==== If the trial numbers of .htb and spsData do not match, something must be wrong. HH20141024 =====
% ==== Upgraded to make sure that each trial condition from TEMPO and Spike2 must exactly match. HH20150409

%  Condition number starts from 48 (ASCII = '0') in TEMPO protocol.  HH20150409
condition_TEMPO = good_data.event_data (good_data.event_data >= 48) - 48;  

if isfield(spsData2(1).spikeInfo,'trialCondition')  % Backward compatibility
    condition_Spike2 = [spsData2(1).spikeInfo.trialCondition];
    
    if sum(isnan(condition_Spike2)) > 0  % Also backward compatibility
        condition_Spike2 = -1 * ones(size(spsData2(1).spikeInfo,2),1);    
    end
else
    condition_Spike2 = -1 * ones(size(spsData2(1).spikeInfo,2),1);    
end

if length(condition_TEMPO) ~= length(condition_Spike2)  % Length does not match, must be wrong
    TEMPO_Spike2_match = 0;
elseif sum(condition_Spike2 == -1) > 0 % Length matches, but not clear about each trial condition
    TEMPO_Spike2_match = 1;
else  % Exact match on each trial
    TEMPO_Spike2_match = all(condition_TEMPO(:) == condition_Spike2(:));
end



if ~TEMPO_Spike2_match
    beep;
    disp('*** CAUTION ***: The trial conditions of htb and Spike2 do NOT match:');
    fprintf(' HTB = %g; Spike 2 = %g\n',length(condition_TEMPO),length(condition_Spike2));
    
    % Show condition list
    condition_list = nan(max( length(condition_TEMPO),length(condition_Spike2)),3);
    condition_list(1:length(condition_TEMPO),1) = condition_TEMPO;
    condition_list(1:length(condition_Spike2),2) = condition_Spike2;
    condition_list(:,3) = condition_list(:,1) == condition_list(:,2);
    disp(condition_list);
    
    disp('*** CAUTION ***: Please check for consistency!');
    
    if strfind(fName,'m5c181r3') % Exceptions
        disp('*** CAUTION ***: Spike2 sorted data STILL loaded');
    else
        disp('*** CAUTION ***: Spike2 sorted data not loaded');
        outputChan = -1;
        edit LoadSortData;
        return;
    end
else
    disp('Condition lists match exactly...');
end

%Find the next channel # to store data in.
chanNumberBegin = 4;  % This makes the first channel from Spike2 be #5. @HH20150206

%The number of bins and slots per bin
try
    binCount=round((spsData2(1).postbuffer+spsData2(1).prebuffer)/spsData2(1).sampleRate*(good_data.htb_header{SPIKE_DB}.speed_units/good_data.htb_header{SPIKE_DB}.speed)/2);%add /2 by ZC 2021/12/5 for 2sample/bin 
catch
    binCount=round((spsData2(1).postbuffer+spsData2(1).prebuffer)/spsData2(1).sampleRate*500);  % Workaround when there is no SPIKE trace in htb but we have spike2 signals. HH20130508
end

slotsperbin=(spsData2(1).postbuffer+spsData2(1).prebuffer)/binCount;
n=zeros(1,binCount);

% Find UnitID for Marker0ID @HH20150207
if isfield(spsData2,'UnitId')
    allID = [spsData2.UnitId];
    Marker0ID = allID(find(diff(allID)>1,1,'last')+1);
else % Backward compatibility
    allID = 1:length(spsData2);
    Marker0ID = Inf;
end

for i=1:length(spsData2)    
    for j=1:length(spsData2(i).spikeInfo)
        %Bin all the spike times
        edges=[spsData2(i).spikeInfo(j).startCodeTime-spsData2(i).prebuffer+1:slotsperbin:spsData2(i).spikeInfo(j).startCodeTime+spsData2(i).postbuffer+slotsperbin];
        combinedSpikes=spsData2(i).spikeInfo(j).SpikeTimes;
        %combinedSpikes = sort([spsData(index).spikeInfo(j).pSpikeTimes, spsData(index).spikeInfo(j).nSpikeTimes, spsData(index).spikeInfo(j).wSpikeTimes]);
        if(isempty(combinedSpikes))
            n=zeros(1,binCount);
        else
            n=histc(combinedSpikes,edges);
            n=n(1:length(n)-1);
        end
              
        % Put the bins in good_data
        if allID(i) < Marker0ID % Real SU
            good_data.spike_data(chanNumberBegin+allID(i),:,j)=n;  % Change it to cell. HH20150207
            outputChan(i) = chanNumberBegin+allID(i);
        else % Marker0 + AllAcrossThreshold
            good_data.spike_data(allID(i),:,j)=n;  % Change it to cell. HH20150207
            outputChan(i)= allID(i);
        end
    end
end


return;