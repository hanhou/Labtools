function     ZCFUNC(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
DataMatrixDef;
trials = 1:size(data.moog_params,2);
select_trials = false(size(trials));
if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
    select_trials(BegTrial:EndTrial) = true;
elseif all(BegTrial > 0) % To be included
    select_trials(BegTrial) = true;
elseif all(BegTrial < 0) % To be excluded
    select_trials(-BegTrial) = true;
    select_trials = ~ select_trials;
else
    disp('Trial selection error...');
    keyboard;
end
% Time information
h = data.htb_header{SPIKE_DB};
spike_timeWin = 1000 * (h.skip + 1) / (h.speed_units / h.speed); % in ms % HH20140520

% Trial information
[temp_a,temp_t1,temp_t2,temp_v,temp_d,Mon_Outcome,RightAns,Mon_Choice,total_trials,One_target_trail,star_Denstity,star_Coherence,star_LifeTime,stim_type,sa,sas,Tcof] = ZC_Trail_Info(data)

distance_per_trial = temp_d(select_trials);
Outcome_per_trail = Mon_Outcome(select_trials);
unique_distance = munique(distance_per_trial');
unique_distance(isnan(unique_distance))=[];
repetitionN = floor(length(distance_per_trial) /length(unique_distance)) ;


% Spike data
spike_in_bin = squeeze(data.spike_data(SpikeChan,:,select_trials));   % 5000 * TrialNum for default

spike_in_bin(spike_in_bin > 1) = 1; % Some MU from Chan21 would have more than one spike per timebin. HH20150210

spikeN_trial = sum(spike_in_bin,1);
spikeN_cum_trial = [0 cumsum(spikeN_trial)];   % Find the start & end spikes number of each trial (for later indexing use)
%spikeT_raw = mod(find(spike_in_bin>0),size(spike_in_bin,1));  % Fast way to find spike times
[spikeT_raw,J] = find(spike_in_bin>0)
% Event data
events_in_bin = squeeze(data.event_data(1,:,select_trials));

unique_stim_type = 2; %？

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
align_markers = [ VSTIM_ON_CD, VSTIM_OFF_CD, SACCADE_BEGIN_CD];  % Desired markers: target onset & target offset & sac onset
psth_windows = [-500 2000; -1000 1500; -2000 500]; % Pre-trigger and Post-trigger PSTH windows (cover all possible range)
plot_markers = [ SACCADE_BEGIN_CD];  % Just plot one figure for simplicity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


binSize = 10; % in ms
stepSize = 10; % in ms
smoothFactor = 50; % in ms (Gaussian kernel)



% Anne 2014 NN: 10 ms bin + 50 ms Gaussian smoothing. HH20141207
binSize = 10; % in ms
stepSize = 10; % in ms
smoothFactor = 50; % in ms (Gaussian kernel)

for j = 1:length(align_markers)    % For each desired marker
    [align_offsets(:,j),~] = find(events_in_bin == align_markers(j));  % Fast way to find offsets of each trial
    
    % --- New Smoothed PSTH for each trial. HH20141207 ---
    % Align spikes    
    spike_aligned{j} = nan(sum(select_trials),ceil((psth_windows(j,2) - psth_windows(j,1)) / spike_timeWin) + 1);
    
    for i = 1:sum(select_trials)
        winBeg = align_offsets(i,j) + ceil(psth_windows(j,1) / spike_timeWin);
        winEnd = align_offsets(i,j) + ceil(psth_windows(j,2) / spike_timeWin);
        if winBeg > 0 && winEnd <= size(spike_in_bin,1) % OK time windows
            spike_aligned{j}(i,:) = spike_in_bin (winBeg : winEnd,i)';   % Align each trial
        else
            beep;
            fprintf('WARNING: PSTH window error at spike_aligned{%g}(%g,:)\n',j,i);
%             error('PSTH window error at spike_aligned: Line 207');
%             winBeg
%             spike_aligned{j}(i, 2-winBeg : end) = spike_in_bin (1 : winEnd,i)';
            % spike_aligned{j}(i,:) = NaN;
        end
    end
    
    for k = 1:length(t_centers{j})
        winBeg = ceil(((k-1) * stepSize) / spike_timeWin) + 1;
        winEnd = ceil(((k-1) * stepSize + binSize) / spike_timeWin) + 1;
        spike_hist{j}(:,k) = sum(spike_aligned{1,j}(: , winBeg:winEnd),2) / binSize*1000 ;  % in Hz
    end
    
    %  Smoothing 
    % Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
    % Note that this only influence PSTH calculation (spike_hist), not CP
    for i = 1:size(spike_hist{j},1) % Each trial
        spike_hist{j}(i,:) = GaussSmooth(t_centers{j},spike_hist{j}(i,:),smoothFactor);
    end    
end

%Long VS Short 
%Only correct trials.
LCidx = distance_per_trial  > Stand_D & Outcome_per_trail == 1;
SCidx = distance_per_trial <= Stand_D & Outcome_per_trail == 1;

for j = 1:length(align_markers)    % For each desired marker
    SpikeLC{j} = spike_hist{j}(LCidx,:)








all_data = loadzcdata(FILE);
[allData{1},allData{2},allData{3}] = SepByStimType(all_data);


%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
SAVPATH = 'E:\Data\Gutou\Raw\';
FILENAME = [FILE(1:end-4) '.mat']; 
data_extract(good_data,FILENAME,SAVPATH);
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

a = [FILE_list.name]