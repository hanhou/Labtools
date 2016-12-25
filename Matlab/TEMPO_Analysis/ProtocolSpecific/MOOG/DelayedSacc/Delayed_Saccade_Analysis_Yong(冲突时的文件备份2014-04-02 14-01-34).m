function Delayed_Saccade_Analysis_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs;

%% Commented by HH201310123
%%{
tic;

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME,:);
temp_spike_data = data.spike_data(SpikeChan,:);
temp_event_data = data.spike_data(2,:);
temp_spike_rates = data.spike_rates(SpikeChan,:);
% temp_accel = data.eye_data(5,:,:);
% temp_diode = data.eye_data(8,:,:);
temp_eyex = data.eye_data(3,:,:);
temp_eyey = data.eye_data(4,:,:);
% accel = reshape(temp_accel, 1000, length(temp_total_trials));
% diode = reshape(temp_diode, 1000, length(temp_total_trials));
eyex=reshape(temp_eyex, 1000, length(temp_total_trials));
eyey=reshape(temp_eyey, 1000, length(temp_total_trials));

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% vector of trial indices
% select_trials = ((trials >= BegTrial) & (trials <= EndTrial));

stim_type = temp_stim_type(trials);
heading = temp_heading(trials);
amplitude= temp_amplitude(trials);
num_sigmas= temp_num_sigmas(trials);
total_trials = temp_total_trials(trials);
spike_rates = temp_spike_rates(trials);
% unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');




timebin=40;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(total_trials);
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

time_step=1;

%Realign data to target onset and saccade onset
for k=1: length(unique_amplitude)
    for i=1:length(unique_heading)
        
        select = logical( (heading==unique_heading(i)) & (amplitude==unique_amplitude(k)) );
        act_found = find( select==1 );
        temp_count = [];

        for repeat=1:length(act_found)
            for n=1:(x_length)
                temp_count(repeat,n)=sum(temp_spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)++n*timebin)));
                time_step=time_step+timebin;
            end
            time_step=1;
        end
        count_y_trial{i,k}(:,:) = temp_count;  % align to the target onset
        count_y_mean{i,k}(1,:) = mean(temp_count);
        max_count(i) = max(count_y_mean{i,k}(1,:));
    end
end

% plot figures
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', 'PSTH','color','w');
orient landscape;
%title(FILE);

stim_start_moog = find(data.event_data(1,:,1) == 4);
stim_over_moog = find(data.event_data(1,:,1) == 5);

x_time = x_time - stim_start_moog/timebin*ones(1,length(x_time));  % in x_time

% x_start = [stim_start_moog/timebin, stim_start_moog/timebin];
% x_stop =  [stim_over_moog/timebin,  stim_over_moog/timebin];

x_start = [0 0];
x_stop =  [(stim_over_moog-stim_start_moog)/timebin,  (stim_over_moog-stim_start_moog)/timebin];
y_marker=[0,max(max_count)];

x_time = x_time * (5 / timebin); % in s (total length 5s, hard coded temporarily. HH20130908)
x_stop = x_stop * (5 / timebin);

positions = [6,3,2,1,4,7,8,9];

for i = 1:8
    
    subplot(3,3,positions(i)); % saccade to right/0
    bar( x_time,count_y_mean{i,1}(1,:) );
    hold on;
    plot( x_start , y_marker, 'r-');
    plot( x_stop,  y_marker, 'r-');
    xlim([-1,5]);
    ylim([0,max(max_count)]);
    box on;
       
end
toc;

%}

%% Added by HH20131023
%%{

% tic;

% Trial information
heading_per_trial = data.moog_params(HEADING,:,MOOG);
unique_heading = munique(heading_per_trial');

% Spike data
spike_in_bin = squeeze(data.spike_data(SpikeChan,:,:));   % 5000 * TrialNum for default
spikeN_trial = sum(spike_in_bin);
spikeN_cum_trial = [0 cumsum(spikeN_trial)];   % Find the start & end spikes number of each trial (for later indexing use)
spikeT_raw = mod(find(spike_in_bin>0),size(spike_in_bin,1));  % Fast way to find spike times

% Event data
events_in_bin = squeeze(data.event_data(1,:,:));
align_markers = [4 7];  % Desired markers: visual onset & sac onset
for j = 1:length(align_markers)    % For each desired marker
    align_offsets(:,j) = mod(find(events_in_bin == align_markers(j)),size(events_in_bin,1));  % Fast way to find offsets of each trial
end

numBins = 50;
result_PSTH = cell(length(align_markers),length(unique_heading));

% Now we align spike time to each marker, for different conditions
for i = 1:length(unique_heading)  % For each heading (delayed-saccade direction)
    select_binary = logical(heading_per_trial == unique_heading(i));
    select_trialN = find(select_binary == 1);
    
    spike_cache = NaN(spikeN_cum_trial(end),length(align_markers));       % Pre-allocation for speed
    spike_cache_numbers = 0;   % Number of spikes already in the cache
    
    for rep = 1:length(select_trialN)  % For each repeats
        this_trialN = select_trialN(rep);
        
        % Put spike times of this trial into cache (with desired offsets)
        spike_cache(spike_cache_numbers+1 : spike_cache_numbers+spikeN_trial(this_trialN) , :) ...
            = repmat(spikeT_raw(spikeN_cum_trial(this_trialN)+1 : spikeN_cum_trial(this_trialN+1)),1,length(align_markers)) ...   % Raw spike time for this trial
            - repmat(align_offsets(this_trialN,:),spikeN_trial(this_trialN),1); % Minus time offsets for j'th marker
        
        spike_cache_numbers = spike_cache_numbers + spikeN_trial(this_trialN);
        
    end
    
    % PSTH
    for j = 1:length(align_markers)    % For each desired marker
        [result_PSTH{j,i}, x(:,j)] = hist(spike_cache(:,j),numBins);
    end
    
end

% Plotting
set(0,'DefaultAxesXTickMode','auto');
for j = 1:length(align_markers)
    
    figure(3+j-1);
    set(3+j-1,'Position', [5,5 1000,700], 'Name', 'PSTH','color','w');
    orient landscape;
    positions = [6,3,2,1,4,7,8,9];
    maxPSTH = 0;
    
    for i = 1:length(unique_heading)
        
        subplot(3,3,positions(i)); % saccade to right/0
        bar( x(:,j), result_PSTH{j,i} );
        maxPSTH = max(maxPSTH, max(result_PSTH{j,i}));
        hold on;
        
    end
%     toc;
    
    for i = 1:length(unique_heading)
        subplot(3,3,positions(i));
        xlim([min(x(:,j)),max(x(:,j))]);
        ylim([0, maxPSTH*1.1]);        
        plot([0 0],[0 maxPSTH*1.1],'r-');
    end
end


%}

%%
return;
