function Memory_Saccade_Analysis_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);

TEMPO_Defs;
Path_Defs;
% ProtocolDefs;

%% Commented by HH201310123
%{
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

trials = 1:size(data.moog_params,2);		% a vector of trial indices

% If length(BegTrial) > 1 and all elements are positive, they are trials to be included.
% Else, if all elements are negative, they are trials to be excluded.
% This enable us to exclude certain trials ** DURING ** the recording more easily. HH20150410
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
heading_per_trial = data.moog_params(HEADING,select_trials,MOOG);

unique_heading = munique(heading_per_trial');
unique_heading(isnan(unique_heading))=[];
repetitionN = floor(length(heading_per_trial) /length(unique_heading)) ;

% Spike data
spike_in_bin = squeeze(data.spike_data(SpikeChan,:,select_trials));   % 5000 * TrialNum for default

spike_in_bin(spike_in_bin > 1) = 1; % Some MU from Chan21 would have more than one spike per timebin. HH20150210

spikeN_trial = sum(spike_in_bin,1);
spikeN_cum_trial = [0 cumsum(spikeN_trial)];   % Find the start & end spikes number of each trial (for later indexing use)
%spikeT_raw = mod(find(spike_in_bin>0),size(spike_in_bin,1));  % Fast way to find spike times
[spikeT_raw,J] = find(spike_in_bin>0)
% Event data
events_in_bin = squeeze(data.event_data(1,:,select_trials));

unique_stim_type = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
align_markers = [ VSTIM_ON_CD, VSTIM_OFF_CD, SACCADE_BEGIN_CD];  % Desired markers: target onset & target offset & sac onset
psth_windows = [-500 2000; -1000 1500; -2000 500]; % Pre-trigger and Post-trigger PSTH windows (cover all possible range)
plot_markers = [ SACCADE_BEGIN_CD];  % Just plot one figure for simplicity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% numBins = 50; % Fixed bin numbers

% Fixed bin size = 50 ms

% SPEED should be 1000 in TEMPO!! %% This is no longer necessary.

% HH20140520
binSize_old = 40; % Old time bin

% Anne 2014 NN: 10 ms bin + 50 ms Gaussian smoothing. HH20141207
binSize = 10; % in ms
stepSize = 10; % in ms
smoothFactor = 50; % in ms (Gaussian kernel)


% spikeSum = sum(spike_in_bin,2);
% spikeBeginTime = find(spikeSum>0,1,'first');
% spikeEndTime = find(spikeSum>0,1,'last');
% numBins = round((spikeEndTime-spikeBeginTime)*spike_bin_width/timeBinSize);

for j = 1:length(align_markers)    % For each desired marker
    align_offsets(:,j) = mod(find(events_in_bin == align_markers(j)),size(events_in_bin,1));  % Fast way to find offsets of each trial
    
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
    
    t_centers_old{j} = (psth_windows(j,1) : binSize_old : psth_windows(j,2) )';
    t_centers{j} = (psth_windows(j,1) + binSize/2 : stepSize : psth_windows(j,2) - binSize/2)'; % XCenters for PSTH. HH20141125

    % Same method as in HD: Anne 2014 NN, 10 ms bin + 50 ms Gaussian smoothing
    spike_hist{j} = zeros(sum(select_trials),length(t_centers{j})); % Preallocation
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

if mean(align_offsets(:,3)-align_offsets(:,2))<= 500 % Not memory saccade, do delayed saccade
    memSacOK = 0;
else
    memSacOK = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define temporal slices = [ beginT endT alignTo Marker color]
if memSacOK
    if mean(align_offsets(:,3)-align_offsets(:,2)) >= 1000  % Memory period = 1000
        temporal_Slice = {
            -500 -200 VSTIM_ON_CD 'Background' 'k';          % Background: -500 ~ -100 ms (04)
            75 400 VSTIM_ON_CD 'LS' 'r';               % LS: 75 ~ 275 ms (04);
            25 900 VSTIM_OFF_CD 'Memory'  'g';          % Sustained: 25 ~ 900 ms (05);
            -300 -50 SACCADE_BEGIN_CD 'pre-S' 'b';         %  pre-S: -200 ~ 0 (07);
            -50 50 SACCADE_BEGIN_CD 'S-Co' 'c';  % Saccade-coincident: -25 ~ 75 (07)
            100 500 SACCADE_BEGIN_CD 'post-S' 'm';     % post-S: 50 ~ 500 (07)
            };
        
    else  % Memory period = 500
        temporal_Slice = {
            -500 -200 VSTIM_ON_CD 'Background' 'k';          % Background: -500 ~ -100 ms (04)
            75 400 VSTIM_ON_CD 'LS' 'r';               % LS: 75 ~ 275 ms (04);
            25 500 VSTIM_OFF_CD 'Memory'  'g';          % Sustained: 25 ~ 500 ms (05);
            -300 -50 SACCADE_BEGIN_CD 'pre-S' 'b';         %  pre-S: -200 ~ 0 (07);
            -50 50 SACCADE_BEGIN_CD 'S-Co' 'c';  % Saccade-coincident: -25 ~ 75 (07)
            100 500 SACCADE_BEGIN_CD 'post-S' 'm';     % post-S: 50 ~ 500 (07)
            };
    end
else
    % Define temporal slices = [ beginT endT alignTo Marker color]
    if mean(align_offsets(:,2)-align_offsets(:,1)) >= 1000  % Visual > 1000
        temporal_Slice = {
            -500 -200 04 'Background' 'k';          % Background: -500 ~ -100 ms (04)
            75 400 VSTIM_ON_CD 'LS' 'r';               % LS: 75 ~ 275 ms (04);
            500 1400 VSTIM_ON_CD 'Sustained'  'g';          % Sustained: 200 ~ 800 ms (04);
            -300 -50 SACCADE_BEGIN_CD 'pre-S' 'b';         %  pre-S: -200 ~ 0 (07);
            -50 50 SACCADE_BEGIN_CD 'S-Co' 'c';  % Saccade-coincident: -25 ~ 75 (07)
            100 500 SACCADE_BEGIN_CD 'post-S' 'm';     % post-S: 50 ~ 500 (07)
            };
    else  % Sustained period = 500
        temporal_Slice = {
            -500 -200 04 'Background' 'k';          % Background: -500 ~ -100 ms (04)
            75 400 VSTIM_ON_CD 'LS' 'r';               % LS: 75 ~ 275 ms (04);
            500 900 VSTIM_ON_CD 'Sustained'  'g';          % Sustained: 200 ~ 800 ms (04);
            -300 -50 SACCADE_BEGIN_CD 'pre-S' 'b';         %  pre-S: -200 ~ 0 (07);
            -50 50 SACCADE_BEGIN_CD 'S-Co' 'c';  % Saccade-coincident: -25 ~ 75 (07)
            100 500 SACCADE_BEGIN_CD 'post-S' 'm';     % post-S: 50 ~ 500 (07)
            };
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preallocating
result_PSTH_old_mean = cell(length(align_markers),length(unique_heading));


% Now we align spike time to each marker, for different conditions
for i = 1:length(unique_heading)  % For each heading (delayed-saccade direction)
    select_binary = logical(heading_per_trial == unique_heading(i));
    select_trialN = find(select_binary == 1);
    
    spike_cache = NaN(spikeN_cum_trial(end),length(align_markers));       % Pre-allocation for speed
    spike_cache_numbers = 0;   % Number of spikes already in the cache
    
    for rep = 1: repetitionN  % For each repeats
        this_trialN = select_trialN(rep);
        
        % Put spike times of this trial into cache (with desired offsets)
        spike_cache_thisTrial  ...
            = repmat(spikeT_raw(spikeN_cum_trial(this_trialN)+1 : spikeN_cum_trial(this_trialN+1)),1,length(align_markers)) ...   % Raw spike time for this trial
            - repmat(align_offsets(this_trialN,:),spikeN_trial(this_trialN),1); % Minus time offsets for j'th marker

        spike_cache(spike_cache_numbers+1 : spike_cache_numbers+spikeN_trial(this_trialN) , :) = spike_cache_thisTrial;
        
        spike_cache_numbers = spike_cache_numbers + spikeN_trial(this_trialN);
        
        % Put spike counts during different phases into cache
        % Note this is not subjected to Gaussian smoothing
        for sliceN = 1:size(temporal_Slice,1)
            resp_trial{1,sliceN}(rep,i) = sum(spike_cache_thisTrial(:, align_markers==temporal_Slice{sliceN,3}) >= temporal_Slice{sliceN,1} & ...
                spike_cache_thisTrial(:,align_markers==temporal_Slice{sliceN,3}) <= temporal_Slice{sliceN,2}) ...
                / (temporal_Slice{sliceN,2}-temporal_Slice{sliceN,1}) * 1000; % in Hz
        end
    end
    
    % PSTH
    for j = 1:length(align_markers)    % For each desired marker
        % --- PSTH I used before (40 ms non-overlapping windows), only averaged ---
        % [result_PSTH{j,i}, x(:,j)] = hist(spike_cache(:,j),numBins);   % This would be wrong if the duration that have spikes are different for each trial!! HH20141125
        [result_PSTH_old_mean{j,i}] = hist(spike_cache(spike_cache(:,j) >= min(t_centers_old{j}) & spike_cache(:,j) <= max(t_centers_old{j}),j), t_centers_old{j});  % xCenters must be fixed! HH20141125
        result_PSTH_old_mean{j,i} = result_PSTH_old_mean{j,i} / repetitionN / binSize_old * 1000;  % in Hz
        
        % --- New Smoothed PSTH for each trial. HH20141207 ---
        % Same method as in HD: Anne 2014 NN, 10 ms bin + 50 ms Gaussian smoothing
        result_PSTH_anne_raw{j,i} = spike_hist{j}(select_binary,:);
        result_PSTH_anne_mean{j,i} = mean(result_PSTH_anne_raw{j,i}(~isnan(result_PSTH_anne_raw{j,i}(:,1)),:),1);
        
        % HH20180609
        result_PSTH_anne_sem{j,i} = nanstd(result_PSTH_anne_raw{j,i},[],1)./sqrt(size(result_PSTH_anne_raw{j,i},1));
    end
    
end


%% Interp PSTH for GROUP_GUI. @HH20150524
MemSac_interp_locations = 0:5:360;  % Resoluation = 5 degree
MemSac_original_locations = [unique_heading' unique_heading(1)+360]; % Make it circular

for j = 1:length(align_markers)
    MemSac_PSTH = reshape([result_PSTH_anne_mean{j,:}],length(t_centers{j}),[]); 
    MemSac_original_PSTH = MemSac_PSTH(:,[1:end 1]);
    
    MemSac_PSTH_sem = reshape([result_PSTH_anne_sem{j,:}],length(t_centers{j}),[]);
    MemSac_original_PSTH_sem = MemSac_PSTH_sem(:,[1:end 1]);
    
    % Interpolate original Memsac traces into higher spatial resolution
    MemSac_interp_PSTH{j} = nan(size(MemSac_original_PSTH,1),length(MemSac_interp_locations));
    MemSac_interp_PSTH_sem{j} = nan(size(MemSac_original_PSTH_sem,1),length(MemSac_interp_locations));
    
    for tttt = 1:size(MemSac_original_PSTH,1)
        MemSac_interp_PSTH{j}(tttt,:) = interp1(MemSac_original_locations, MemSac_original_PSTH(tttt,:), MemSac_interp_locations);
        MemSac_interp_PSTH_sem{j}(tttt,:) = interp1(MemSac_original_locations, MemSac_original_PSTH_sem(tttt,:), MemSac_interp_locations);
    end
end

%% Spatial tuning
vectSum = zeros(1,size(temporal_Slice,1));
vectAmp = zeros(1,size(temporal_Slice,1));

for sliceN = 1:size(temporal_Slice,1)
    % --- Calculating spatial selectivity ---
    % Statistics
    resp_mean{1,sliceN} = mean(resp_trial{sliceN},1);
    resp_err{1,sliceN} = std(resp_trial{sliceN},0,1) / sqrt(repetitionN);
    p(1,sliceN) =  anova1(resp_trial{sliceN},'','off');
    
    % --- Vector Sum ----
    [vectSum(sliceN), ~ , vectAmp(sliceN)] = vectorsumAngle(resp_mean{sliceN}, unique_heading);
    
    % --- DDI ----
    % HTI{sliceN} = vectAmp(sliceN) / sum(abs( resp(sliceN) - mean(resp{1})))
    DDI(1,sliceN) = ( max(resp_mean{sliceN})-min(resp_mean{sliceN}) ) / ( max(resp_mean{sliceN})-min(resp_mean{sliceN})+ ...
        2 * sqrt( sum(resp_err{sliceN}.^2 * repetitionN * (repetitionN-1))/(length(heading_per_trial) - length(unique_heading))));    % delta / (delta + 2 * sqrt(SSE/(N-M)))

    % --- PREF_NULL DI ---  @HH20150524
    [~,pref] = min(abs(vectSum(4) - MemSac_interp_locations)); % Use pre period
    pref = mod(pref-1,length(MemSac_interp_locations)-1)+1;
    null = mod(pref + (length(MemSac_interp_locations)-1)/2 -1, length(MemSac_interp_locations)-1)+1; % The opposite position
    
    j = align_markers == temporal_Slice{sliceN,3};
    DI_time_ind = temporal_Slice{sliceN,1} <= t_centers{j} & t_centers{j} <= temporal_Slice{sliceN,2};
    
    PREF_mean = mean(MemSac_interp_PSTH{j}(DI_time_ind,pref),1);
    NULL_mean = mean(MemSac_interp_PSTH{j}(DI_time_ind,null),1);
    
    PREF_NULL_DI(1,sliceN) = (PREF_mean - NULL_mean) / (PREF_mean + NULL_mean);

    
    % --- Calculating activity index (Barash, et al 1991) ---
    if sliceN == 1
        background = mean(resp_mean{1});
        activityIndex(1) = 0;
        continue; % No activity index for background
    end
    
    if background >= 1
        if max(resp_mean{sliceN}) >= 1.2 * background % Excitatory
            activityIndex(sliceN) = (max(resp_mean{sliceN})-background)/sqrt(background);
        else  % Inhibitory
            activityIndex(sliceN) = (min(resp_mean{sliceN})-background)/sqrt(background);
        end
    else
        if max(resp_mean{sliceN}) >= 1.2 * background % Excitatory
            activityIndex(sliceN) = (max(resp_mean{sliceN})-background);
        else  % Inhibitory
            activityIndex(sliceN) = (min(resp_mean{sliceN})-background);
        end
    end
end


% Plotting
positions = [8,7,4,1,2,3,6,9];
unique_heading_all = 0:45:359;

% PSTHs
for jj = 1:length(plot_markers)
    figure(4+jj-1); clf;
    if ~isempty(batch_flag)
        set(4+jj-1,'Position', [-1200,50 1000,700], 'Name', 'MemSac PSTH','color','w');    
    else
        set(4+jj-1,'Position', [50,50 1000,700], 'Name', 'MemSac PSTH','color','w');
    end
    
    orient landscape;
    set(0,'DefaultAxesXTickMode','auto');

    maxPSTH = 0;
    
    ax = tight_subplot(3,3,0.08,[0.05 0.05],[0.05 0.22]);

    for i = 1:length(unique_heading)
        axes(ax(positions(unique_heading(i)==unique_heading_all)));
        hold on;
        
        % Old-PSTH (40 ms)
        set(bar(t_centers_old{align_markers==plot_markers(jj)}, result_PSTH_old_mean{align_markers==plot_markers(jj),i},1),'FaceColor',[0.4 0.4 0.4],'EdgeColor','none');
        maxPSTH = max(maxPSTH, max(result_PSTH_old_mean{align_markers==plot_markers(jj),i}));
        
        % New-PSTH (10 ms bin + 50 ms Gaussian smooth)
        plot(t_centers{align_markers==plot_markers(jj)}, result_PSTH_anne_mean{align_markers==plot_markers(jj),i},'k','linewidth',2.5);
        
        if i>=6
            set(gca,'xticklabelmode','auto')
        end
        if i<=6 && i>=4
            set(gca,'yticklabelmode','auto')
        end
    end
    %     toc;
    
    axes(ax(1)); hold on;
    if plot_markers(jj) == 4 % Targ On
        title([FILE 'unit' num2str(SpikeChan) '  Target On (04)']);
        plotXLim = [-500 2000];
    elseif plot_markers(jj) == 7 % Sac On
        title([FILE 'unit' num2str(SpikeChan) '   Saccade On (07)']);
        plotXLim = [-2000 500];
    elseif plot_markers(jj) == 5 % Targ Off
        title([FILE 'unit' num2str(SpikeChan) '   Target Off (05)']);
        plotXLim = [-1000 1500];
    end
    
    for i = 1:length(unique_heading)
        axes(ax(positions(unique_heading(i)==unique_heading_all))); hold on;
        xlim(plotXLim);
        ylim([0, maxPSTH*1.1]);
        plot([0 0],[0 maxPSTH*1.1],'k-','LineWidth',2);
        
        % Annotating temporal slice windows
        for sliceN = 1:size(temporal_Slice,1)
            if temporal_Slice{sliceN,3} == plot_markers(jj) % Precise windows
                p_ = patch([temporal_Slice{sliceN,1} temporal_Slice{sliceN,1} temporal_Slice{sliceN,2} temporal_Slice{sliceN,2}],[maxPSTH*1.05 maxPSTH*1.1 maxPSTH*1.1 maxPSTH*1.05],temporal_Slice{sliceN,5});
                set(p_,'edgecolor','none');
                %                 set(p_,'FaceAlpha',0.4,'EdgeAlpha',0);
            else  % Windows that is not so precise
                % Mean shift
                meanShift = mean(align_offsets(:,align_markers == temporal_Slice{sliceN,3})-align_offsets(:,align_markers==plot_markers(jj)),1);
                plot( [meanShift meanShift],[0 maxPSTH*1.1],'k--','LineWidth',1);
                p_ = patch(meanShift+[temporal_Slice{sliceN,1} temporal_Slice{sliceN,1} temporal_Slice{sliceN,2} temporal_Slice{sliceN,2}],[maxPSTH*1.05 maxPSTH*1.1 maxPSTH*1.1 maxPSTH*1.05],temporal_Slice{sliceN,5});
                set(p_,'edgecolor','none');
                %                 set(p_,'FaceAlpha',0.4,'EdgeAlpha',0);
            end
            
        end
        
    end
    
end


% SetFigure
SetFigure(12)
% box off;
% set(gcf,'color','w');
% set(findall(gcf,'fontsize',10),'fontsize',12);
% set(findall(gcf,'tickdir','i'),'tickdir','o');
% set(findall(gcf,'type','axes','linewidth',0.5),'linewidth',1);


% Statistics
% figure(3); clf;
% set(3,'Position', [5,5 600,700], 'Name', 'PSTH','color','w');

% axis_plot = axes('Position',[0.05 0.05 0.9 0.6]);
% axis_typing = axes('Position',[0.05 0.66 0.9 0.27]);
axis_plot = subplot_tight(3,3,5,-0.1,[0.05 0.05 0.05 0.22]);
axis_typing = axes('Position',[0.8 0.66 0.4 0.27]);

% The largest goes first in the polar plot
[~, polarOrder] = sort(max(cell2mat(resp_mean'),[],2),1,'descend');

for sliceN = polarOrder(:)'
    
    axes(axis_plot);
    if strcmp(temporal_Slice{sliceN,5},''); continue; end
    
    if p(sliceN) < 0.05
        sigMarker = '-'; wid = 2;
    else
        sigMarker = '-'; wid = 1;
    end;
    
    polarwitherrorbar([unique_heading/180*pi; unique_heading(1)/180*pi], ...
        [resp_mean{sliceN}'; resp_mean{sliceN}(1)], (p(sliceN) < 0.05) * [resp_err{sliceN}' ; resp_err{sliceN}(1)],[temporal_Slice{sliceN,5} sigMarker],wid);
    hold on;
    
    if p(sliceN) < 0.05
        %         h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,[0 vectAmp(sliceN)],[temporal_Slice{sliceN,5} '-']);
        h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,[max(cell2mat(resp_mean))*0.9 max(cell2mat(resp_mean))*1.3],[temporal_Slice{sliceN,5} '-']);
        set(h_p,'LineWidth',3);
    end
    %
    axes(axis_typing); axis off;
    text(0,- sliceN * 0.15+1.1,sprintf('%s:  \\itp_ = %4.2g',temporal_Slice{sliceN,4},p(sliceN)),'color',temporal_Slice{sliceN,5},'FontSize',15 * (p(sliceN) < 0.05) + 9 * (p(sliceN) >= 0.05 || isnan(p(sliceN))));
    
    
end

% title(FILE);

% SetFigure
% box off;
% set(gcf,'color','w');
% set(findall(gcf,'fontsize',10),'fontsize',15);
% set(findall(gcf,'tickdir','i'),'tickdir','o');
% set(findall(gcf,'type','axes','linewidth',0.5),'linewidth',2);


%}


%% Data Saving

% Reorganized. HH20141124
config.batch_flag = batch_flag;

% Output information for test. HH20160415
if isempty(batch_flag)
    config.batch_flag = 'test.m';
    disp('Saving results to \batch\test\ ');
end

%%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = PackResult(FILE, SpikeChan, repetitionN, unique_stim_type, ... % Obligatory!!
                    p, vectSum, vectAmp, DDI, activityIndex, PREF_NULL_DI, ...
                    align_markers, plot_markers, binSize, temporal_Slice, unique_heading, ...
                    t_centers, align_offsets, result_PSTH_old_mean,  result_PSTH_anne_mean, result_PSTH_anne_raw, result_PSTH_anne_sem,...
                    MemSac_interp_locations, MemSac_interp_PSTH, MemSac_interp_PSTH_sem,...
                    resp_trial, resp_mean, resp_err);

config.suffix = 'MemSac';
config.xls_column_begin = 'p_bkgnd';
config.xls_column_end = 'AI_post';

% Figures to save
h_gcf = gcf;
config.save_figures = h_gcf.Number;

% Only once
config.sprint_once_marker = 'gggg';
config.sprint_once_contents = 'result.p(:), aziToHeading(result.vectSum), result.DDI(:),result. activityIndex(:)';
% Loop across each stim_type
config.sprint_loop_marker = {};
config.sprint_loop_contents = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveResult(config, result);
                


% %%%%%%%%%%%%%%%%%%%%%%%  Output   HH20140510 / HH20140621 %%%%%%%%%%%%%%%%%
% 
% if ~isempty(batch_flag)  % Figures
%     
%     %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     suffix = ['MemSac'];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     outpath = ['Z:\Data\Tempo\Batch\' batch_flag(1:end-2) '\'];
%    
%     if ~exist(outpath,'dir')
%         mkdir(outpath);    
%     end
%     
%     % Save figures
%     orient landscape;
%     savefilename = [outpath [FILE '_' num2str(SpikeChan)] '_' suffix '.png'];
% 
%     if exist(savefilename,'file')
%         delete(savefilename);
%     end
%     saveas(gcf,savefilename,'png');
%         
% end
% 
% 
% % Print results
% 
% %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Only once
% sprint_once_marker_temp = 'gggg';
% sprint_once_contents = 'p{:}, aziToHeading(cell2mat(vectSum)), DDI{:},activityIndex{:}';
% % Loop across each stim_type
% sprint_loop_marker_temp = '''''';  
% sprint_loop_contents = ''''''; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% sprint_once_marker = [];
% for i = 1:length(sprint_once_marker_temp)
%     sprint_once_marker = [sprint_once_marker '%' sprint_once_marker_temp(i) '\t '];
% end
% 
% sprint_loop_marker = [];
% for i = 1:length(sprint_loop_marker_temp)
%     sprint_loop_marker = [sprint_loop_marker '%' sprint_loop_marker_temp(i) '\t '];
% end
% 
% if ~isempty(batch_flag)  % Print to file
%     
%     outfile = [outpath suffix '.dat'];
%     printHead = 0;
%     if (exist(outfile, 'file') == 0)   % file does not yet exist
%         printHead = 1;
%     end
%     
%     fid = fopen(outfile, 'a');
%     % This line controls the output format
% 
%     if (printHead)
%         fprintf(fid, ['FILE\t ' sprint_once_contents '\t' sprint_loop_contents]);
%         fprintf(fid, '\r\n');
%     end
%     
%     fprintf(fid,'%s\t',[FILE '_' num2str(SpikeChan)]);
%     
% else  % Print to screen
%     fid = 1;  
% end
% 
% toClip = [];
% 
% % Print once
% if ~isempty(sprint_once_marker_temp)
%     eval(['buff = sprintf(sprint_once_marker,' sprint_once_contents ');']);
%     fprintf(fid, '%s', buff);
%     toClip = [toClip sprintf('%s', buff)];
% end
% 
% % Print loops
% for conditions = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
%     if sum(unique_stim_type == conditions)==0
%         buff = sprintf(sprint_loop_marker,ones(1,sum(sprint_loop_marker=='%'))*NaN);
%     else
%         k = find(unique_stim_type == conditions);
%         eval(['buff = sprintf(sprint_loop_marker,' sprint_loop_contents ');']);
%     end
%     fprintf(fid, '%s', buff);
%     toClip = [toClip sprintf('%s', buff)];
% end
% 
% fprintf(fid, '\r\n');
% toClip = [toClip sprintf('\r\n')];
% clipboard('copy',toClip);
% 
% if ~isempty(batch_flag)  % Close file
%     fclose(fid);
% end


% %%%%%%%%%%%%%%%%%%%%%%%  Batch Output.   %%%%%%%%%%%%%%%%%
% 
% if ~isempty(batch_flag)
%     
%     outpath = ['Z:\Data\Tempo\Batch\' batch_flag(1:end-2) '\'];
%    
%     if ~exist(outpath,'dir')
%         mkdir(outpath);    
%     end
%     
%     % Save figures
%     orient landscape;
%     
%     %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     suffix = ['MemSac'];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     savefilename = [outpath [FILE '_' num2str(SpikeChan)] '_' suffix '.png'];
%     if exist(savefilename)
%         delete(savefilename);
%     end
%     
%     saveas(4+jj-1,savefilename,'png');
    
    % Print results
    
%     %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sprint_txt_temp = 'ss';  % For each stim_type
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     sprint_txt = [];
%     for i = 1:length(sprint_txt_temp)
%         sprint_txt = [sprint_txt '%' sprint_txt_temp(i) '\t '];
%     end
%     
%     outfile = [outpath suffix '.dat'];
%     printHead = 0;
%     if (exist(outfile, 'file') == 0)   % file does not yet exist
%         printHead = 1;
%     end
%     
%     fid = fopen(outfile, 'a');
%     if (printHead)
%         fprintf(fid, 'FILE\t  reps, vest PREF, vest NULL, vis PREF, vis NULL, comb PREF, comb NULL ');
%         fprintf(fid, '\r\n');
%     end
%     
%     fprintf(fid,'%s\t %g\t %s\t',[FILE '_' num2str(SpikeChan)],repetitionN,num2str(result{1,1,2}.ts));
% 
%     for conditions = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
%         if sum(unique_stim_type == conditions)==0
%             buff = sprintf(sprint_txt, ones(1,length(sprint_txt_temp))*NaN);
%         else
%             k = find(unique_stim_type == conditions);
%             %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             buff = sprintf(sprint_txt, num2str(result{1,1,2}.ys((k-1)*2+1,:)), num2str(result{1,1,2}.ys((k-1)*2+2,:)));  % Fig = 1, row = 1, column(alignmarker) = 2
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
%         fprintf(fid, '%s', buff);
%     end
%     
%     fprintf(fid, '\r\n');
%     fclose(fid);
%     
% end;
%%
return;
