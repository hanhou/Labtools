% Rescue all available data from CED in case where .htb/.log files are lost.
% This is doable given the trial conditions are delivered to CED at the end of each trial.
% HH20150624

function good_data = RescueFromCED(PATH,FILE,rescue_info)
TEMPO_Defs;		% Reads in some definitions that we need

% --- Read spike2 data ---
if ~isempty(strfind(FILE,'.htb'))
    FILE = FILE(1:end-4);
end

fileName = [PATH FILE '.smr'];
spike2file=fopen(fileName);

[CHAN32, ~]=SONGetMarkerChannel(spike2file,32);
num_good_trials = sum(CHAN32.markers(:,1) == SUCCESS_CD);  % Number of good trials

% --- Initiation ---
% We will recover these data trial by trial
good_data.event_data = zeros(1, 5000, num_good_trials);
good_data.moog_params = nan(HEADING, num_good_trials, MOOG); % We only need STIM_TYPE (11) & HEADING (12) right now
good_data.misc_params = nan(OUTCOME, num_good_trials);

% We don't change the followings any more.
good_data.spike_data = zeros(3, 5000, num_good_trials);
good_data.lfp_data = [];
good_data.eye_data = nan(8, 1000, num_good_trials);

good_data.htb_header{EYE_DB}.skip = 4;
good_data.htb_header{EYE_DB}.speed_units = 100000;
good_data.htb_header{EYE_DB}.speed = 100;
good_data.htb_header{SPIKE_DB} = good_data.htb_header{EYE_DB};
good_data.htb_header{SPIKE_DB}.skip = 0;
good_data.htb_header{EVENT_DB} = good_data.htb_header{SPIKE_DB};
good_data.htb_header{EYE_DB}.nchannels = 8;

% --- Recover trial condition list ---
stim_type = rescue_info{2};
heading = rescue_info{3};

if strcmp(rescue_info{1},'HD')
    good_data.one_time_params(PROTOCOL) = HEADING_DISCRIM;
    
    condition_list = [reshape(repmat(stim_type,sum(heading <=0),1),[],1) repmat(heading(heading <=0)',length(stim_type),1)];
    condition_list = [condition_list; [reshape(repmat(fliplr(stim_type),sum(heading > 0),1),[],1) ...
                                       repmat(heading(heading >0)',length(stim_type),1)]];
elseif strcmp(rescue_info{1},'MemSac')
    good_data.one_time_params(PROTOCOL) = MEMORY_SACCADE;
    condition_list = [repmat(stim_type,length(heading),1) heading'];
elseif strcmp(rescue_info{1},'DelSac')
    good_data.one_time_params(PROTOCOL) = DELAYED_SACCADE;
    condition_list = [repmat(stim_type,length(heading),1) heading'];
end
    

% --- Recover data for each trial ---
marker = CHAN32.markers(:,1);
marker_t = CHAN32.timings * 1000; % in ms 

marker_success = find(marker == SUCCESS_CD);
marker_start = find(marker == TRIAL_START_CD);
marker_end = find(marker == TRIAL_END_CD);
marker_stim_on = find(marker == VSTIM_ON_CD);

for i = 1:num_good_trials
    
    % --- Fill good_data.event_data ---
    this_success_ind = marker_success(i);
    this_start_ind = marker_start(find(this_success_ind > marker_start,1,'last')); % The nearest TRIAL_START_CD before this_success
    this_end_ind = marker_end(find(this_success_ind < marker_end,1,'first'));  % The end of this trial
    this_stim_on_ind = marker_stim_on(find(this_success_ind > marker_stim_on,1,'last')); % The nearest TRIAL_START_CD before this_success

    this_marker_to_fill = marker(this_start_ind:this_end_ind);
    this_marker_t = round(marker_t(this_start_ind:this_end_ind) - marker_t(this_stim_on_ind) + 991); % VSTIM_ON fixed at 991 (as in other normal .htb)
    
    positive_this_marker_t = this_marker_t > 0;
    if sum(positive_this_marker_t) < length(this_marker_t)
        fprintf('Neglect markers in #%g trial\n',i);
    end
    
    % Fill event_data
    good_data.event_data(1,this_marker_t(positive_this_marker_t),i) = this_marker_to_fill(positive_this_marker_t);    
    
    % --- Recover trial condition ---
    this_trial_condition = this_marker_to_fill(this_marker_to_fill >= 48) - 47; % Trial condition
    good_data.moog_params([STIM_TYPE HEADING COHERENCE], i, MOOG) = [condition_list(this_trial_condition,:) rescue_info{4}];
    
    % --- Trial outcome ---
    if sum(REWARD_CD == this_marker_to_fill) > 0
        this_outcome = CORRECT;
    else
        this_outcome = ERR_WRONG_CHOICE;
    end
    good_data.misc_params(OUTCOME,i) = this_outcome;
    
end

