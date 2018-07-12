%-----------------------------------------------------------------------------------------------------------------------
% PSTH for choice-related activity
% Last modified HH20140526
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cum_PSTH_HH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);

if_figure = isempty(batch_flag);

TEMPO_Defs;
Path_Defs;


%% Commented by HH20140523
%
%% Added by HH20140523
%%{

%% Get data

%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%
method = 0; % 0: Maximum likelihood; 1: Square error
tolerance = 10;

% Define align markers and offsets
align_markers = {
    % Marker    Before(ms)    After(ms)    Notes    Which other markers are presented?   Their notes
    VSTIM_ON_CD,  -500, 2200, 'Stim On' , [VSTIM_OFF_CD SACCADE_BEGIN_CD], {'Stim Off','Sac On'};
    SACCADE_BEGIN_CD,  -1000, 500, 'Saccade On' , [VSTIM_ON_CD VSTIM_OFF_CD], {'Stim On', 'Stim Off'};
    %     VSTIM_ON_CD,  -500, 2200, 'Stim On' , [VSTIM_OFF_CD SACCADE_BEGIN_CD], {'Stim Off','Sac On'};
    %     SACCADE_BEGIN_CD,  -500, 700, 'Saccade On' , [VSTIM_ON_CD VSTIM_OFF_CD], {'Stim On', 'Stim Off'};
    };

% Permutation numbers
CD_perm_N = -1;   % Light version
CP_perm_N = -1;   % Light version
ChoicePref_perm_N = 1000;
ModalityPref_perm_N = ChoicePref_perm_N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;


% Override the default eye channel settings. HH20150722
if data.one_time_params(LEFT_EYE_X_CHANNEL) > 0 % NOT (Is NaN (not overriden) or is zero (rescued from CED))
    LEYE_H = data.one_time_params(LEFT_EYE_X_CHANNEL);
end

if data.one_time_params(LEFT_EYE_Y_CHANNEL) > 0
    LEYE_V = data.one_time_params(LEFT_EYE_Y_CHANNEL);
end

if data.one_time_params(RIGHT_EYE_X_CHANNEL) > 0
    REYE_H = data.one_time_params(RIGHT_EYE_X_CHANNEL);
end

if data.one_time_params(RIGHT_EYE_Y_CHANNEL) > 0
    REYE_V = data.one_time_params(RIGHT_EYE_Y_CHANNEL);
end


stim_type_names = {'All','Vest','Vis','Comb'}; % stim_type = 0, 1, 2, 3
stim_type_colors = [0 0 1; 1 0 0; 0 0.8 0.4];

% -- Trial information

trials = 1:size(data.moog_params,2);		% a vector of trial indices

% If length(BegTrial) > 1 and all elements are positive, they are trials to be included.
% Else, if all elements are negative, they are trials to be excluded.
% This enable us to exclude certain trials ** DURING ** the recording more easily. HH20150410
select_trials = false(length(trials),1);
if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
    select_trials(BegTrial:EndTrial) = true;
elseif all(BegTrial > 0) % To be included
    select_trials(BegTrial) = true;
elseif all(BegTrial < 0) % To be excluded
    select_trials(-BegTrial) = true;
    select_trials = ~ select_trials;
else
    disp('Trial selection error...');
    error(sprintf('Trial selection error! You chose: %g, %g',BegTrial,EndTrial));
    % keyboard;
end

stim_type_per_trial = data.moog_params(STIM_TYPE,select_trials,MOOG);
heading_per_trial   = data.moog_params(HEADING, select_trials, MOOG);
outcome_per_trial = data.misc_params(OUTCOME, select_trials)';

unique_stim_type = munique(stim_type_per_trial');
unique_heading = munique(heading_per_trial');

repetitionN = floor(sum(select_trials) / length(unique_heading) / length(unique_stim_type)) ;

% -- Time information
eye_timeWin = 1000/(data.htb_header{EYE_DB}.speed_units/data.htb_header{EYE_DB}.speed/(data.htb_header{EYE_DB}.skip+1)); % in ms
spike_timeWin = 1000/(data.htb_header{SPIKE_DB}.speed_units/data.htb_header{SPIKE_DB}.speed/(data.htb_header{SPIKE_DB}.skip+1)); % in ms
event_timeWin = 1000/(data.htb_header{EVENT_DB}.speed_units/data.htb_header{EVENT_DB}.speed/(data.htb_header{EVENT_DB}.skip+1)); % in ms

% -- Spike data
spike_in_bin = squeeze(data.spike_data(SpikeChan,:,select_trials))';   % TrialNum * 5000
spike_in_bin( spike_in_bin > 100 ) = 1; % something is absolutely wrong

% -- Event data
event_in_bin = squeeze(data.event_data(:,:,select_trials))';  % TrialNum * 5000

% -- Eye data --
% Fix the bug that if the first selected trial > 1, the eye traces got messed up. @HH20160906
eye_data = data.eye_data(:,:,select_trials);

% Monkey's choice
LEFT = 1;
RIGHT = 2;

% The previous one was awful. HH20140522
choice_per_trial = LEFT * squeeze(sum(event_in_bin == IN_T2_WIN_CD,2)) + RIGHT * squeeze(sum(event_in_bin == IN_T1_WIN_CD,2));
if length(unique(choice_per_trial)) > 2  % This is safer
    disp('Neither T1 or T2 chosen / More than one target chosen.  This should not happen! File must be bogus.');
    
    %%%%%%%%%%%%  CAUTION HERE !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(FILE,'m5c174r1') % This file is ugly (eye trace noise!! All LEFT actually)  HH20141224
        for ii = find(choice_per_trial==3)'
            %             figure(ii);
            %             plot(data.eye_data(1,:,ii)); hold on;
            m7 = find(event_in_bin(ii,:)==7)/5;
            m8 = find(event_in_bin(ii,:)==8)/5;
            m9 = find(event_in_bin(ii,:)==9)/5;
            %             plot([m7 m7],[-20 20],'r',[m8 m8],[-20 20],'k',[m9 m9],[-20 20],'k');
            disp(find(event_in_bin(ii,:)==8 |event_in_bin(ii,:)==9))
            
            choice_per_trial(ii) = LEFT;            % Overwrite all LEFT !!  HH20141224
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        keyboard; beep;
    end
    
end

%% 1. Align data


% Save configuration
% result.align_markers = align_markers;

for j = 1:size(align_markers,1)    % For each desired marker
    align_offsets(:,j) = mod(find(event_in_bin' == align_markers{j,1}),size(event_in_bin,2));  % Fast way to find offsets for each trial
    
    spike_aligned{1,j} = zeros(sum(select_trials), ceil((align_markers{j,3} - align_markers{j,2}) / spike_timeWin) + 1); % allocation
    spike_aligned{2,j} = align_markers{j,2} : spike_timeWin : align_markers{j,3} ;% Time stamps
    
    eyeX_aligned{1,j} = zeros(sum(select_trials), ceil((align_markers{j,3} - align_markers{j,2}) / eye_timeWin) + 1); % Preallocation
    eyeX_aligned{2,j} = align_markers{j,2} : eye_timeWin : align_markers{j,3} ; % Time stamps
    eyeY_aligned{1,j} = zeros(sum(select_trials), ceil((align_markers{j,3} - align_markers{j,2}) / eye_timeWin) + 1); % Preallocation
    eyeY_aligned{2,j} = align_markers{j,2} : eye_timeWin : align_markers{j,3} ; % Time stamps
    
    for i = 1:sum(select_trials)
        winBeg = align_offsets(i,j) + ceil(align_markers{j,2} / spike_timeWin);
        winEnd = align_offsets(i,j) + ceil(align_markers{j,3} / spike_timeWin);
        spike_aligned{1,j}(i,:) = spike_in_bin (i, winBeg : winEnd);   % Align each trial
        
        % Note that sampling rate of eye traces (200 Hz) is lower than spike data (1000 Hz), so the alignment would be not perfect.
        winBeg_eye = round(align_offsets(i,j) / eye_timeWin) + ceil(align_markers{j,2} / eye_timeWin);
        
        % -- The below two lines are wrong because I failed to select trials for eye data!! @HH20160906
        % eyeX_aligned{1,j}(i,:) = squeeze(data.eye_data (LEYE_H, winBeg_eye : winBeg_eye + size(eyeX_aligned{1,j},2) -1 , i));   % Align each trial
        % eyeY_aligned{1,j}(i,:) = squeeze(data.eye_data (LEYE_V, winBeg_eye : winBeg_eye + size(eyeX_aligned{2,j},2) -1, i));   % Align each trial
        eyeX_aligned{1,j}(i,:) = squeeze(eye_data (LEYE_H, winBeg_eye : winBeg_eye + size(eyeX_aligned{1,j},2) -1 , i));   % Align each trial
        eyeY_aligned{1,j}(i,:) = squeeze(eye_data (LEYE_V, winBeg_eye : winBeg_eye + size(eyeX_aligned{2,j},2) -1, i));   % Align each trial
        
        % Calibrate eye offset (average 200 ms after 04, then we calculate eye_offset when j = 1)
        if j == 1
            eye_offsetX(i) = mean(eyeX_aligned{1,1}(i,winBeg_eye : winBeg_eye + 200/eye_timeWin));
            eye_offsetY(i) = mean(eyeY_aligned{1,1}(i,winBeg_eye : winBeg_eye + 200/eye_timeWin));
        end
        
        eyeX_aligned{1,j}(i,:) = eyeX_aligned{1,j}(i,:) - eye_offsetX(i);
        eyeY_aligned{1,j}(i,:) = eyeY_aligned{1,j}(i,:) - eye_offsetY(i);
    end
    
    % Other time markers that are presented
    for jj = 1:length(align_markers{j,5})
        align_offsets_others{j}(:,jj) = mod(find(event_in_bin' == align_markers{j,5}(jj)),size(event_in_bin,2))-align_offsets(:,j);  % Fast way to find offsets for each trial
    end
end

%% 1.1 Preferred Direction
%%%%%%%%  Define "Preferred direction" %%%%%%%%%%%
% Use a global "preferred direction" for all conditions.
% But see pref. direction for CP below.  HH20141126  < XXX @HH20150417>

% To prevent from getting a weird result that three conditions have different preferred directions,
% now I decide to use the "global preferred direction" for CD and CP as well @HH20150417
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% left_all = sum(sum(spike_aligned{1,2}(choice_per_trial == LEFT, 1: round(abs(align_markers{2,2})/(align_markers{2,3} - align_markers{2,2})*end))))/sum(choice_per_trial == LEFT); % Roughly pre-saccade time widow
% right_all = sum(sum(spike_aligned{1,2}(choice_per_trial == RIGHT, 1: round(abs(align_markers{2,2})/(align_markers{2,3} - align_markers{2,2})*end))))/sum(choice_per_trial == RIGHT); % Roughly pre-saccade time widow
left_all = sum(mean(spike_aligned{1,2}( choice_per_trial == LEFT & outcome_per_trial == CORRECT, -500 <= spike_aligned{2,2} &  spike_aligned{2,2} <= -50))); % Roughly pre-saccade time widow
right_all = sum(mean(spike_aligned{1,2}( choice_per_trial == RIGHT & outcome_per_trial == CORRECT, -500 <= spike_aligned{2,2} &  spike_aligned{2,2} <= -50 ))); % Roughly pre-saccade time widow

% Debugging.  HH20150125
% figure();  hold on;
% plot(t_centers{2},mean(spike_hist{2}(...
%     choice_per_trial == LEFT & outcome_per_trial == CORRECT & stim_type_per_trial' >= 1,:)));
% plot(t_centers{2},mean(spike_hist{2}(...
%     choice_per_trial == RIGHT & outcome_per_trial == CORRECT & stim_type_per_trial' >= 1,:)),'r');
%
% figure();  hold on;
% plot(t_centers{2},mean(spike_hist{2}(...
%     choice_per_trial == LEFT & outcome_per_trial == ERR_WRONG_CHOICE & stim_type_per_trial' <= 3,:)),'Linew',2);
% plot(t_centers{2},mean(spike_hist{2}(...
%     choice_per_trial == RIGHT & outcome_per_trial == ERR_WRONG_CHOICE & stim_type_per_trial' <= 3,:)),'r','Linew',2);


if left_all >= right_all
    PREF = LEFT;
    NULL = RIGHT;
    choose_names = {'Chs PREF(L)' ;'Chs NULL(R)'};
else
    PREF = RIGHT;
    NULL = LEFT;
    choose_names = {'Chs PREF(R)' ;'Chs NULL(L)'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------ Override of PREF direction -----------------
% In the following cells, I manually switch the PREF and NULL directions, because:
%   1. PREF determined above is opposite to their PREF in MemSac task (which is unreasonable)
%   2. Their ramping activities are noisy per se (even PREFs of different modalities do not match), which makes PREF incorrect
%   3. Examining the preference of post-sac activity further justifies my manual switch
% This is a temporary workaround (the deadline of my progress report is coming man!!)
% Maybe someday I will rewrite it in a more formal (soft-coded) way.
% @HH20150524

flip_PREF_cell = {
    % File name, SpikeChan
    'm5c211r3_4_5',5;
    'm5c145r1',5};

for ff = 1:size(flip_PREF_cell)
    if strcmp(flip_PREF_cell{ff,1},FILE) && SpikeChan == flip_PREF_cell{ff,2}
        PREF = NULL;
        NULL = setdiff([LEFT RIGHT],PREF);
        beep;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.2 Get actual target location for better MemSac measurement. @HH20150524

j = 2; % Saccade onset
eye_target_period = 50 <= eyeX_aligned{2,j} & eyeX_aligned{2,j} <= 150 ; % The period when eyes are at targets

% Eye positions when monkey chooses PREF
eyeXs_pref = mean(eyeX_aligned{1,j}(choice_per_trial == PREF, eye_target_period),2);
eyeYs_pref = mean(eyeY_aligned{1,j}(choice_per_trial == PREF, eye_target_period),2);

% Get mean location of PREF_target
% This would let PREF_target_location be [-90,270] ... HH20160906
PREF_target_location = median(atan(eyeYs_pref./eyeXs_pref)/pi * 180 + (eyeXs_pref < 0)*180);


%% 1.5 Raster plot
if if_figure
    
    set(figure(59),'pos',[90 563 1525 398],'name','Raster plot'); clf;
    
    j = 1;
    j_other = 2;
    num_trials_for_raster = 20;
    
    face_col = {'none','k'};
    spike_height = 0.3;
    
    for stim_type = 1:3
        find_trial_pref = find((stim_type_per_trial' == stim_type) & (choice_per_trial == PREF) & (outcome_per_trial == CORRECT));
        find_trial_null = find((stim_type_per_trial' == stim_type) & (choice_per_trial == NULL) & (outcome_per_trial == CORRECT));
        
        if length(find_trial_pref) >= num_trials_for_raster && length(find_trial_null) >= num_trials_for_raster
            
            % Find and sort the "other marker" (saccade onset for j = 1)
            other_time_marker_pref = align_offsets_others{j}(find_trial_pref,j_other);
            other_time_marker_null = align_offsets_others{j}(find_trial_null,j_other);
            
            [~,sort_pref] = sort(abs(other_time_marker_pref));
            [~,sort_null] = sort(abs(other_time_marker_null));
            
            % Choose who are gonna to be plotted
            select_for_raster_pref = sort_pref(round(linspace(round(length(sort_pref)*0.05),round(length(sort_pref)*0.95),num_trials_for_raster)));
            %         select_for_raster_null = sort_null(round(linspace(round(length(sort_pref)*0.05),round(length(sort_pref)*0.95),num_trials_for_raster)));   % Bug found by Yuchen 20160321
            select_for_raster_null = sort_null(round(linspace(round(length(sort_null)*0.05),round(length(sort_null)*0.95),num_trials_for_raster)));
            
            other_time_marker_all = [other_time_marker_null(select_for_raster_null); other_time_marker_pref(select_for_raster_pref)];
            raster_spike_time_all = spike_aligned{1,j}([find_trial_null(select_for_raster_null);find_trial_pref(select_for_raster_pref)],:);
            raster_t = spike_aligned{2,j}(1,:);
            
            % Plotting
            subplot(1,3,stim_type); hold on
            for tr = 1:2*num_trials_for_raster
                if sum(raster_spike_time_all(tr,:))>0
                    plot([raster_t(logical(raster_spike_time_all(tr,:))); raster_t(logical(raster_spike_time_all(tr,:)))],...
                        [tr-spike_height tr+spike_height],'-','color',stim_type_colors(stim_type,:),'linew',1.2);
                    plot(other_time_marker_all(tr),tr + spike_height,'vk','markerfacecol',face_col{(tr > num_trials_for_raster) + 1},'markersize',5);
                end
            end
            set(gca,'ytick',[]);
            ylim([0 2*num_trials_for_raster+1]); xlim([-200 2500]);
            plot([0 0],ylim,'k-'); plot([1500 1500],ylim,'k-'); plot(xlim,0.5 + [num_trials_for_raster num_trials_for_raster],'k--');
        end
    end
    SetFigure(15);
end

%% 2. Calculate PSTH using sliding windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time windows
% binSize_rate = 80  % in ms
% stepSize_rate = 20; % in ms
% smoothFactor = 1; % in bin

% Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
binSize_rate = 10;  % in ms
stepSize_rate = 10; % in ms
smoothFactor = 50; % in ms !!
p_critical = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% result.binSize_rate = binSize_rate;
% result.stepSize_rate = stepSize_rate;


for j = 1:size(align_markers,1)    % For each desired marker
    
    t_centers{j} = align_markers{j,2} + binSize_rate/2 : stepSize_rate : align_markers{j,3} - binSize_rate/2; % Centers of PSTH time windows
    
    spike_hist{j} = zeros(sum(select_trials),length(t_centers{j})); % Preallocation
    for k = 1:length(t_centers{j})
        winBeg = ceil(((k-1) * stepSize_rate) / spike_timeWin) + 1;
        winEnd = ceil(((k-1) * stepSize_rate + binSize_rate) / spike_timeWin) + 1;
        spike_hist{j}(:,k) = sum(spike_aligned{1,j}(: , winBeg:winEnd),2) / binSize_rate*1000 ;  % in Hz
    end
    
    % --- Smoothing ---
    % Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
    % Note that this only influence PSTH calculation (spike_hist), not CP
    if smoothFactor > 0
        for i = 1:size(spike_hist{j},1) % Each trial
            spike_hist{j}(i,:) =  GaussSmooth(t_centers{j},spike_hist{j}(i,:),smoothFactor);
        end
    end
end


%% 3. Sort trials into different categories

% Define colormap for headings
colors = colormap(cool); % This is really cool.
color_for_headings = colors(end-round(linspace(1,64,sum(unique_heading<0)))+1,:);
if sum(unique_heading ==0)
    color_for_headings = [color_for_headings; 0 0 0; flipud(color_for_headings)]; % Zero heading: black
else
    color_for_headings = [color_for_headings; flipud(color_for_headings)];
end

style_for_headings = vertcat(repmat({'-'},sum(unique_heading<=0),1), repmat({'--'},sum(unique_heading>0),1)); % Preferred LEFT by default
if PREF == RIGHT
    style_for_headings = flipud(style_for_headings);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define sort information (The syntax is a little bit intricate, but when you understand it, you'll love it.)
VAR = 1;
LOGLS = 2;
VALUES = 3;
COLORS = 4;
LINESTYLES = 5;
NOTES = 6;
ERRORBAR = 7;

sort_info = {
    % Each cell for differnt bases of classification (generate separated figures)
    
    % Fig.0: Rows: stim types; Sort according to: Choose Right and Coose Left, all headings
    %     {
    %         % Rows in subplot (stim types)
    %         {   % stim_types (0 = all)   % Notes
    %              0, stim_type_names(0+1)
    %         };
    %
    %         % In each subplot (max nested level: 2)
    %         %  Variable(s) , Logical,  Values  , Colors,  LineStyles,    Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    %         {
    %             'choice_per_trial', {'=='}, [PREF NULL], [1 0 0; 0 0 0], {'-'; '-'}, choose_names, 2;
    %         }
    %     };
    
    % Fig.1: Rows: none; Sort according to: stim type & Choice
    {
    % Rows in subplot (stim types)
    {   % stim_types (0 = all)   % Notes
    0, stim_type_names(0+1)
    };
    
    % In each subplot (max nested level: 2)
    {    %  Variable(s) , Logical,  Values  , Colors,  LineStyles,    Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    'stim_type_per_trial''', {'=='}, unique_stim_type, stim_type_colors(unique_stim_type,:) ,{},stim_type_names(unique_stim_type+1), 999;
    'choice_per_trial', {'=='}, [PREF NULL], [], {'-'; '--'}, choose_names, 2;
    }
    }
    
    % Fig.2: Rows: stim types; Sort according to: Abs(Heading angles) && choices
    % This is equivalent to conventional "sorted by angles but only correct trials"
    {
    % Rows in subplot (stim types)
    {   % stim_types (0 = all)   % Notes
    unique_stim_type, stim_type_names(unique_stim_type+1)
    };
    
    % In each subplot (max nested level: 2)
    {    %  Variable(s) ,  Logical,  Values  , Colors,  LineStyles,    Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    'abs(heading_per_trial)''',{'=='}, unique(abs(unique_heading)), color_for_headings(ceil(end/2):end,:), {}, cellstr(num2str(unique(abs(unique_heading)))) , 999;
    'choice_per_trial', {'=='}, [PREF NULL], [], {'-'; '--'}, choose_names, 0;
    }
    }
    
    %     % Fig.3: Rows: stim types; Sort according to: Large/small angles && choices (difficulty)
    {
    % Rows in subplot (stim types)
    {   % stim_types (0 = all)   % Notes
    unique_stim_type, stim_type_names(unique_stim_type+1)
    };
    
    % In each subplot (max nested level: 2)
    {    %  Variable(s) ,  Logical,  Values  , Colors,  LineStyles,    Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    'abs(heading_per_trial)''',{'<','>='}, [abs(unique_heading(fix(length(unique_heading)/2)/2)) abs(unique_heading(fix(length(unique_heading)/2)/2))], [0 0 0; 1 0 0], {}, {['|Heading| < ' num2str(abs(unique_heading(fix(length(unique_heading)/2)/2)))],['|Heading| >= ' num2str(abs(unique_heading(fix(length(unique_heading)/2)/2)))]} , 999;
    'choice_per_trial', {'=='}, [PREF NULL], [], {'-'; '--'}, choose_names, 2;
    }
    }
    
    % Fig.4: Rows: stim types; Sort according to: outcomes
    {
    % Rows in subplot (stim types)
    {   % stim_types (0 = all)   % Notes
    unique_stim_type, stim_type_names(unique_stim_type+1)
    };
    
    % In each subplot (max nested level: 2)
    {    %  Variable(s) ,  Logical,  Values  , Colors,  LineStyles,    Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    'outcome_per_trial',{'=='}, [CORRECT ERR_WRONG_CHOICE], [1 0 0; 0 0 0], {'-';'-'}, {'Correct','Wrong'} , 999;
    'choice_per_trial', {'=='}, [PREF NULL], [], {'-'; '--'}, choose_names, 2;
    }
    }
    
    % Fig.5: Rows: stim types; Sort according to: Different angles but wrong trials. @HH20150523
    {
    % Rows in subplot (stim types)
    {   % stim_types (0 = all)   % Notes
    unique_stim_type, stim_type_names(unique_stim_type+1)
    };
    
    % In each subplot (max nested level: 2)
    {    %  Variable(s) ,  Logical,  Values  , Colors,  LineStyles,    Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    'outcome_per_trial',{'=='}, [ERR_WRONG_CHOICE], [0 0 0], {'-';'-'}, {'Wrong'} , 999;
    'heading_per_trial''',{'=='}, unique_heading(unique_heading~=0), color_for_headings, style_for_headings(unique_heading~=0), cellstr(num2str(unique_heading(unique_heading~=0))) , 0;
    }
    }
    
    % Fig.6: Rows: stim types; Sort according to: Different angles but ALL (correct + wrong) trials. @HH20180711 (For GuYong's Hangzhou ppt)
    {
    % Rows in subplot (stim types)
    {   % stim_types (0 = all)   % Notes
    unique_stim_type, stim_type_names(unique_stim_type+1)
    };
    
    % In each subplot (max nested level: 2)
    {    %  Variable(s) ,  Logical,  Values  , Colors,  LineStyles,    Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    'heading_per_trial''',{'=='}, unique_heading, color_for_headings, style_for_headings, cellstr(num2str(unique_heading)) , 0;
    }
    }

    
    % %     % Fig.5: Rows: stim types; Sort according to: Heading angles
    %     {
    %         % Rows in subplot (stim types)
    %         {   % stim_types (0 = all)   % Notes
    %             unique_stim_type, stim_type_names(unique_stim_type+1)
    %         };
    %
    %         % In each subplot (max nested level: 2)
    %         {    %  Variable(s) ,  Logical,  Values  , Colors,  LineStyles,  Notes,   Errorbar? (0:Nothing; 1:95% CI; 2: p<0.05; 3:both)
    % %             'heading_per_trial''', {'=='}, unique_heading, color_for_headings, style_for_headings, cellstr(num2str(unique_heading)) , 0;
    %             'heading_per_trial''', {'<' '>='}, [0 0], [0 0 0; 0 0 1], {'-';'-'}, {'Heading < 0','Heading >=0'} , 2;
    %         }
    %     };
    
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Organize data...\n')
% Save configuration
% result.sort_info = sort_info;

% Calculate temporal duration ratio for each align_marker to keep the plots
% in scale

temp_duration_ratio = zeros(1,size(align_markers,1));
for j = 1:size(align_markers,1)
    temp_duration_ratio(j) = align_markers{j,3} - align_markers{j,2};
end

% Plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transparent = 1;
outcome_mask = CORRECT;
% outcome_mask = ERR_WRONG_CHOICE;
outcome_mask_enable = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save configuration
% result.smoothFactor = smoothFactor;
% result.outcome_enable = outcome_enable;

for sortInd = 1:length(sort_info) % For each figure
    
    stim_type_to_plot = sort_info{sortInd}{1}{1}; % [unique_stim_type 4]
    
    if if_figure
        set(figure(60+sortInd),'position',[187-1300*~isempty(batch_flag)  91   1100   600]); clf;
        set(61,'position',[-1110 80 1100 600]);
        h_subplot = tight_subplot(length(stim_type_to_plot),size(align_markers,1),[0.01 0.02],[0.1 0.1],[0.08 0.03],[],temp_duration_ratio);
    end
    
    nest_levels = size(sort_info{sortInd}{2},1);
    
    % Preallocation of legend
    if nest_levels > 2 ;
        disp('Too many nest levels...'); keyboard;
        return;
    elseif nest_levels == 1
        catsOut = sort_info{sortInd}{2}{1,VALUES};
        logOut = sort_info{sortInd}{2}{1,LOGLS};
        if length(logOut) == 1; logOut = repmat(logOut,1,length(catsOut)); end
        
        h_legend = zeros(length(catsOut),1);
        txt_legend = cell(size(h_legend));
    elseif nest_levels == 2
        catsOut = sort_info{sortInd}{2}{1,VALUES};
        logOut = sort_info{sortInd}{2}{1,LOGLS};
        if length(logOut) == 1; logOut = repmat(logOut,1,length(catsOut)); end
        
        catsIn = sort_info{sortInd}{2}{2,VALUES};
        logIn = sort_info{sortInd}{2}{2,LOGLS};
        if length(logIn) == 1; logIn = repmat(logIn,1,length(catsIn)); end
        
        h_legend = zeros(length(catsOut)*length(catsIn),1);
        txt_legend = cell(size(h_legend));
    end
    
    for k = 1:length(stim_type_to_plot)  % The last one: all conditions
        for j = 1:size(align_markers,1)
            
            if if_figure
                l = sub2ind([length(stim_type_to_plot) size(align_markers,1)],k,j);
                %             axes(h_subplot(l));
                set(gcf,'CurrentAxes',h_subplot(l));
            end
            
            if nest_levels == 1 % 1 nested condition
                
                for catNum_Out = 1:length(catsOut)
                    catNum_In = 1;
                    
                    selected_condition = eval(sprintf('%s %s %s',sort_info{sortInd}{2}{1,VAR},logOut{catNum_Out}, num2str(catsOut(catNum_Out))));
                    
                    if sort_info{sortInd}{1}{1}(k) > 0  % Not all conditions
                        selected_condition = selected_condition & (stim_type_per_trial' == unique_stim_type(k));
                    end
                    
                    % If the sorting type is not outcome per se (or "ALL" override @HH20180711), we only
                    % choose the correct trials by default, except 0 heading.
                    if ~ (strcmp(sort_info{sortInd}{2}{1},'outcome_per_trial') || sortInd == 6)
                        selected_condition = selected_condition & (~outcome_mask_enable | outcome_per_trial == outcome_mask | heading_per_trial' == 0); %% HH20140825
                    end
                    
                    lineColor = sort_info{sortInd}{2}{1,COLORS}(catNum_Out,:);
                    lineStyle = sort_info{sortInd}{2}{1,LINESTYLES}{catNum_Out,:};
                    
                    
                    % ------- Data smoothing has been moved to the previous section --------- HH20141207
                    % ys = smooth(mean(spike_hist{j}(selected_condition,:),1),smoothFactor);
                    
                    % Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
                    % ys = GaussSmooth(t_centers{j},mean(spike_hist{j}(selected_condition,:),1),smoothFactor);
                    
                    ys = mean(spike_hist{j}(selected_condition,:),1);
                    sem = std(spike_hist{j}(selected_condition,:),0,1)/sqrt(sum(selected_condition)); % @HH20160906
                    
                    % Save data
                    PSTH{j,sortInd,k}.raw{catNum_Out,1} = spike_hist{j}(selected_condition,:); % Save all raw data for further processing. HH20141118
                    
                    PSTH{j,sortInd,k}.ys(catNum_Out,:) = ys;
                    PSTH{j,sortInd,k}.sem(catNum_Out,:) = sem; % @HH20160906
                    PSTH{j,sortInd,k}.ts = t_centers{j};
                    
                    if if_figure
                        if sort_info{sortInd}{2}{1,ERRORBAR} == 1 || sort_info{sortInd}{2}{1,ERRORBAR} == 3   % 95% CI
                            % Mean with shaded error bar
                            CI = sem*1.96;
                            h=shadedErrorBar(t_centers{j},ys ,CI,{'Color',lineColor,'LineStyle',lineStyle},transparent);
                            set(h.mainLine,'LineWidth',2);  hold on;
                        else
                            % Mean
                            h.mainLine=plot(t_centers{j},ys,'Color',lineColor,'LineStyle',lineStyle);
                            set(h.mainLine,'LineWidth',2);  hold on;
                        end
                    end
                    if_indicator = sort_info{sortInd}{2}{1,ERRORBAR} == 2 || sort_info{sortInd}{2}{1,ERRORBAR} == 3 && length(catsOut) == 2;
                    
                    % Add significance indicator if needed (only if we have and only have 2 inner conditions)
                    if if_indicator
                        if catNum_Out == 1 % Cache the first data
                            firstData = spike_hist{j}(selected_condition,:);
                        else % Draw indicators
                            secondData = spike_hist{j}(selected_condition,:);
                            
                            try
                                [~,pp]=ttest2(firstData,secondData); % two-sample t-test
                            catch
                                pp = NaN(1,size(firstData,2));
                            end
                            
                            % Cache the p values and sign for plot indicators later
                            if if_figure
                                pp_for_indicator(1,:) = pp;
                                data_diff_for_indicator(1,:) = mean(firstData,1) - mean(secondData,1);
                            end
                            
                            % Save data
                            PSTH{j,sortInd,k}.ps(1,:) = pp;
                        end
                        
                    end
                    
                    if if_figure
                        h_legend(catNum_Out) = h.mainLine;
                        txt_legend{catNum_Out} = [sort_info{sortInd}{2}{1,NOTES}{catNum_Out} ',' num2str(sum(selected_condition))];
                        axis tight
                    end
                    
                end % sort_conditions
                
                % Plotting indicators here. HH20150412
                if if_figure && if_indicator
                    
                    ylims = ylim;
                    
                    catNum_Out = 1;
                    pp = pp_for_indicator(catNum_Out,:);
                    data_diff = data_diff_for_indicator(catNum_Out,:);
                    
                    indicator_pos = zeros(size(data_diff)) + (catNum_Out-1)*ylims(2)/70;
                    
                    if PREF == LEFT  % LEFT is always at the bottom (in line with on-line figure). HH20150412
                        indicator_pos (data_diff < 0) = indicator_pos (data_diff < 0) + ylims(2);
                    else
                        indicator_pos (data_diff >= 0) = indicator_pos (data_diff >= 0) + ylims(2);
                    end
                    
                    lineColor = sort_info{sortInd}{2}{1,COLORS}(catNum_Out,:);
                    plot(t_centers{j}(pp<p_critical),indicator_pos(pp<p_critical),'s','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',5);
                    
                    
                    pp_for_indicator = [];
                    data_diff_for_indicator = [];
                    
                end
                
            elseif nest_levels == 2 % 2 nested conditions
                
                for catNum_Out = 1:length(catsOut)
                    for catNum_In = 1:length(catsIn)
                        
                        selected_condition = eval(sprintf('%s %s %s & %s %s %s',...
                            sort_info{sortInd}{2}{1,VAR},logOut{catNum_Out},num2str(catsOut(catNum_Out)),...
                            sort_info{sortInd}{2}{2,VAR},logIn{catNum_In},num2str(catsIn(catNum_In))));
                        
                        if sort_info{sortInd}{1}{1}(k) > 0  % Not all conditions
                            selected_condition = selected_condition & (stim_type_per_trial' == unique_stim_type(k));
                        end
                        
                        % If the sorting type is not outcome per se, we only
                        % choose the correct trials by default, except 0 heading.
                        if ~strcmp(sort_info{sortInd}{2}{1},'outcome_per_trial')
                            selected_condition = selected_condition & (~outcome_mask_enable | outcome_per_trial == outcome_mask | heading_per_trial' == 0); %% HH20140825
                        end
                        lineColor = sort_info{sortInd}{2}{1,COLORS}(catNum_Out,:);
                        lineStyle = sort_info{sortInd}{2}{2,LINESTYLES}{catNum_In,:};
                        
                        % ------- Data smoothing has been moved to the previous section --------- HH20141207
                        % ys = smooth(mean(spike_hist{j}(selected_condition,:),1),smoothFactor);
                        
                        % Anne Churchland NN, bin = 10 ms with Gaussian filter (sigma = 50 ms)
                        % ys = GaussSmooth(t_centers{j},mean(spike_hist{j}(selected_condition,:),1),smoothFactor);
                        
                        ys = mean(spike_hist{j}(selected_condition,:),1);
                        sem = std(spike_hist{j}(selected_condition,:),0,1)/sqrt(sum(selected_condition)); % @HH20160906
                        
                        % Save data
                        PSTH{j,sortInd,k}.raw{(catNum_Out-1)*length(catsIn)+catNum_In,1} = spike_hist{j}(selected_condition,:); % Save all raw data for further processing. HH20141118
                        
                        PSTH{j,sortInd,k}.ys((catNum_Out-1)*length(catsIn)+catNum_In,:) = ys;
                        PSTH{j,sortInd,k}.sem((catNum_Out-1)*length(catsIn)+catNum_In,:) = sem;
                        PSTH{j,sortInd,k}.ts = t_centers{j};
                        
                        if if_figure
                            
                            if sort_info{sortInd}{2}{2,ERRORBAR} == 1 || sort_info{sortInd}{2}{2,ERRORBAR} == 3 % 95% CI
                                % Mean with shaded error bar
                                CI = sem*1.96;
                                h=shadedErrorBar(t_centers{j},ys,CI,{'Color',lineColor,'LineStyle',lineStyle},transparent);
                                set(h.mainLine,'LineWidth',2);  hold on;
                            else
                                % Mean
                                h.mainLine=plot(t_centers{j},ys,'Color',lineColor,'LineStyle',lineStyle);
                                set(h.mainLine,'LineWidth',3);  hold on;
                            end
                        end
                        % Add significance indicator if needed (only if we have and only have 2 inner conditions)
                        
                        if_indicator = sort_info{sortInd}{2}{2,ERRORBAR} == 2 || sort_info{sortInd}{2}{2,ERRORBAR} == 3 && length(catsIn) == 2;
                        
                        if if_indicator
                            if catNum_In == 1 % Cache the first data
                                firstData = spike_hist{j}(selected_condition,:);
                            else % Draw indicators
                                secondData = spike_hist{j}(selected_condition,:);
                                try
                                    [~,pp]=ttest2(firstData,secondData); % two-sample t-test
                                catch
                                    pp = NaN(1,size(firstData,2));
                                end
                                
                                % Cache the p values and sign for plot indicators later
                                if if_figure
                                  pp_for_indicator(catNum_Out,:) = pp;
                                  data_diff_for_indicator(catNum_Out,:) = mean(firstData,1) - mean(secondData,1);
                                end
                                
                                % Save data
                                PSTH{j,sortInd,k}.ps(catNum_Out,:) = pp;
                            end
                            
                        end
                        if if_figure
                            try
                                h_legend((catNum_Out-1)*length(catsIn)+catNum_In) = h.mainLine;
                            catch  err
                                keyboard;
                            end
                            
                            txt_legend{(catNum_Out-1)*length(catsIn)+catNum_In} = [sort_info{sortInd}{2}{1,NOTES}{catNum_Out} ',' sort_info{sortInd}{2}{2,NOTES}{catNum_In} ', ' num2str(sum(selected_condition))];
                            
                            axis tight
                        end
                    end % In sort_conditions
                    
                end % Out sort_conditions
                
                
                % Plotting indicators here. HH20150412
                if if_figure && if_indicator
                    
                    ylims = ylim;
                    
                    for catNum_Out = 1: length(catsOut)
                        pp = pp_for_indicator(catNum_Out,:);
                        data_diff = data_diff_for_indicator(catNum_Out,:);
                        
                        indicator_pos = zeros(size(data_diff)) + (catNum_Out-1)*ylims(2)/70;
                        
                        if PREF == LEFT  % LEFT is always at the bottom (in line with on-line figure). HH20150412
                            indicator_pos (data_diff < 0) = indicator_pos (data_diff < 0) + ylims(2);
                        else
                            indicator_pos (data_diff >= 0) = indicator_pos (data_diff >= 0) + ylims(2);
                        end
                        
                        lineColor = sort_info{sortInd}{2}{1,COLORS}(catNum_Out,:);
                        plot(t_centers{j}(pp<p_critical),indicator_pos(pp<p_critical),'s','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',5);
                    end
                    
                    pp_for_indicator = [];
                    data_diff_for_indicator = [];
                end
            end
            
            debug = 0; % Under construction
            if debug && j == 2 && sortInd == 1  % Replot PSTH using SeriesComparison for better illustration. @HH20150527
                
                time_markers{1} = [0 mean(align_offsets_others{1})];
                time_markers{2} = [mean(align_offsets_others{2}) 0];
                
                % Wrong because SeriesComparison cannot yet receive data with different Ns
                h = SeriesComparison({PSTH{1,sortInd,k}.ys' PSTH{2,sortInd,k}.ys'},...
                    {PSTH{1,sortInd,k}.ts PSTH{2,sortInd,k}.ts time_markers},...
                    'Colors',{'k','k','m','m'},'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing');
            end
            
            if if_figure
                % Post-plot stuffs
                if k < length(stim_type_to_plot) ;set(gca,'xticklabel',''); end
                if k == 1 ;
                    if j==2
                        title([FILE 'unit' num2str(SpikeChan) ', reps = '  num2str(repetitionN) '     ' align_markers{j,4}],'color','k');
                    else
                        title(align_markers{j,4},'color','k');
                    end
                end
                
                if j > 1
                    set(gca,'yticklabel','');
                else
                    ylabel(sort_info{sortInd}{1}{2}{k});
                end
                
                ylims(l,:) = ylim;
                
                % Other time markers
                for jj = 1:length(align_markers{j,5})
                    mean_minus_std = mean(align_offsets_others{j}(:,jj)) - std(align_offsets_others{j}(:,jj));
                    mean_plus_std = mean(align_offsets_others{j}(:,jj)) + std(align_offsets_others{j}(:,jj));
                    if (mean_plus_std-mean_minus_std)>=10
                        plot([mean_minus_std mean_minus_std],[0 max(ylims(:,2))*1.05],'k--','linewidth',1.5);
                        plot([mean_plus_std mean_plus_std],[0 max(ylims(:,2))*1.05],'k--','linewidth',1.5);
                    else
                        temp_duration = round((mean_plus_std+mean_minus_std)/2/50)*50;
                        plot([temp_duration temp_duration],[0 max(ylims(:,2))*1.05],'k-','linewidth',2);
                    end
                    text((mean_minus_std+mean_plus_std)/2,0,align_markers{j,6}{jj},'VerticalAlignment','middle','Rotation',90,'color','k');
                end
                
                if j == 2
                    h_legends(k) = legend(h_legend,txt_legend,'Location','NorthWest','color','none');
                end
            end
            
        end % align_markers
    end  % stim_type
    
    if if_figure
        % Post-plot stuffs
        for l = 1: k*j
            %         axes(h_subplot(l));
            set(gcf,'CurrentAxes',h_subplot(l));
            
            ylim([0 max(ylims(:,2))*1.05]);
            plot([0 0],[0 max(ylims(:,2))*1.05],'k','linewidth',2);
            %         grid on;
            
        end
        SetFigure(13);
        set(h_legends,'box','on','FontSize',9);
    end
    
end % sort_bases

drawnow;

%% 4. Choice Divergence and Preference

fprintf('Choice Divergence...\n')

% Depends on "sort_info"
ALL_CHOICE = 1; CHOICE_ANGLE = 2; CHOICE_DIFFICULT = 3;

%{
% XXX @HH20150417.
% Re-determine the preferred direction according to Anne's method (see above)
% Note here the preferred direction is computed from j=2 (saccade-aligned)
mean_rates = PSTH{2,ALL_CHOICE,1}.ys(2*k-1:2*k,:);      % Note for sortInd = 1 (ALL_CHOICE), all conditions are put together
% first_larger = mean(mean_rates(1, rate_ts > -500 & rate_ts <= -200)) > mean(mean_rates(2,rate_ts > -500 & rate_ts <= -200));
first_larger = mean(mean_rates(1, rate_ts > -500 & rate_ts <= -50)) > mean(mean_rates(2,rate_ts > -500 & rate_ts <= -50));
%}

% To prevent from getting a weird result that three conditions have different preferred directions,
% now I decide to use the "global preferred direction" for CD and CP as well. @HH20150417
first_larger = 1; % Follow the preferred direction defined by PSTH before ("first" = PREF found at ~ LINE 370)

if if_figure
    set(figure(1899),'pos',[175-1600*~isempty(batch_flag)   90         902         394]); clf;
    if ~isempty(batch_flag)
        set(1899,'Visible','off');
    end
    
    set(figure(1999),'pos',[175-1600*~isempty(batch_flag)   448    902     394]); clf;
    if ~isempty(batch_flag)
        set(1999,'Visible','off');
    end
    
    % figure: pref - null
    set(figure(2009),'pos',[175-1600*~isempty(batch_flag)   600    902     394]); clf;
    if ~isempty(batch_flag)
        set(2009,'Visible','off');
    end
end

for j = 1:size(align_markers,1) % Include two align methods. @HH20150417
    
    %     j = 1;
    rate_ts{j} = PSTH{j,1,1}.ts;
    ChoiceDivergence_ALL{j} = NaN(3,length(rate_ts{j}));
    ChoiceDivergence_ALL_perm{j}.std = NaN(3,length(rate_ts{j})); % HH20180608
    ChoiceDivergence_ALL_perm{j}.p = NaN(3,length(rate_ts{j}));
    
    ChoiceDivergence_Difficult{j} = NaN(3,length(rate_ts{j}));
    ChoiceDivergence_Easy{j} = NaN(3,length(rate_ts{j}));
    Pref_Null_ALL{j} = NaN(3,length(rate_ts{j}));
    
    
    for stim_type = 1:length(unique_stim_type)  % Always output three conditions
        
        fprintf('>');
        
        k = find(stim_type == unique_stim_type);
        
        if ~isempty(k)   % We have this condition
            
            % Calculate auROC for each time bin
            for tt = 1:length(rate_ts{j})
                
                [auROC, ~, perm] = rocN(PSTH{j,ALL_CHOICE,1}.raw{2*k-1}(:,tt),...
                    PSTH{j,ALL_CHOICE,1}.raw{2*k}(:,tt), [], CD_perm_N);      % Default: first larger

                ChoiceDivergence_ALL{j}(stim_type,tt) = 2 * (auROC - 0.5);
                ChoiceDivergence_ALL_perm{j}.std(stim_type,tt) = 2 * perm.std;
                ChoiceDivergence_ALL_perm{j}.p(stim_type,tt) = perm.pValue;
                
                ChoiceDivergence_Difficult{j}(stim_type,tt) = 2 * (rocN(PSTH{j,CHOICE_DIFFICULT,k}.raw{1}(:,tt),...
                    PSTH{j,CHOICE_DIFFICULT,k}.raw{2}(:,tt)) - 0.5); % Fixed the factor of 2. HH20180608
                
                ChoiceDivergence_Easy{j}(stim_type,tt) = 2 * (rocN(PSTH{j,CHOICE_DIFFICULT,k}.raw{3}(:,tt),...
                    PSTH{j,CHOICE_DIFFICULT,k}.raw{4}(:,tt)) - 0.5);
                
                Pref_Null_ALL{j}(stim_type,tt) = PSTH{j,ALL_CHOICE,1}.ys(2*k-1,tt) - PSTH{j,ALL_CHOICE,1}.ys(2*k,tt);
                
            end
            
            
            if ~first_larger % If not, flip (Never flip now, see above @HH20150417)
                ChoiceDivergence_ALL{j}(stim_type,:) = - ChoiceDivergence_ALL{j}(k,:);
                ChoiceDivergence_Difficult{j}(stim_type,:) = - ChoiceDivergence_Difficult{j}(k,:);
                ChoiceDivergence_Easy{j}(stim_type,:) = - ChoiceDivergence_Easy{j}(k,:);
            end
            
        end  % if ~isempty(k)
        
    end
    
   
    if if_figure
        
        set(0,'defaultAxesColorOrder',[0 0 1; 1 0 0; 0 1 0; 1 0.8 0; 0 1 1;]);
        
        set(0,'currentfig',1899);
        subplot(1,2,j);
        plot(rate_ts{j}, ChoiceDivergence_ALL{j}','Linewidth',2); SetFigure();  hold on;
        plot(rate_ts{j}, ChoiceDivergence_ALL_perm{j}.std','linestyle','--','linew',1)
        axis tight;
        
        if j == 2 ;title([FILE 'unit' num2str(SpikeChan) ', reps = '  num2str(repetitionN) ', Choice Divergence']); end
        %     print(1899,'-dbitmap',[mat_file_fullname{i} '_ChoiceDivergenceAll.bmp']);
        
        set(0,'currentfig',1999);
        subplot(1,2,j);
        plot(rate_ts{j}, ChoiceDivergence_Easy{j}'- ChoiceDivergence_Difficult{j}','Linewidth',2); SetFigure(); axis tight;
        if j == 2 ;title([FILE 'unit' num2str(SpikeChan) ', reps = '  num2str(repetitionN) ', Easy - Difficult']); end
        %     print(1999,'-dbitmap',[mat_file_fullname{i} '_ChoiceDivergenceEasyMinusDifficult.bmp']);
        
        set(0,'currentfig',2009);
        subplot(1,2,j);
        plot(rate_ts{j}, Pref_Null_ALL{j}','Linewidth',2); SetFigure(); axis tight;
        
        if j == 2 ;title([FILE 'unit' num2str(SpikeChan) ', reps = '  num2str(repetitionN) ', Pref-Null']); end
        
        % Other time markers
        hold on; ylims = ylim;
        for jj = 1:length(align_markers{j,5})
            mean_minus_std = mean(align_offsets_others{j}(:,jj)) - std(align_offsets_others{j}(:,jj));
            mean_plus_std = mean(align_offsets_others{j}(:,jj)) + std(align_offsets_others{j}(:,jj));
            if (mean_plus_std-mean_minus_std)>=10
                temp_duration = round((mean_plus_std+mean_minus_std)/2/50)*50;
                plot([temp_duration temp_duration],[0 ylims(2)*1.05],'k--','linewidth',2);
            else
                temp_duration = round((mean_plus_std+mean_minus_std)/2/50)*50;
                plot([temp_duration temp_duration],[0 ylims(2)*1.05],'k-','linewidth',2);
            end
        end
        plot([0 0],[0 ylims(2)*1.05],'k','linewidth',2);
        
    end
    
end;

fprintf('\n');

%% Choice preference (related to "PREF" of this cell).  @HH20150418
%  Will be transformed to be related to "Contralateral" in GROUP_GUI
%  @HH20160915

% j = 2;  % This is not too much time-sensitive, so I choose j = 2 (aligned to sac onset)
% choice_or_mod_pref_timewin = {
%     mean(align_offsets_others{2}(:,1)) <= rate_ts{2} & rate_ts{2} <= 0;   % Stimlus onset - saccade onset
%     0 < rate_ts{2} & rate_ts{2} <= inf;   % Postsaccade period
%     };
choice_or_mod_pref_timewin = {
    2, mean(align_offsets_others{2}(:,1)) <= rate_ts{2} & rate_ts{2} <= 0;   % Stimlus onset - saccade onset
    2, 0 < rate_ts{2} & rate_ts{2} <= inf;   % Postsaccade period
    1, 0 <= rate_ts{1} & rate_ts{1} <= mean(align_offsets_others{1}(:,1)); % I added a new choice preference which only includes stim-on to stim-off to better select out the ramping cells. HH20160918
    };

ChoicePreference = nan(length(choice_or_mod_pref_timewin),3);
ChoicePreference_pvalue = nan(length(choice_or_mod_pref_timewin),3);

for stim_type = 1:3  % Always output three conditions
    k = find(stim_type == unique_stim_type,1);
    
    if ~isempty(k)   % We have this condition
        
        %         2 * (rocN(mean(PSTH{j,ALL_CHOICE,1}.raw{2*k-1}(:,choice_pref_timewin),2),...
        %             mean(PSTH{j,ALL_CHOICE,1}.raw{2*k}(:,choice_pref_timewin),2)) - 0.5)
        
        for cmpt = 1:length(choice_or_mod_pref_timewin)
            
            % Here I let CP_HH calculate the auROC AND p_value (permutation method) for me
            %             fake_spk_pref = mean(PSTH{choice_or_mod_pref_timewin{cmpt,1},ALL_CHOICE,1}.raw{2*k-1}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            %             fake_spk_null = mean(PSTH{choice_or_mod_pref_timewin{cmpt,1},ALL_CHOICE,1}.raw{2*k}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            %             fake_heading = zeros(length(fake_spk_pref)+length(fake_spk_null),1);
            %             fake_choice = [ones(length(fake_spk_pref),1) * RIGHT; ones(length(fake_spk_null),1) * LEFT];
            %             temp = CP_HH(fake_heading, fake_choice, [fake_spk_pref ; fake_spk_null],1000,0);
            
            % Not anymore. Use the new rocN(with permutation)'s functionality. HH20180607
            
            fake_spk_pref = mean(PSTH{choice_or_mod_pref_timewin{cmpt,1},ALL_CHOICE,1}.raw{2*k-1}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            fake_spk_null = mean(PSTH{choice_or_mod_pref_timewin{cmpt,1},ALL_CHOICE,1}.raw{2*k}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            
            [auROC, ~, perm] = rocN (fake_spk_pref, fake_spk_null , [], ChoicePref_perm_N);
            
            ChoicePreference(cmpt,stim_type) = 2*(auROC-0.5);
            ChoicePreference_pvalue(cmpt,stim_type) =  perm.pValue ;
            
        end
    end
end

if if_figure
    set(0,'currentfig',1899);  xlims = xlim; ylims = ylim;
    text(xlims(1)*0.9,ylims(1)*0.7,sprintf('ChoicePref = %s',num2str(ChoicePreference(1,:))));
end

%% 5. Modality Divergence and Modality Preference  @HH20150418

modality_pair = {[1 2],[1 3],[2 3]};  % The later is set to be "Pref modality"

for j = 1:size(align_markers,1) % Include two align methods.
    
    ModalityDivergence{j} = NaN(3,length(rate_ts{j}));
    ModalityDivergence_choice_separate = NaN(2,length(rate_ts{j}));
    
    if length(unique_stim_type) == 3 % If we have all three modalities
        
        raw_for_md = PSTH{j,ALL_CHOICE,1}.raw;
        
        % Now calculate auROC
        for mp = 1:length(modality_pair)
            
            pref_mod = modality_pair{mp}(2);
            null_mod = modality_pair{mp}(1);
            
            % Calculate modality divergence for Pref and Null choices separately and then average them. (Anne 2014)
            for choice = 1:2
                % Calculate auROC for each time bin
                for tt = 1:length(rate_ts{j})
                    ModalityDivergence_choice_separate(choice,tt) = ...
                        2 * (rocN(raw_for_md{2*(pref_mod-1)+choice}(:,tt), ...
                        raw_for_md{2*(null_mod-1)+choice}(:,tt)) - 0.5);      % Default: first larger
                end
            end
            
            % Averaged modality divergence
            ModalityDivergence{j}(mp,:) = mean(ModalityDivergence_choice_separate);
        end
    end
    
    if if_figure
        set(0,'currentfig',1899);  subplot(1,2,j); hold on;
        plot(rate_ts{j},ModalityDivergence{j}(1,:),'color',[0.5 0.5 0.5],'linew',2);   % Visual - vest
        plot(rate_ts{j},ModalityDivergence{j}(3,:),'color',[0.5 0.5 0.5],'linew',2,'linestyle',':');   % Comb - visual
    end
end

%% --- Modality Preference ---

j = 2;  % This is not too much time-sensitive, so I choose j = 2 (aligned to sac onset)

ModalityPreference = nan(length(choice_or_mod_pref_timewin),3);
ModalityPreference_pvalue = nan(length(choice_or_mod_pref_timewin),3);


if length(unique_stim_type) >= 3 % If we have all three modalities
    
    for mp = 1:length(modality_pair)
        
        pref_mod = modality_pair{mp}(2);
        null_mod = modality_pair{mp}(1);
        
        for cmpt = 1:length(choice_or_mod_pref_timewin)
            
            raw_for_md = PSTH{choice_or_mod_pref_timewin{cmpt,1},ALL_CHOICE,1}.raw;
            
            % I let CP_HH.m calculate the auROC AND p_value (permutation method) for me.
            % Note here I use the grand CP because I combine PREF and NULL trials together
            % by assigning PREF and NULL choices to 0 and 1 fake heading, respectively.
            fake_spk_pref_mod_pref_choice = mean(raw_for_md{2*(pref_mod-1)+1}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            fake_spk_null_mod_pref_choice = mean(raw_for_md{2*(null_mod-1)+1}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            fake_spk_pref_mod_null_choice = mean(raw_for_md{2*(pref_mod-1)+2}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            fake_spk_null_mod_null_choice = mean(raw_for_md{2*(null_mod-1)+2}(:,choice_or_mod_pref_timewin{cmpt,2}),2);
            
            fake_0_heading = zeros(length(fake_spk_pref_mod_pref_choice)+length(fake_spk_null_mod_pref_choice),1);
            fake_1_heading = ones(length(fake_spk_pref_mod_null_choice)+length(fake_spk_null_mod_null_choice),1);
            
            fake_0_heading_choice = [ones(length(fake_spk_pref_mod_pref_choice),1) * RIGHT; ones(length(fake_spk_null_mod_pref_choice),1) * LEFT];
            fake_1_heading_choice = [ones(length(fake_spk_pref_mod_null_choice),1) * RIGHT; ones(length(fake_spk_null_mod_null_choice),1) * LEFT];
            
            temp = CP_HH([fake_0_heading; fake_1_heading], [fake_0_heading_choice; fake_1_heading_choice],...
                [fake_spk_pref_mod_pref_choice; fake_spk_null_mod_pref_choice; fake_spk_pref_mod_null_choice; fake_spk_null_mod_null_choice],...
                ModalityPref_perm_N,0);
            
            ModalityPreference(cmpt,mp) = 2*(temp.CP_grand-0.5);
            ModalityPreference_pvalue(cmpt,mp) =  temp.CP_grand_p_perm ;
        end
    end
end

if if_figure
    set(0,'currentfig',1899); xlims = xlim; ylims = ylim;
    text(xlims(1)*0.9,ylims(1)*0.9,sprintf('ModPref = %s',num2str(ModalityPreference(1,:))));
    
    drawnow;
end

%% 6. Calculate sliding CP and Neuro-threshold (Rewrite, with CP_HH function. @HH20140925)
%  %{

fprintf('CP...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time windows
binSize_CP = 250;  % in ms  % According to Shadlen 2001
stepSize_CP = 10; % in ms  % Match PSTH, CD, and others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving
% result.binSize_CP = binSize_CP;
% result.stepSize_CP = stepSize_CP;


%------- Preferred direction for CP ---------
% Here I decide to use a fixed pref. direction for each cell when
% calculating CP (*** same across time ***)
% E.g.: Anne 2014 paper: 100-200ms before decision. @ HH20141126

% To prevent from getting a weird result that three conditions have different preferred directions,
% now I decide to use a "global preferred direction" (*** same across time AND modalities***)
% for CD and CP as well. @HH20150417

PREF_CP = PREF;

% Nevertheless, I still keep this part to calculate a PREF_CP_obsolete for
% each modalities

for k = 1: length(unique_stim_type)   % For each stim type
    selected_condition = stim_type_per_trial' == unique_stim_type(k);
    for hh = 1:length(unique_heading)
        curr_heading = (selected_condition) & (heading_per_trial' == unique_heading(hh));
        resp_mean(hh) = sum(mean(spike_aligned{1,2}(curr_heading, -1000 <= spike_aligned{2,2} &  spike_aligned{2,2} <= -100)));
    end
    
    [rr,pp] = corrcoef(unique_heading, resp_mean);
    PREF_CP_obsolete(k) = (rr(1,2) <= 0) * LEFT + (rr(1,2) > 0) * RIGHT;
end

%------- END Preferred direction for CP ---------

CP = cell(size(align_markers,1),length(unique_stim_type));

for j = 1:2
    CP_ts{j} = align_markers{j,2} + binSize_CP/2 : stepSize_CP : align_markers{j,3} - binSize_CP/2; % Centers of CP time windows
end

% Fisher info
FIAngles = {1:length(unique_heading),3:length(unique_heading)-2};  % {all headings, +/- 2 degrees}
fisherSimpleGu = nan(length(CP_ts{1}),3,length(FIAngles),5);  % Time, Stim type, FIangles, FIs


for k = 1: length(unique_stim_type)   % For each stim type
    selected_condition = stim_type_per_trial' == unique_stim_type(k);
    headings = heading_per_trial(selected_condition);
    choices = choice_per_trial(selected_condition);
    
    %------- Calculate Psycho function once for each condition
    for hh = 1:length(unique_heading)
        num_headings(hh,1) = sum(headings == unique_heading(hh));
        rightward_prop(hh,1) = sum(choices(headings == unique_heading(hh)) == RIGHT)/ sum(headings == unique_heading(hh));
    end
    
    % Fitting
    [Psy_bias,Psy_thres] = cum_gaussfit_max1([unique_heading, rightward_prop, num_headings],method,0);
    [Psy_bias_tol,Psy_thres_tol] = cum_gaussfit_max1([unique_heading, rightward_prop, num_headings],method,tolerance);
    
    %------- Time sensitive stuffs
    
    for j = 1:size(align_markers,1)  % For each temporal alignment

        fprintf('.')
        
        % Saving Psy to both of the time alignements
        CP{j,k}.Psy_func = [unique_heading, rightward_prop, num_headings];
        CP{j,k}.Psy_para = [Psy_bias,Psy_thres];
        CP{j,k}.Psy_para_tol = [Psy_bias_tol,Psy_thres_tol];
        
        % CP_ts{j} = align_markers{j,2} + binSize_CP/2 : stepSize_CP : align_markers{j,3} - binSize_CP/2; % Centers of CP time windows
        
        % Preallocation
        CP{j,k}.CP_grand = zeros(1,length(CP_ts{j}));
        CP{j,k}.CP_p = zeros(1,length(CP_ts{j}));
        CP{j,k}.CP_std = zeros(1,length(CP_ts{j}));
        CP{j,k}.Neu_thres = zeros(1,length(CP_ts{j}));
        CP{j,k}.ts = CP_ts{j};
        

        for tt = 1:length(CP_ts{j})    % For each CP time window
            
            winBeg = ceil(((tt-1) * stepSize_CP) / spike_timeWin) + 1;
            winEnd = ceil(((tt-1) * stepSize_CP + binSize_CP) / spike_timeWin) + 1;
            
            % Time sliced spike counts
            spike_counts = sum(spike_aligned{1,j}(selected_condition,winBeg:winEnd),2) / binSize_CP * 1000;
            
            if ~isempty(batch_flag)
                CP_result = CP_HH(headings,choices,spike_counts,CP_perm_N,0);    % Permutation for Batch
                
                if CP_perm_N > 0
                    CP{j,k}.CP_p(tt) = CP_result.CP_grand_p_perm;
                else
                    CP{j,k}.CP_p(tt) = CP_result.CP_grand_p_ttest;
                end
                
                CP{j,k}.CP_std(tt)  = CP_result.CP_grand_std_perm;
            else
                CP_result = CP_HH(headings,choices,spike_counts,-1,0);    % Fast version
                CP{j,k}.CP_p(tt) = CP_result.CP_grand_p_ttest;
                CP{j,k}.CP_std(tt)  = NaN;
            end
            
            % Keep all data first
            CP{j,k}.raw_CP_result{tt} = CP_result;
            
            % For quick access
            CP{j,k}.CP_grand(tt) = CP_result.CP_grand;
            CP{j,k}.Neu_thres(tt) = CP_result.Neu_para_anti(2);
            CP{j,k}.pref_CP(tt) = CP_result.pref;
            
            % ============ Calculate Fisher information (Simple method like Gu2010) ========== HH20180619
            if j == 1 % Only care about j = 1
                % Get data
                thisTuning = CP_result.Neu_tuning;
                
                if PREF == LEFT % To keep PREF = RIGHT (Actually no effect because we have slope^2 afterwards in Fisher information)
                    thisTuning(:,2:3) = flipud(thisTuning(:,2:3));
                end
                
                for AA = 1:length(FIAngles)  % Large and small delta_theta
                    linearFit = polyfit(thisTuning(FIAngles{AA},1),thisTuning(FIAngles{AA},2),1);
                    slopeInRad = linearFit(1)*(180/pi);
                    varPoisson = mean(thisTuning(FIAngles{AA},2)); % Assuming Poisson with fano = 1
                    varReal = mean((thisTuning(FIAngles{AA},3) .* ...
                        sqrt(sum(CP_result.Neu_tuning_Dora_matrix_n(:,FIAngles{AA}))')).^2); % Using real variance (usually much larger)
                    
                    % Compute simple Fisher (simple method like Gu)
                    thisFisherVarPoisson = slopeInRad^2/varPoisson;
                    thisFisherVarReal = slopeInRad^2/varReal;
                    
                    % Save
                    fisherSimpleGu(tt,unique_stim_type(k),AA,:) = [thisFisherVarPoisson, thisFisherVarReal, slopeInRad, varPoisson, varReal];
                end
           end
            
        end
        
       
        if ~isempty(CP_ts{j}) % In case there are no CP windows
            
            % Flip CP if local pref direction (of each time bin) is not aligned with global pref. direction
            % That is, now the preferred direction is stable across *** time AND modality *** (@HH20150417)
            need_flip = (CP{j,k}.pref_CP ~= PREF_CP);
            CP{j,k}.CP_grand(need_flip) = 1 - CP{j,k}.CP_grand(need_flip);
            
            % Center CP (For saccade alignment: Similar to the CP we called before)
            center_t = mean(align_offsets_others{j}(:,2) + align_offsets_others{j}(:,1))/2; % Other markers 2 - 1
            [~,center_t_ind] = min(abs(CP_ts{j}-center_t));
            CP{j,k}.CP_grand_center = CP{j,k}.CP_grand(center_t_ind);
            
            % Pre-alignment CP (For saccade alignment: The first CP window that has the saccade onset timepoint)
            sac_t = - binSize_CP/2;
            [~,sac_t_ind] = min(abs(CP_ts{j}-sac_t));
            CP{j,k}.CP_grand_sac = CP{j,k}.CP_grand(sac_t_ind);
            
            %             % Other non-time sensitive
            %             CP{j,k}.Psy_para = CP_result.Psy_para;
            %             CP{j,k}.Psy_func = CP_result.Psy_func;
            
        else
            CP{j,k}.CP_grand_center = NaN;
            CP{j,k}.CP_grand_sac = NaN;
        end
        
        
    end
end

fprintf('\n')

%% Neuro tuning of different time duration (center, pre, post) % HH20150403
% Note: This should be modified so that only the correct trials are
% included !! (Already done. see below)
if if_figure
    
    set(figure(2000),'position',[383-1600*~isempty(batch_flag) 180 935 548]); clf;
    if ~isempty(batch_flag)
        set(2000,'Visible','off');
    end
    
    j = 2;
    
    for k = 1:length(unique_stim_type)   % For each stim type
        real_k = unique_stim_type(k);
        
        subplot(2,3,1); hold on;
        center_tuning = CP{j,k}.raw_CP_result{center_t_ind}.Neu_tuning ;
        errorbar(center_tuning(:,1),center_tuning(:,2),center_tuning(:,3),'color',stim_type_colors(real_k,:),'LineWid',2);
        axis tight; title('Center'); ylabel('All');
        
        subplot(2,3,2); hold on;
        sac_tuning = CP{j,k}.raw_CP_result{sac_t_ind}.Neu_tuning;
        errorbar(sac_tuning(:,1),sac_tuning(:,2),sac_tuning(:,3),'color',stim_type_colors(real_k,:),'LineWid',2);
        axis tight; title('Pre-Sac');
        
        subplot(2,3,3); hold on;
        post_tuning = CP{j,k}.raw_CP_result{end}.Neu_tuning;
        errorbar(post_tuning(:,1),post_tuning(:,2),post_tuning(:,3),'color',stim_type_colors(real_k,:),'LineWid',2);
        axis tight; title('Post-Sac');
        
        subplot(2,3,4); hold on; ylabel('Correct only');
        center_tuning = CP{j,k}.raw_CP_result{center_t_ind}.Neu_tuning_correctonly ;
        errorbar(center_tuning(:,1),center_tuning(:,2),center_tuning(:,3),'color',stim_type_colors(real_k,:),'LineWid',2);
        axis tight; title('Center');
        
        subplot(2,3,5); hold on;
        sac_tuning = CP{j,k}.raw_CP_result{sac_t_ind}.Neu_tuning_correctonly;
        errorbar(sac_tuning(:,1),sac_tuning(:,2),sac_tuning(:,3),'color',stim_type_colors(real_k,:),'LineWid',2);
        axis tight; title('Pre-Sac');
        
        subplot(2,3,6); hold on;
        post_tuning = CP{j,k}.raw_CP_result{end}.Neu_tuning_correctonly;
        errorbar(post_tuning(:,1),post_tuning(:,2),post_tuning(:,3),'color',stim_type_colors(real_k,:),'LineWid',2);
        axis tight; title('Post-Sac');
        
        
        SetFigure(15);
    end
    
    %% Plotting Psychometric, Neurometric functions, and CP
    set(figure(161),'position',[45-1300*~isempty(batch_flag)  126   1111 704]); clf;  clf
    if ~isempty(batch_flag)
        set(161,'Visible','off');
    end
    subplot_tight(2,3,1,0.07,[0.1 0.15 0.05 0.05]);
    
    for k = 1:length(unique_stim_type)
        Psy_func = CP{2,k}.Psy_func(:,2);
        Psy_para = CP{2,k}.Psy_para;
        plot(unique_heading,Psy_func,'o','color',stim_type_colors(unique_stim_type(k),:),'markerfacecolor',stim_type_colors(unique_stim_type(k),:)); hold on
        xx = min(unique_heading):0.01:max(unique_heading);
        plot(xx,cum_gaussfit(Psy_para,xx),'color',stim_type_colors(unique_stim_type(k),:));
        
        text(max(unique_heading)*0.5,0.25-k*0.07,sprintf('%6.3g',Psy_para(2)),'color',stim_type_colors(unique_stim_type(k),:));
    end
    xlim([min(unique_heading),max(unique_heading)]);
    
    
    for j = 1:2
        
        subplot_tight(2,3,[2+3 * (j-1): 3+3 * (j-1)],0.07,[0.1 0.1 0.1 0.05]);
        
        % j = 2;   % Align to sac
        
        for k = 1:length(unique_stim_type)
            % Plot neurometric and CP
            xx = CP{j,k}.ts;
            
            if k==1
                [ax h_(k,1),h_(k,2)] = plotyy(xx, CP{j,k}.CP_grand, xx, CP{j,k}.Neu_thres);
            else
                %         axes(ax(1));
                set(gcf,'CurrentAxes',ax(1));
                hold on; h_(k,1) = plot(xx,CP{j,k}.CP_grand);
                %         axes(ax(2));
                set(gcf,'CurrentAxes',ax(2));
                hold on; h_(k,2) = plot(xx,CP{j,k}.Neu_thres);
            end
            
            set(h_(k,1),'linewidth',1,'marker','o','markersize',10,'color',stim_type_colors(unique_stim_type(k),:))
            set(h_(k,2),'linewidth',1,'linestyle','--','marker','^','markersize',10,'color',stim_type_colors(unique_stim_type(k),:))
            
            set(gcf,'CurrentAxes',ax(1));
            hold on;
            
            plot(xx(CP{j,k}.CP_p < 0.05),CP{j,k}.CP_grand(CP{j,k}.CP_p < 0.05),...
                'o','color',stim_type_colors(unique_stim_type(k),:),'markerfacecolor',stim_type_colors(unique_stim_type(k),:),'markersize',10);   % Sig. CP
            
            % Plot the "preferred direction" for each window. HH20140602
            %     plot(xx,ylimCP(1)+(line_re_shift{k}>0)*(range(ylimCP)-0.06)+(3-k)*0.02,['s' colors{unique_stim_type(k)}],'markerfacecolor',colors{unique_stim_type(k)},'markersize',5);
            
            % Plot Psychometric threshold
            psy_thres = CP{1,k}.Psy_para(2);
            set(gcf,'CurrentAxes',ax(2)); hold on;
            plot([xx(1) xx(end)],[psy_thres psy_thres],'-','color',stim_type_colors(unique_stim_type(k),:));
            
        end
        
        set(ax,'xlim',[xx(1) xx(end)]);
        
        ylimCP = [0 1];
        set(gcf,'CurrentAxes',ax(1));  hold on;
        lims = axis; hold on;
        plot([lims(1) lims(2)],[0.5 0.5],'k--');
        ylim(ylimCP);
        set(ax(1),'xtick',[],'ytick',[ylimCP(1):0.1:ylimCP(2)],'ycolor','k');
        ylabel(ax(1),'Grand CP (O)')
        
        set(gcf,'CurrentAxes',ax(2));ylim([1 200]);
        
        if j == 1
            % Center and sac markers for CP
            plot(ax(1),[CP_ts{j}(center_t_ind) CP_ts{j}(center_t_ind)],ylimCP,'k--');
            plot(ax(1),[CP_ts{j}(sac_t_ind) CP_ts{j}(sac_t_ind)],ylimCP,'k--');
            
            title([FILE 'unit' num2str(SpikeChan) ', reps = '  num2str(repetitionN) '     ' align_markers{j,4}],'color','k');
        else
            % Center and sac markers for CP
            plot(ax(1),[CP_ts{j}(center_t_ind) CP_ts{j}(center_t_ind)],ylimCP,'k--');
            plot(ax(1),[CP_ts{j}(sac_t_ind) CP_ts{j}(sac_t_ind)],ylimCP,'k--');
            
            xlabel(ax(2),['Center of ' num2str(binSize_CP) ' ms time window (ms)']);
        end
        
        ylabel(ax(2),sprintf('Neuronal Threshold (\\Delta anti-model)'))
        set(ax(2),'yscale','log','yticklabel',[1 10 100],'ytick',[1 10 100],'ycolor','k');
        
    end
    
    SetFigure(15);
    
end
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
if exist('CP','var')  % Full version
    result = PackResult(FILE, SpikeChan, repetitionN, unique_stim_type, ... % Obligatory!!
        stim_type_per_trial, heading_per_trial, choice_per_trial,... % Trial info
        align_markers, align_offsets_others, sort_info, smoothFactor ,...
        PREF, PREF_CP_obsolete, PREF_target_location, outcome_mask_enable, ...
        binSize_rate, stepSize_rate, rate_ts, binSize_CP, stepSize_CP, CP_ts ,...
        spike_aligned, spike_hist, CP, PSTH, fisherSimpleGu,...
        ChoiceDivergence_ALL, ChoiceDivergence_ALL_perm, ChoiceDivergence_Difficult, ChoiceDivergence_Easy,ChoicePreference,ChoicePreference_pvalue,...
        ModalityDivergence, ModalityPreference, ModalityPreference_pvalue);
    % Figures to save
    if if_figure
        config.save_figures = [60 + (1:length(sort_info)) ,161, 1899, 1999, 2000, 59,2009];
    else
        config.save_figures = [];
    end
    
else
    result = PackResult(FILE, SpikeChan, repetitionN, unique_stim_type, ... % Obligatory!!
        stim_type_per_trial, heading_per_trial, choice_per_trial,... % Trial info
        align_markers, align_offsets_others, sort_info, smoothFactor ,...
        PREF, PREF_target_location, outcome_mask_enable, ...
        binSize_rate, stepSize_rate, rate_ts,...
        spike_aligned, spike_hist, PSTH, fisherSimpleGu, ...
        ChoiceDivergence_ALL, ChoiceDivergence_ALL_perm, ChoiceDivergence_Difficult, ChoiceDivergence_Easy,ChoicePreference,ChoicePreference_pvalue,...
        ModalityDivergence, ModalityPreference, ModalityPreference_pvalue);
    
    %     result = PackResult(FILE, SpikeChan, repetitionN, unique_stim_type, ... % Obligatory!!
    %                         PREF_target_location);
    
    % Figures to save
    config.save_figures = [];
    
end

config.suffix = 'PSTH';
config.xls_column_begin = 'HD_rep';
% config.xls_column_end = 'HD_comb_p';
config.xls_column_end = 'HD_comb_ChoicePref_p';

% Only once
config.sprint_once_marker = 'gs';
config.sprint_once_contents = 'result.repetitionN,num2str(result.PSTH{2,1,1}.ts)';
% Loop across each stim_type
% config.sprint_loop_marker = {'gg';
%                            'gg';
%                            'sss'};
% config.sprint_loop_contents = {'result.CP{2,k}.Psy_para(2), result.CP{2,k}.Psy_para(1)';
%                         'result.CP{2,k}.CP_grand_center, result.CP{2,k}.CP_grand_sac';
%                         'num2str(result.PSTH{2,1,1}.ys((k-1)*2+1,:)), num2str(result.PSTH{2,1,1}.ys((k-1)*2+2,:)), num2str(result.PSTH{2,1,1}.ps(k,:))'};

% Replace the meaningless CPs with ChoicePreference (surprisingly, they have the same abbreviation!). HH20160419
config.sprint_loop_marker = {'gg';
    'gg';
    'gg'};
config.sprint_loop_contents = {'result.CP{2,k}.Psy_para(2), result.CP{2,k}.Psy_para(1)';
    'result.CP{2,k}.CP_grand_center, result.CP{2,k}.CP_grand_sac';
    'result.ChoicePreference(1,k), result.ChoicePreference_pvalue(1,k)'};

config.append = 1; % Overwrite or append
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveResult(config, result);


%%%%%%%%%%%%%%%%%%%%%  Output   HH20140510 / HH20140621 / HH20141003 %%%%%%%%%%%%%%%%%

% if ~isempty(batch_flag)  % Figures and raw data (always in "result" structure)
%
%     %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     suffix = ['PSTH'];
%
%     % Figures to save
%     save_figures = [60 + (1:length(sort_info)) ,161];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     outpath = ['Z:\Data\Tempo\Batch\' batch_flag(1:end-2) '\'];
%
%     % Check directory
%     if ~exist(outpath,'dir')
%         mkdir(outpath);
%     end
%     savefilename = [outpath [FILE '_' num2str(SpikeChan)] '_' suffix];
%
%     % Delete existing data files
%     if exist([savefilename '.mat'],'file')
%         delete([savefilename '*.*']);
%     end
%
%     % Save raw data
%     save(savefilename,'result');
%
%     % Save figures
%     for ff = 1:length(save_figures)
% %         orient landscape;
%         set(save_figures(ff),'Visible','on');
%         print(save_figures(ff),'-dbitmap',[savefilename '_fig_' num2str(ff) '.bmp']);
%         close(save_figures(ff));
% %         saveas(save_figures(ff),[savefilename '_fig_' num2str(ff)],'bmp');
%     end
%
% end
%
%
% % Print part of data to texts (clipboard or .dat file)
%
% %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Only once
% sprint_once_marker_temp = 'gs';
% sprint_once_contents = 'repetitionN,num2str(PSTH{2,1,1}.ts)';
% % Loop across each stim_type
%     % sprint_loop_marker_temp = {'sss';
%     %                        };
%     % sprint_loop_contents = {'num2str(PSTH{2,1,1}.ys((k-1)*2+1,:)), num2str(PSTH{2,1,1}.ys((k-1)*2+2,:)), num2str(PSTH{2,1,1}.ps(k,:))';
%     %                        };
%
% % sprint_loop_marker_temp = {};
% % sprint_loop_contents = {};
% sprint_loop_marker_temp = {'gg';
%                            'gg';
%                            'sss'};
% sprint_loop_contents = {'CP{2,k}.Psy_para(2), CP{2,k}.Psy_para(1)';
%                         'CP{2,k}.CP_grand_center, CP{2,k}.CP_grand_sac';
%                         'num2str(PSTH{2,1,1}.ys((k-1)*2+1,:)), num2str(PSTH{2,1,1}.ys((k-1)*2+2,:)), num2str(PSTH{2,1,1}.ps(k,:))'};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % HD_rep	HD_vest_psy_thres	HD_vis_psy_thres	HD_comb_psy_thres	HD_vest_CP_grand_center	HD_vest_CP_grand_sac	HD_vis_CP_grand_center	HD_vis_CP_grand_sac	HD_comb_CP_grand_center	HD_comb_CP_grand_sac
%
%
% sprint_once_marker = [];
% for i = 1:length(sprint_once_marker_temp)
%     sprint_once_marker = [sprint_once_marker '%' sprint_once_marker_temp(i) '\t '];
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
%         fprintf(fid, ['FILE\t ' sprint_once_contents '|\t']);
%
%         for ll = 1:length(sprint_loop_contents)
%             fprintf(fid,[sprint_loop_contents{ll} '|\t']);
%         end
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
% for ll = 1:length(sprint_loop_contents)
%
%     sprint_loop_marker = [];
%     for i = 1:length(sprint_loop_marker_temp{ll})
%         sprint_loop_marker = [sprint_loop_marker '%' sprint_loop_marker_temp{ll}(i) '\t '];
%     end
%
%     for conditions = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
%         if sum(unique_stim_type == conditions)==0
%             buff = sprintf(sprint_loop_marker,ones(1,sum(sprint_loop_marker=='%'))*NaN);
%         else
%             k = find(unique_stim_type == conditions);
%             eval(['buff = sprintf(sprint_loop_marker,' sprint_loop_contents{ll} ');']);
%         end
%         fprintf(fid, '%s', buff);
%         toClip = [toClip sprintf('%s', buff)];
%     end
%
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
%     suffix = ['PSTH'];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     for sortInd = 1:length(sort_info)
%         savefilename = [outpath [FILE '_' num2str(SpikeChan)] '_' suffix '_' num2str(sortInd) '.png'];
%         if exist(savefilename)
%             delete(savefilename);
%         end
%         saveas(60+sortInd,savefilename,'png');
%     end
%
%     % Print PSTHs
%
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
%     fprintf(fid,'%s\t %g\t %s\t',[FILE '_' num2str(SpikeChan)],repetitionN,num2str(PSTH{2,1,1}.ts));
%
%     for conditions = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
%         if sum(unique_stim_type == conditions)==0
%             buff = sprintf(sprint_txt, ones(1,length(sprint_txt_temp))*NaN);
%         else
%             k = find(unique_stim_type == conditions);
%             %%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             buff = sprintf(sprint_txt, num2str(PSTH{2,1,1}.ys((k-1)*2+1,:)), num2str(PSTH{2,1,1}.ys((k-1)*2+2,:)));  % Fig = 1, row = 1, column(alignmarker) = 2
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