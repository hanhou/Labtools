%-----------------------------------------------------------------------------------------------------------------------
%-- Firing rate analysis for LIP neurons during Reaction Time Heading
%Discrimination task
%--	05/16/10 Adhira
%-----------------------------------------------------------------------------------------------------------------------

function Reaction_Time_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME,:);
temp_spike_data = data.spike_data(SpikeChan,:);
temp_spike_rates = data.spike_rates(SpikeChan,:);

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials = ((trials >= BegTrial) & (trials <= EndTrial));

stim_type = temp_stim_type(select_trials);
heading = temp_heading(select_trials);
amplitude= temp_amplitude(select_trials);
num_sigmas= temp_num_sigmas(select_trials);
total_trials = temp_total_trials(select_trials);
spike_rates = temp_spike_rates(select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');

% spike_data(find(spike_data>100)) = 1;
% event_data = data.event_data(1, :, :);
% event_data(find(event_data>100)) = 1;

%--------------------------------------------------------------------------
% 1. Split up spike data by trials (fixed at 5000 events/trial)
% 2. N.B.: spike_data collection is 1000ms before stim start and 4000ms
%after. So, need to replace everything before trial start and after trial
%end with NaNs so they can be ignored during averaging later (eNaN func).
% 3. Spike rate in bins for each trial - 50ms.
% 4. Pick only correct rials
% 5. Separate trials by in/out of RF, diff. headings and stim type. Then
% average across trials accordingly

%--------------------------------------------------------------------------

%1. Split up spike data into different trials (5000xtotal_trials). Sampling
%rate of 1/ms - total data length for each trial (irrespective of actual
%trial durations) is 5000 events.

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
data_duration = length(temp_spike_data)/length(temp_azimuth); %usually 5000
Discard_trials = find(trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*data_duration+1) :  Discard_trials(i)*data_duration ) = 9999;
end
spike_data_good = temp_spike_data(temp_spike_data~=9999);
%Each column is a new trial
spike_data = reshape(spike_data_good, data_duration, length(trials)); %length(total_trials)
size(spike_data)
% spike_data_test = spike_data(1:data_duration,1:length(total_trials));

% timebin for plot PSTH
timebin=35; %ms
% length of data for each trial - usually 5000
data_length=length(temp_spike_data)/length(temp_total_trials);
% x-axis for PSTH
num_time_bins=1:(data_length/timebin);

% Calculate trial durations and remove bad data


for i = 1:length(total_trials)
    trial_begin = find(data.event_data(1,:,i) == TRIAL_START_CD);
    trial_over = find(data.event_data(1,:,i) == TRIAL_END_CD);
    fix_pt_begin = find(data.event_data(1,:,i) == FP_ON_CD);
    stim_begin = find(data.event_data(1,:,i) == VSTIM_ON_CD ); %constant at 996...
    stim_over = find(data.event_data(1,:,i) == VSTIM_OFF_CD );
    monk_in_start = find(data.event_data(1,:,i) == IN_FIX_WIN_CD );
    target_begin = find(data.event_data(1,:,i) == TARGS_ON_CD );

    if (isempty(trial_begin|trial_over|target_begin|stim_begin|stim_over) == 1)
        i;
        %         trial_start
        %         trial_end
        %         fix_pt_on
        %         stim_start
        %         stim_end
        %         monk_in_start_window
        trial_over = 5000;
        trial_begin = 500;
%         fix_pt_begin = 10;
%         stim_begin = 991;
%         stim_over = 30*timebin; % Just so the min bin is not 1
%         monk_in_start = 0;
%         %         trial_durations(i) = 0;
%         target_begin = 10*timebin; % Just so the min bin is not 1
%         spike_data_adj(:,i) = NaN; % Get rid of the "bad data"
    end


    trial_start(i)= trial_begin;
    trial_end(i)=trial_over;
%     fix_pt_on(i)=fix_pt_begin;
    stim_start(i)=stim_begin;
    stim_end(i)=stim_over;
%     monk_in_start_window(i)=monk_in_start;
    target_onset(i) = target_begin;

    % Find the bin in which target_onset and saccade_begin
    %     saccade_onset(i) = saccade_beg;
    stim_onset_bin(i) = floor(stim_start(i)/timebin) + 1;
    target_onset_bin(i) = floor(target_onset(i)/timebin) + 1;
    saccade_onset_bin(i) = floor(stim_end(i)/timebin) + 1;

    RT(i) = (stim_end(i) - stim_start(i)); %RT for each trial

%     %     Plot saccades:
%     figure(2)
%     plot(data.eye_data(1,floor(fix_pt_on(i)/5):floor(trial_end(i)/5),i), data.eye_data(2,floor(fix_pt_on(i)/5):floor(trial_end(i)/5),i))
%     pause

    %--------------------------------------------------------------------------
    %     start_to_fix(i) = (trial_start(i) - fix_pt_on(i));
    %     trial_durations(i) = (trial_end(i) - trial_start(i));
    %     monk_fixates(i) = (monk_in_start_window(i) - fix_pt_on(i));
    %     fixation_to_target(i) = (target_onset(i) - monk_in_start_window(i));
    %     stim_durations(i) = (stim_end(i) - stim_start(i));
    %     pre_stim_delay(i) = (stim_start(i) - monk_in_start_window(i));
    %--------------------------------------------------------------------------

    % 2. Replace everything before and after trial start/end with NaNs

    %More efficient way, but isn't working with selected trials. Will fix
    %later. Until then use the alternate way below.
            spike_data(1:(trial_start-1),i) = NaN;
            spike_data((trial_end+1):end,i) = NaN;
    
            spike_data_adj = spike_data;

    %     % Alternate way...
%     for j = 1:data_duration
%         %         spike_data_matrix = NaN([total_trials,5000]);
%         if ((j < trial_start(i)) || (j > trial_end(i)))
%             spike_data_adj(j,i) = NaN;
%         else
%             spike_data_adj(j,i) = spike_data(j,i);
%         end
%     end
end


% monkey's choice - Left(2) or Right(1)
for i= 1 : length(total_trials)
    %     choice(i) = find(data.event_data(1,:,i) == IN_T1_WIN_CD);
    if (sum(data.event_data(1,:,i) == IN_T1_WIN_CD)>0)
        choice(i) = 1; %Right
    else
        choice(i) = 2; %Left
    end
end

% 3. Choose only the correct trials. Replace all wrong trials with NaN values
% that will be averaged out later using eNaN. (Include all 0 heading
% trials)
for i = 1:length(total_trials)
    if (heading(i) < 0 && choice(i) == 2)
        spike_data_adj(:,i) = spike_data_adj(:,i);
    elseif (heading(i) > 0 && choice(i) == 1)
        spike_data_adj(:,i) = spike_data_adj(:,i);
    elseif (heading(i) == 0)
        spike_data_adj(:,i) = spike_data_adj(:,i);
    else
        spike_data_adj(:,i) = NaN;
    end
end


for k=1: length(unique_stim_type)
    for i=1:length(unique_heading)

        select = logical( (heading==unique_heading(i)) & (stim_type==unique_stim_type(k)) );
        act_found = find( select==1 );

        targ_spike_bins{i,k}(:,:) = NaN(length(num_time_bins), length(act_found));
        targ_bin_start = min(target_onset_bin); %align all data to this bin number later
        targ_bin_end = max(target_onset_bin);

        sacc_spike_bins{i,k}(:,:) = NaN(length(num_time_bins), length(act_found));
        sacc_bin_start = min(saccade_onset_bin); %align all data to this bin number later
        sacc_bin_end = max(saccade_onset_bin);

        stim_spike_bins{i,k}(:,:) = NaN(length(num_time_bins), length(act_found)); % No need to realign... Stim on is aligned in the raw data
        stim_bin_start = min(stim_onset_bin); %align all data to this bin number later
        stim_bin_end = max(stim_onset_bin);


        for l = 1:length(act_found)
            % 4. Spike count in __ms time bins in each trial
            for j = 1:length(num_time_bins)
                spike_rate_bins{i,k}(j,l) = (sum(spike_data_adj(timebin*(j-1)+1:timebin*j,act_found(l))))/(timebin/1000); %converted to spikes/s
            end
            % 5. Realign data by target onset and saccade onset
            targ_spike_bins{i,k}(1:targ_bin_start,l) = spike_rate_bins{i,k}((target_onset_bin(act_found(l))-targ_bin_start)+1:target_onset_bin(act_found(l)),l);
            targ_spike_bins{i,k}(targ_bin_start+1:targ_bin_start+length(num_time_bins)-targ_bin_end,l) = spike_rate_bins{i,k}(target_onset_bin(act_found(l))+1:target_onset_bin(act_found(l))+length(num_time_bins)-targ_bin_end,l);

            sacc_spike_bins{i,k}(1:sacc_bin_start,l) = spike_rate_bins{i,k}((saccade_onset_bin(act_found(l))-sacc_bin_start)+1:saccade_onset_bin(act_found(l)),l);
            sacc_spike_bins{i,k}(sacc_bin_start+1:sacc_bin_start+length(num_time_bins)-sacc_bin_end,l) = spike_rate_bins{i,k}(saccade_onset_bin(act_found(l))+1:saccade_onset_bin(act_found(l))+length(num_time_bins)-sacc_bin_end,l);
            
        %%% Instead, use pdf for smoothing %%%
            gauss2 = normpdf(1:100,50,10); % smooth data at a SD of 50ms (10 points. Probably too much)
            spike_rate_pdf_temp = conv(spike_data_adj(:,act_found(l)), gauss2);
            spike_rate_pdf{i,k}(:,l) = spike_rate_pdf_temp(50:end-50);
            
            % Align the smoothed data to target and saccade onset
            targ_spike{i,k}(:,l) = spike_rate_pdf{i,k}(
            
        end

        % Seperate dataset of just long RTs
        find_long_RT = logical(RT(act_found) < 750 );
        long_RT{i,k} = find( find_long_RT == 1 );

        for m = 1:length(long_RT{i,k})
            long_targ_spike_bins{i,k}(:,m) = targ_spike_bins{i,k}(:,long_RT{i,k}(:,m));
            long_sacc_spike_bins{i,k}(:,m) = sacc_spike_bins{i,k}(:,long_RT{i,k}(:,m));
            long_stim_spike_bins{i,k}(:,m) = spike_rate_bins{i,k}(:,long_RT{i,k}(:,m));
        end


        % 6. Average the spike rates over trials
        targ_average_spikes{i,k} = nanmean(targ_spike_bins{i,k}(:,:),2);
        sacc_average_spikes{i,k} = nanmean(sacc_spike_bins{i,k}(:,:),2);
        stim_average_spikes{i,k} = nanmean(spike_rate_bins{i,k}(:,:),2);

        %--- Smoothing filter - weighted moving average
        a=1;
        b=[1/2 1/4 1/8 1/8];

        %--- Average of spike rates for long RTs
        %       long_targ_average_spikes{i,k} = filter(b,a,nanmean(long_targ_spike_bins{i,k}(:,:),2));
        %       long_sacc_average_spikes{i,k} = filter(b,a,nanmean(long_sacc_spike_bins{i,k}(:,:),2));
        %       long_stim_average_spikes{i,k} = filter(b,a,nanmean(long_stim_spike_bins{i,k}(:,:),2));

        smooth_targ_average_spikes{i,k} = filter(b,a,targ_average_spikes{i,k});
        smooth_sacc_average_spikes{i,k} = filter(b,a,sacc_average_spikes{i,k});
        smooth_stim_average_spikes{i,k} = filter(b,a,stim_average_spikes{i,k});


    end

end



%----- Averaging the middle headings,.25 and .5 together.
% for i = length(unique_heading)
%     if i == 4
%         x1 = [targ_average_spikes{4,1}, targ_average_spikes{5,1}];
%         x2 = [targ_average_spikes{4,2}, targ_average_spikes{5,2}];
%         x3 = [targ_average_spikes{4,3}, targ_average_spikes{5,3}];
%         targ_average_spikes{i,1} = nanmean(x1,2);
%         targ_average_spikes{i,2} = nanmean(x2,2);
%         targ_average_spikes{i,3} = nanmean(x3,2);
%
%         x1 = [sacc_average_spikes{4,1}, sacc_average_spikes{5,1}];
%         x2 = [sacc_average_spikes{4,2}, sacc_average_spikes{5,2}];
%         x3 = [sacc_average_spikes{4,3}, sacc_average_spikes{5,3}];
%         sacc_average_spikes{i,1} = nanmean(x1,2);
%         sacc_average_spikes{i,2} = nanmean(x2,2);
%         sacc_average_spikes{i,3} = nanmean(x3,2);
%
%         x1 = [stim_average_spikes{4,1}, stim_average_spikes{5,1}];
%         x2 = [stim_average_spikes{4,2}, stim_average_spikes{5,2}];
%         x3 = [stim_average_spikes{4,3}, stim_average_spikes{5,3}];
%         stim_average_spikes{i,1} = nanmean(x1,2);
%         stim_average_spikes{i,2} = nanmean(x2,2);
%         stim_average_spikes{i,3} = nanmean(x3,2);
%
%     elseif i == 7
%         x1 = [targ_average_spikes{6,1}, targ_average_spikes{7,1}];
%         x2 = [targ_average_spikes{6,2}, targ_average_spikes{7,2}];
%         x3 = [targ_average_spikes{6,3}, targ_average_spikes{7,3}];
%         targ_average_spikes{i,1} = nanmean(x1,2);
%         targ_average_spikes{i,2} = nanmean(x2,2);
%         targ_average_spikes{i,3} = nanmean(x3,2);
%
%         x1 = [sacc_average_spikes{6,1}, sacc_average_spikes{7,1}];
%         x2 = [sacc_average_spikes{6,2}, sacc_average_spikes{7,2}];
%         x3 = [sacc_average_spikes{6,3}, sacc_average_spikes{7,3}];
%         sacc_average_spikes{i,1} = nanmean(x1,2);
%         sacc_average_spikes{i,2} = nanmean(x2,2);
%         sacc_average_spikes{i,3} = nanmean(x3,2);
%
%         x1 = [stim_average_spikes{6,1}, stim_average_spikes{7,1}];
%         x2 = [stim_average_spikes{6,2}, stim_average_spikes{7,2}];
%         x3 = [stim_average_spikes{6,3}, stim_average_spikes{7,3}];
%         stim_average_spikes{i,1} = nanmean(x1,2);
%         stim_average_spikes{i,2} = nanmean(x2,2);
%         stim_average_spikes{i,3} = nanmean(x3,2);
%     end
% end
%-----

%Analyse data for 300-600ms after stimulus onset
for h = 1: length(unique_heading)
    for k = 1:length(unique_stim_type)
        stim_firing(h,k) = nanmean(stim_average_spikes{h,k}(stim_bin_start-floor(300/timebin):stim_bin_start+floor(600/timebin)));
    end
end

% PLOT
for h = 1: length(unique_heading)
    fig = [0 0 0 0 0 0 2 4 6 8]; % temporarily hard-coded for 10 directions
    if unique_heading(h) < 0 % Left (out of RF for now)
        figure (h+1)
        orient landscape;
        axis tight;
        for k = 1:length(unique_stim_type)
            if k == 1
                color = 'b';
                color2 = 'b:';
            elseif k == 2
                color = 'r';
                color2 = 'r:';
            else
                color = 'g';
                color2 = 'g:';
            end
            
            % Plot smoothed out data
            subplot(2,3,1), plot(timebin*(-targ_bin_start+1:floor(500/timebin)),smooth_targ_average_spikes{h,k}(1:targ_bin_start+floor(500/timebin),1), color,'LineWidth',2) %targ_bin_start before target and 300ms after target on
            hold on
            subplot(2,3,2), plot(timebin*(-floor(400/timebin):floor(800/timebin)),smooth_stim_average_spikes{h,k}(stim_bin_start-floor(400/timebin):stim_bin_start+floor(800/timebin),1), color, 'LineWidth',2)
            title ('out of RF');
            hold on
            subplot(2,3,3), plot(timebin*(-floor(400/timebin):floor(300/timebin)),smooth_sacc_average_spikes{h,k}(sacc_bin_start-floor(400/timebin):sacc_bin_start+floor(300/timebin),1), color,'LineWidth',2)
            title (unique_heading(h));
            hold on
            
%             %Plot each trial - Just to check.
%             for i = 1:length(act_found)
%             subplot(2,3,1), plot(timebin*(1:targ_bin_start+floor(300/timebin)),targ_spike_bins{h,k}(1:targ_bin_start+floor(300/timebin),i), color2, targ_bin_start*timebin,1:.5:50,'k') %targ_bin_start before target and 300ms after target on
%             hold on
%             subplot(2,3,2), plot(timebin*(stim_bin_start-200/timebin:stim_bin_start+floor(800/timebin)),spike_rate_bins{h,k}(stim_bin_start-floor(200/timebin):stim_bin_start+floor(800/timebin),i), color2, stim_bin_start*timebin,1:.5:50,'k')
%             hold on
%             subplot(2,3,3), plot(timebin*(sacc_bin_start-400/timebin:sacc_bin_start+floor(300/timebin)),sacc_spike_bins{h,k}(sacc_bin_start-floor(400/timebin):sacc_bin_start+floor(300/timebin),i), color2, sacc_bin_start*timebin,1:.5:50,'k')
%             hold on
%             end

%             % Plot average spikes (no filtering)
            %             subplot(2,3,1), plot(timebin*(1:75),targ_average_spikes{h,k}(1:75,1), color, targ_bin_start*timebin,1:.5:100,'k') %targ_bin_start before target and 300ms after target on
            %             title ('out of RF');
            %             xlabel('time (ms)');
            %             hold on
            %             subplot(2,3,2), plot(timebin*(1:75),stim_average_spikes{h,k}(1:75,1), color, stim_onset_bin*timebin,1:.5:100,'k')
            %             xlabel('time (ms)');
            %             hold on
            %             subplot(2,3,3), plot(timebin*(1:75),sacc_average_spikes{h,k}(1:75,1), color, sacc_bin_start*timebin,1:.5:100,'k')
            %             xlabel('time (ms)');
            %             title (unique_heading(h));
            %             hold on
        end
    end

    if unique_heading(h) > 0 % Right (in RF for now)
        figure (h-fig(h))
        for k = 1:length(unique_stim_type)
            if k == 1
                color = 'b';
                color2 = 'b:';
            elseif k == 2
                color = 'r';
                color2 = 'r:';
            else
                color = 'g';
                color2 = 'g:';
            end

            subplot(2,3,4), plot(timebin*(-targ_bin_start+1:floor(500/timebin)),smooth_targ_average_spikes{h,k}(1:targ_bin_start+floor(500/timebin),1), color,'LineWidth',2) %targ_bin_start before target and 300ms after target on
            hold on
            subplot(2,3,5), plot(timebin*(-floor(400/timebin):floor(800/timebin)),smooth_stim_average_spikes{h,k}(stim_bin_start-floor(400/timebin):stim_bin_start+floor(800/timebin),1), color,'LineWidth',2)
            title ('into RF');
            hold on
            subplot(2,3,6), plot(timebin*(-floor(400/timebin):floor(300/timebin)),smooth_sacc_average_spikes{h,k}(sacc_bin_start-floor(400/timebin):sacc_bin_start+floor(300/timebin),1), color,'LineWidth',2)
            title (unique_heading(h));
            hold on
            
%             %Plot each trial - Just to check.
%             for i = 1:length(act_found)
%             subplot(2,3,4), plot(timebin*(1:targ_bin_start+floor(300/timebin)),targ_spike_bins{h,k}(1:targ_bin_start+floor(300/timebin),i), color2, targ_bin_start*timebin,1:.5:50,'k') %targ_bin_start before target and 300ms after target on
%             hold on
%             subplot(2,3,5), plot(timebin*(stim_bin_start-floor(200/timebin):stim_bin_start+800/timebin),spike_rate_bins{h,k}(stim_bin_start-floor(200/timebin):stim_bin_start+floor(800/timebin),i), color2, stim_bin_start*timebin,1:.5:50,'k')
%             hold on
%             subplot(2,3,6), plot(timebin*(sacc_bin_start-floor(400/timebin):sacc_bin_start+300/timebin),sacc_spike_bins{h,k}(sacc_bin_start-floor(400/timebin):sacc_bin_start+floor(300/timebin),i), color2, sacc_bin_start*timebin,1:.5:50,'k')
%             hold on
%             end
            
            %             subplot(2,3,4), plot(timebin*(1:75), targ_average_spikes{h,k}(1:75,1), color, targ_bin_start*timebin,1:.5:100,'k') %targ_bin_start before target and 300ms after target on
            %             title ('in the RF');
            %             xlabel('time (ms)');
            %             hold on
            %             subplot(2,3,5), plot(timebin*(1:75), stim_average_spikes{h,k}(1:75,1), color, stim_onset_bin*timebin,1:.5:100,'k')
            %             xlabel('time (ms)');
            %             hold on
            %             subplot(2,3,6), plot(timebin*(1:75), sacc_average_spikes{h,k}(1:75,1), color, sacc_bin_start*timebin,1:.5:100,'k')
            %             xlabel('time (ms)');
            %             title (unique_heading(h));
            %             hold on

        end
    end

    if unique_heading(h) == 0 % 0 heading
        figure (7)
        title (unique_heading(h));
        for k = 1:length(unique_stim_type)
            if k == 1
                color = 'b';
            elseif k == 2
                color = 'r';
            else
                color = 'g';
            end

            subplot(2,3,4), plot(timebin*(1:targ_bin_start+300/timebin),smooth_targ_average_spikes{h,k}(1:targ_bin_start+300/timebin,1), color, targ_bin_start*timebin,1:.5:100,'k') %targ_bin_start before target and 300ms after target on
            title ('into RF');
            hold on
            subplot(2,3,5), plot(timebin*(stim_bin_start-200/timebin:stim_bin_start+800/timebin),smooth_stim_average_spikes{h,k}(stim_bin_start-200/timebin:stim_bin_start+800/timebin,1), color, stim_bin_start*timebin,1:.5:100,'k')
            hold on
            subplot(2,3,6), plot(timebin*(sacc_bin_start-400/timebin:sacc_bin_start+300/timebin),smooth_sacc_average_spikes{h,k}(sacc_bin_start-400/timebin:sacc_bin_start+300/timebin,1), color, sacc_bin_start*timebin,1:.5:100,'k')
            title (unique_heading(h));
            hold on

            %             subplot(2,3,4), plot(timebin*(1:75), targ_average_spikes{h,k}(1:75,1), color, targ_bin_start*timebin,1:.5:100,'k') %targ_bin_start before target and 300ms after target on
            %             title ('in the RF');
            %             xlabel('time (ms)');
            %             hold on
            %             subplot(2,3,5), plot(timebin*(1:75), stim_average_spikes{h,k}(1:75,1), color, stim_onset_bin*timebin,1:.5:100,'k')
            %             xlabel('time (ms)');
            %             hold on
            %             subplot(2,3,6), plot(timebin*(1:75), sacc_average_spikes{h,k}(1:75,1), color, sacc_bin_start*timebin,1:.5:100,'k')
            %             xlabel('time (ms)');timebin

            %             title (unique_heading(h));
            %             hold on

        end
    end


end

%---------------------------------------------------------------------------------------
return;