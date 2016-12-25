%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_1D_Az_VF_headeye.m -- For MOOG 1D Vary Fixation expt with eye and head:
%-- Plots response as a function of azimuth for each stimulus condition and each eye/head position  CRF 03/2010
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_1D_Az_VF_headeye(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% time = clock;
% disp('START TIME = ');
% disp(time(4:6));
% tic

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
temp_head_fix_x = data.moog_params(HEAD_FIX_X,:,MOOG);
temp_head_fix_y = data.moog_params(HEAD_FIX_Y,:,MOOG);
%now, get the firing rates and spike data for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             
temp_spike_data = data.spike_data(SpikeChan, :);
% temp_event_data = data.event_data(SpikeChan, :);
temp_event_data = squeeze(data.event_data);

%--------------------------------------------------------------------------
% Based on a reviewer comment, this will use a 50 ms analysis window and
% re-compute DI at multiple intervals throughout the stimulus period.  (CRF, 9/2006)
%--------------------------------------------------------------------------

% for bin_num = 1:40  % the END of this loop is at the very end of the script

% bin_num
% clear temp_spike_rates
% for a = BegTrial:EndTrial
%     bin_incr = ((a-1)*5000 + 996 + (bin_num-1)*50 + 1);
%     if bin_num == 40
%         temp_spike_rates(a) = (sum(temp_spike_data( bin_incr : bin_incr + 50 ))) * 20;
%     else
%         temp_spike_rates(a) = (sum(temp_spike_data( bin_incr : bin_incr + 100 ))) * 10;
%     end
% end

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
spon_found = find(null_trials==1); 
%%spon_resp = mean(temp_spike_rates(spon_found));
%     This won't work -- evidently there are no null trials in this protocol
%     Instead, compute spontaneous from pre-trial activity -- CRF 3-5-10
for i = 1:EndTrial-BegTrial+1
    trial_event_times{i} = find(temp_event_data(:,i)>0);
    trial_event_codes{i} = temp_event_data(trial_event_times{i},i);
    if trial_event_codes{i}(1) == 1
        trial_start_times(i) = trial_event_times{i}(1);
    else  % sometimes (rarely) the 'start trial' event is missing, so just manually pick 500 for those
        trial_start_times(i) = 500;
    end
    trial_spon_rate(i) = sum(temp_spike_data( i*5000-4999 : i*5000-(5000-trial_start_times(i)) )) / (trial_start_times(i)/1000);
%     trial_spon_rate(i) = sum(temp_spike_data( (BegTrial+i-1)*5000-4999 : (BegTrial+i-1)*5000-(5000-trial_start_times(i)) )) / (trial_start_times(i)/1000); ?????
end
spon_resp = mean(trial_spon_rate);

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_trials ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
fix_x     = temp_fix_x(~null_trials & select_trials);
fix_y     = temp_fix_y(~null_trials & select_trials);
head_fix_x = temp_head_fix_x(~null_trials & select_trials);
head_fix_y = temp_head_fix_y(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 9999;
end
spike_data(1,:) = temp_spike_data( 1, find(temp_spike_data(1,:)~=9999) );
spike_data(1, find(spike_data>10) ) = 1; % something is absolutely wrong 

unique_azimuth  = munique(azimuth');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_fix_x    =  munique(fix_x');
unique_fix_y    =  munique(fix_y');
unique_head_fix_x    =  munique(head_fix_x');
unique_head_fix_y    =  munique(head_fix_y');


if length(unique_fix_y) == 1
   condition = fix_x;
   temp_condition = temp_fix_x;
   unique_condition = unique_fix_x;
else
   condition = fix_y; 
   temp_condition = temp_fix_y;
   unique_condition = unique_fix_y;
end
if length(unique_head_fix_y) == 1
   condition_2 = head_fix_x;
   temp_condition_2 = temp_head_fix_x;
   unique_condition_2 = unique_head_fix_x;
else
   condition_2 = head_fix_y; 
   temp_condition_2 = temp_head_fix_y;
   unique_condition_2 = unique_head_fix_y;
end

% add titles
titles{1} = 'Vestibular, ';
titles{2} = 'Visual, ';
for n=1: length(unique_stim_type)
    for k=1: length(unique_condition)
        h_title{k,n} = [titles{unique_stim_type(n)}, num2str(unique_condition(k))];
    end
end

% number of vectors in each condition 
vector_num = length(unique_azimuth);
% repetition = floor( length(spike_rates) / (length(unique_azimuth)*length(unique_condition)*length(unique_condition_2)*length(unique_stim_type)) );
% ^^ doesn't work because not all combinations of eye and head are tested. Kluge for now:
repetition = floor( length(spike_rates) / (length(unique_azimuth)*5*length(unique_stim_type)) );
vari_straight = (1+length(unique_condition))/2;


% for Xiaodong: below here, must first change all loops to cycle through
% each combination of head and eye position, which is NOT
% length(unique_condition) * length(unique_condition_2)...


% response matrix
unique_azimuth_8=[0,45,90,135,180,225,270,315];
StartEventBin(1)=996;
windowlength = 50; % 50 ms
for w = 1:1
windowlength = windowlength + (w-1)*50;
%count = 0;
for s = 1:1
    % replace spike_rates with spike_data based on analize window set    
%     for ss =  1 : length(spike_rates) % ss marks the index of trial
% %         spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+windowlength*(s-1)+5000*(ss-1) : StartEventBin(1)+windowlength+windowlength*(s-1)+5000*(ss-1)) ) ; % 996~3006 every 200ms
%           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+900+5000*(ss-1) : StartEventBin(1)+900+200+5000*(ss-1)) );%Enlarge the window
%     end    
    resp_mat = [];
    resp_mat_std = []; 
    count = 0;
%    for k=1:length(unique_condition)
    for k=1:1        
        for n=1:length(unique_stim_type)
            count = 0;
         %   for i=1:length(unique_azimuth_8)    
             for i=1:length(unique_azimuth) 
        %      select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
         %       select = logical( (azimuth==unique_azimuth_8(i)) & (condition==unique_condition(vari_straight)) & (stim_type==unique_stim_type(n)) );
                select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(vari_straight)) & (stim_type==unique_stim_type(n)) );
                if (sum(select) > 0) 
                    count = count+1;
                    spike_temp = spike_rates(select);
                    resp_mat_anova_horizontal{n}(1:repetition,i) = spike_temp(1:repetition);
                    
                    resp_mat(k, n, i) = mean(spike_rates(select));
                    resp_mat_std(k, n, i) = std(spike_rates(select));    % calculate std between trials for later HTI usage                
                    resp_mat_count(n,i) = resp_mat(1, n, i);
                    resp_mat_std_count(n,i) = resp_mat_std(1, n, i);
                    
                    % z-score data for spike count correlation analysis
                    z_dist = spike_rates(select);
                    if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                       z_dist = (z_dist - mean(z_dist))/std(z_dist);
                    else
                        z_dist = 0;
                    end
                    Z_Spikes(select) = z_dist;  
            
                else
                    resp_mat(k, n, i) = 0;                
                    resp_mat_std(k, n, i) = 0;
                end
            end     
            [p_anova, table, stats] = anova1(resp_mat_anova_horizontal{n},[],'off');
            P_anova_horizontal(n) = p_anova;
        end
    end

    % create 'wrap-around' version of resp_mat such that the center of the plot 
    % will correspond to forward and both lateral edges to backward
    resp_mat_tran = [];
    resp_mat_std_tran = [];
    for t = 1:length(unique_azimuth)  
        if t < 3       %%%%%(not completely general)
            resp_mat_tran(:,:,t) = resp_mat(:,:,t + (length(unique_azimuth) - 2));
            resp_mat_std_tran(:,:,t) = resp_mat_std(:,:,t + (length(unique_azimuth) - 2));
            unique_az_tran(t) = unique_azimuth(t + (length(unique_azimuth) - 2));
        else
            resp_mat_tran(:,:,t) = resp_mat(:,:,t-2);
            resp_mat_std_tran(:,:,t) = resp_mat_std(:,:,t-2);
            unique_az_tran(t) = unique_azimuth(t-2);
        end
    end
    resp_mat_tran(:,:,end+1) = resp_mat_tran(:,:,1);
    resp_mat_std_tran(:,:,end+1) = resp_mat_std_tran(:,:,1);
    unique_az_tran(end+1) = unique_az_tran(1);
    unique_az_tran = unique_az_tran';
    
    % calculate spontaneous firing rate -- CORRECTED (now only includes Begtrial:Endtrial)
    temp_spike_rates_true = temp_spike_rates(BegTrial:EndTrial);
    for k= 1:length(unique_condition)
        spon_resp(k) = mean( temp_spike_rates_true( find( null_trials(BegTrial:EndTrial)==1 & temp_condition(BegTrial:EndTrial)==unique_condition(k) ) ) );
    end
    
    trials_per_rep = (length(unique_azimuth)*length(unique_condition)*length(unique_stim_type)+3);
    repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
    
    % % %------------------------------------------------------------------
    % % %--------------TEMP - OUTPUT BASIC TUNING CURVE DATA---------------
    % % %------------------------------------------------------------------
    % 
    % % response matrix
    % resp_mat_basic = [];
    % resp_mat_ste_basic = [];
    % for k=1:length(unique_condition)
    % 	for n=1:length(unique_stim_type)
    %         for i=1:length(unique_azimuth)
    %             select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
    %             if (sum(select) > 0)                
    %                 resp_mat_basic(i,n,k) = mean(spike_rates(select));
    %                 resp_mat_ste_basic(i,n,k) = std(spike_rates(select))./sqrt(repetitions);    % calculate std between trials for later HTI usage                
    %             else
    %                 resp_mat_basic(i,n,k) = 0;                
    %                 resp_mat_ste_basic(i,n,k) = 0;
    %             end
    %         end        
    % 	end
    % end
    % 
    % for i=1:length(unique_azimuth)
    %     buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', ...
    %     FILE, unique_azimuth(i), resp_mat_basic(i,1,1), resp_mat_basic(i,2,1), resp_mat_basic(i,3,1), resp_mat_ste_basic(i,1,1), resp_mat_ste_basic(i,2,1), resp_mat_ste_basic(i,3,1), resp_mat_basic(i,1,2), resp_mat_basic(i,2,2), resp_mat_basic(i,3,2), resp_mat_ste_basic(i,1,2), resp_mat_ste_basic(i,2,2), resp_mat_ste_basic(i,3,2), resp_mat_basic(i,1,3), resp_mat_basic(i,2,3), resp_mat_basic(i,3,3), resp_mat_ste_basic(i,1,3), resp_mat_ste_basic(i,2,3), resp_mat_ste_basic(i,3,3));
    %     outfile = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\justin_cong-opp.dat'];
    %     printflag = 0;
    %     if (exist(outfile, 'file') == 0)    %file does not yet exist
    %         printflag = 1;
    %     end
    % 	fid = fopen(outfile, 'a');
    % 	if (printflag)
    %         fprintf(fid, 'FILE\t Azimuth\t Ves_mean_-22.5\t Vis_mean_-22.5\t Comb_mean_-22.5\t Ves_StErr_-22.5\t Vis_StErr_-22.5\t Comb_StErr_-22.5\t Ves_mean_0\t Vis_mean_0\t Comb_mean_0\t Ves_StErr_0\t Vis_StErr_0\t Comb_StErr_0\t Ves_mean_+22.5\t Vis_mean_+22.5\t Comb_mean_+22.5\t Ves_StErr_+22.5\t Vis_StErr_+22.5\t Comb_StErr_+22.5\t');
    %         fprintf(fid, '\r\n');
    % 	end
    % 	fprintf(fid, '%s', buff);
    % 	fprintf(fid, '\r\n');
    % 	fclose(fid);
    % end
    % 
    % return;
    % 
    % % %------------------------------------------------------------------
    % % %--------------END OUTPUT BASIC TUNING CURVE DATA------------------
    % % %------------------------------------------------------------------
    
    
    % % x_azimuth = [1 2 3 4 4.5 5 5.5 6 7 8 9];
    % x_azimuth = [270 225 180 135 112.5 90 67.5 45 0 -45 -90];
    % 
    % data_ves(:,1) = x_azimuth';
    % data_ves(:,2) = squeeze(resp_mat_tran(1,1,:));
    % data_ves(:,3) = squeeze(resp_mat_std_tran(1,1,:)/sqrt(repetitions));
    % data_ves(:,4) = squeeze(resp_mat_tran(2,1,:));
    % data_ves(:,5) = squeeze(resp_mat_std_tran(2,1,:)/sqrt(repetitions));
    % data_ves(:,6) = squeeze(resp_mat_tran(3,1,:));
    % data_ves(:,7) = squeeze(resp_mat_std_tran(3,1,:)/sqrt(repetitions));
    % 
    % data_vis(:,1) = x_azimuth';
    % data_vis(:,2) = squeeze(resp_mat_tran(1,2,:));
    % data_vis(:,3) = squeeze(resp_mat_std_tran(1,2,:)/sqrt(repetitions));
    % data_vis(:,4) = squeeze(resp_mat_tran(2,2,:));
    % data_vis(:,5) = squeeze(resp_mat_std_tran(2,2,:)/sqrt(repetitions));
    % data_vis(:,6) = squeeze(resp_mat_tran(3,2,:));
    % data_vis(:,7) = squeeze(resp_mat_std_tran(3,2,:)/sqrt(repetitions));
    % 
    % data_comb(:,1) = x_azimuth';
    % data_comb(:,2) = squeeze(resp_mat_tran(1,3,:));
    % data_comb(:,3) = squeeze(resp_mat_std_tran(1,3,:)/sqrt(repetitions));
    % data_comb(:,4) = squeeze(resp_mat_tran(2,3,:));
    % data_comb(:,5) = squeeze(resp_mat_std_tran(2,3,:)/sqrt(repetitions));
    % data_comb(:,6) = squeeze(resp_mat_tran(3,3,:));
    % data_comb(:,7) = squeeze(resp_mat_std_tran(3,3,:)/sqrt(repetitions));
    
    
    % calculate min and max firing rate, standard deviation, HTI, Vectorsum
    for n=1:length(unique_stim_type)
        for k=1: length(unique_condition)
            % calculate min and max firing rate, standard deviation, HTI, Vectorsum
            Min_resp(k,n) = min( resp_mat_tran(k,n,:) );
            Max_resp(k,n) = max( resp_mat_tran(k,n,:) );
            resp_std(k,n) = sum( sum(resp_mat_std(k,n,:)) ) / vector_num; % average standard deviation?
            M = squeeze(resp_mat(k,n,:));
            M = M';
            M(3) = [];  % for now, remove the extra sampled points at front of sphere
            M(4) = [];  % (can't do vectorsum and still keep them in there)
            [Azi, Ele, Amp] = vectorsum(M);
            Vec_sum{k,n} = [Azi, Ele, Amp];
            Vec_sum_azi(k,n) = Azi;
            HTI_temp(k,n) = HTI(M,spon_resp(k));
        end
    end
    
    % -------------------------------------------------------------------------
    % ANOVA to test tuning significance
    
    % first parse raw data into repetitions, including null trials
    for q = 1:repetitions
        azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        stim_type_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        condition_rep{q} = temp_condition(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        spike_data_rep{q} = temp_spike_data( (BegTrial-1)*5000 + ((q-1)*trials_per_rep*5000) + 1 : q*trials_per_rep*5000 + (BegTrial-1)*5000 );
        event_data_rep{q} = temp_event_data( (BegTrial-1)*5000 + ((q-1)*trials_per_rep*5000) + 1 : q*trials_per_rep*5000 + (BegTrial-1)*5000 );
    end
    
    for n = 1:length(unique_stim_type)
        for k = 1:length(unique_condition)
            for i = 1:length(unique_azimuth)
                clear select_rep;
                for q = 1:repetitions
                    select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & condition_rep{q}==unique_condition(k) & stim_type_rep{q}==unique_stim_type(n) );
                    resp_anova{k,n}(q,i) = spike_rates_rep{q}(select_rep{q});
                end
            end
            [p, table, stats] = anova1(resp_anova{k,n}(:,:),[],'off');
            p_anova(k,n) = p;  anova_table{k,n} = table;  anova_stats{k,n} = stats;
            F_val(k,n) = anova_table{k,n}(2,5);
        end
    end
    F_val = cell2mat(F_val);
    
    % -------------------------------------------------------------------------
    % Gain Field test (ANOVA/K-W of individual trial data across gaze angles)
    % Added GF slope and linear regression r/p values - crf 11/05
    % Added monotonic test - crf 02/06
    
    for n=1:length(unique_stim_type)
        for k=1:length(unique_condition)
            % max and min of response matrix (means) at each stim/gaze condition
            whereis_max = logical(resp_mat == Max_resp(k,n));
            whereis_max = squeeze(whereis_max(k,n,:));
            max_azi{k,n} = unique_azimuth(find(whereis_max));
            whereis_min = logical(resp_mat == Min_resp(k,n));
            whereis_min = squeeze(whereis_min(k,n,:));
            min_azi{k,n} = unique_azimuth(find(whereis_min));
            if length(max_azi{k,n}) > 1   % need these to be unique, and believe it or not they sometimes aren't
                max_azi{k,n} = max_azi{k,n}(1);
            end
            if length(min_azi{k,n}) > 1
                min_azi{k,n} = min_azi{k,n}(1);
            end
    
            % now for each rep, calculate max-min and max-spontaneous
            clear q;
            for q = 1: length(spike_rates_rep)
                max_index = find(azimuth_rep{q}==max_azi{k,n} & condition_rep{q}==unique_condition(k) & stim_type_rep{q}==unique_stim_type(n));
                min_index = find(azimuth_rep{q}==min_azi{k,n} & condition_rep{q}==unique_condition(k) & stim_type_rep{q}==unique_stim_type(n));
                spon_index = find(stim_type_rep{q} == -9999 & condition_rep{q} == unique_condition(k));
                max_allreps{n}(q,k) = spike_rates_rep{q}(max_index);
                spon_allreps{n}(q,k) = spike_rates_rep{q}(spon_index);
                max_min_allreps{n}(q,k) = spike_rates_rep{q}(max_index) - spike_rates_rep{q}(min_index);
                max_spon_allreps{n}(q,k) = spike_rates_rep{q}(max_index) - spike_rates_rep{q}(spon_index);
                x(q,k) = unique_condition(k);
            end
        end
        
        % simple monotonicity test: if max-spon at zero gaze is sig > or < than at both
        % eccentric gaze angles, then GF is non-monotonic
        clear resp_minus resp_zero resp_plus;
        resp_minus = max_spon_allreps{n}(:,1); resp_zero = max_spon_allreps{n}(:,2); resp_plus = max_spon_allreps{n}(:,3);
        if ((mean(resp_zero) > mean(resp_plus) & mean(resp_zero) > mean(resp_minus)) | (mean(resp_zero) < mean(resp_plus) & mean(resp_zero) < mean(resp_minus))) & (ttest2(resp_zero, resp_plus) == 1 & ttest2(resp_zero, resp_minus) == 1)
            monotonic(n) = 0;
        else
            monotonic(n) = 1;
        end
        
        p_anova_max(n) = anova1(max_allreps{n},[],'off');
        p_anova_spon(n) = anova1(spon_allreps{n},[],'off');
        p_anova_maxmin(n) = anova1(max_min_allreps{n},[],'off');
        p_anova_maxspon(n) = anova1(max_spon_allreps{n},[],'off');
        gf_anovas{n} = [p_anova_max(n) p_anova_spon(n) p_anova_maxmin(n) p_anova_maxspon(n)];
        
        p_KW_max(n) = kruskalwallis(max_allreps{n},[],'off');
        p_KW_spon(n) = kruskalwallis(spon_allreps{n},[],'off');
        p_KW_maxmin(n) = kruskalwallis(max_min_allreps{n},[],'off');
        p_KW_maxspon(n) = kruskalwallis(max_spon_allreps{n},[],'off');
        gf_KWs{n} = [p_KW_max(n) p_KW_spon(n) p_KW_maxmin(n) p_KW_maxspon(n)];
    
        y = max_allreps{n};
        p = polyfit(x,y,1);
        f = polyval(p,x);
        [c,P] = corrcoef(f,y);
        GFr_max(n) = c(1,2);
        GFp_max(n) = P(1,2);
        GFslope_max(n) = p(1);
        
        y = spon_allreps{n};
        p = polyfit(x,y,1);
        f = polyval(p,x);
        [c,P] = corrcoef(f,y);
        GFr_spon(n) = c(1,2);
        GFp_spon(n) = P(1,2);
        GFslope_spon(n) = p(1);
        
        y = max_min_allreps{n};
        p = polyfit(x,y,1);
        f = polyval(p,x);
        [c,P] = corrcoef(f,y);
        GFr_maxmin(n) = c(1,2);
        GFp_maxmin(n) = P(1,2);
        GFslope_maxmin(n) = p(1);
        
        y = max_spon_allreps{n};
        p = polyfit(x,y,1);
        f = polyval(p,x);
        [c,P] = corrcoef(f,y);
        GFr_maxspon(n) = c(1,2);
        GFp_maxspon(n) = P(1,2);
        GFslope_maxspon(n) = p(1);
        
    end
    clear x y P f c;
    
    % %------------------------------------------------------------------
    % %----------------------TEMP - GF OUTPUT----------------------------
    % %------------------------------------------------------------------
    % 
    % % buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %3.0f\t', ...
    % %        FILE, p_anova(:), spon_resp, Min_resp(:), Max_resp(:), gf_anovas{:}, gf_KWs{:}, GFr_max, GFp_max, GFslope_max, GFr_spon, GFp_spon, GFslope_spon, GFr_maxmin, GFp_maxmin, GFslope_maxmin, GFr_maxspon, GFp_maxspon, GFslope_maxspon, EndTrial-(BegTrial-1));
    %        % currently 92 fields (11/05)
    % buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t', FILE, monotonic);
    % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\Dir1D_Az_VaryFix_Sum.dat'];
    % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\monotonic.dat'];
    % printflag = 0;
    % if (exist(outfile, 'file') == 0)    %file does not yet exist
    %     printflag = 1;
    % end
    % fid = fopen(outfile, 'a');
    % if (printflag)
    %     fprintf(fid, 'FILE\t Monotonic1\t Monotonic2\t Monotonic3\t');
    %     fprintf(fid, '\r\n');
    % end
    % fprintf(fid, '%s', buff);
    % fprintf(fid, '\r\n');
    % fclose(fid);
    % 
    % return;
    % %------------------------------------------------------------------
    % %----------------------TEMP - GF OUTPUT----------------------------
    % %------------------------------------------------------------------
    
    
    % % -------------------------------------------------------------------------
    % % Response Latency (sliding ANOVA method -- Avillac et al. 2005)
    % 
    % clear baseline_spikes stim_spikes latency;
    % for n = 1:length(unique_stim_type)
    %     for k = 1:length(unique_condition)
    %         for q = 1:repetitions
    %             clear max_index max_events max_spikes;
    %             max_index = find(azimuth_rep{q}==max_azi{k,n} & condition_rep{q}==unique_condition(k) & stim_type_rep{q}==unique_stim_type(n));
    %             max_events = event_data_rep{q}(((max_index-1)*5000)+1 : max_index*5000);
    %             max_spikes = spike_data_rep{q}(((max_index-1)*5000)+1 : max_index*5000);
    %         % for baseline, take 100 ms window before stimulus onset (event code '4')
    %             baseline_spikes{k,n}(q) = sum(max_spikes(find(max_events==4)-100 : find(max_events==4)-1));
    %             stim_spikes{k,n}(q,:) = max_spikes(find(max_events==4) : end);
    %         end
    %         % now compare baseline rate to rate in a sliding 20 ms window moving through 'stim_spikes'
    %         p_anova_latency = 1;
    %         T = 1;
    %         while p_anova_latency >= 0.05 | mean((sum(stim_spikes{k,n}(:,T:19+T)')'*50)) <= mean((baseline_spikes{k,n}')*10)
    %             p_anova_latency = anova1([sum(stim_spikes{k,n}(:,T:19+T)')'*50 (baseline_spikes{k,n}')*10], [] , 'off');
    %             T = T + 1;
    %             if T == 1981
    %                 break;
    %             end
    %             
    % %             if floor(T/20) == T/20
    % %                 T
    % %                 [sum(stim_spikes{k,n}(:,T:19+T)')'*50 (baseline_spikes{k,n}')*10]
    % %                 p_anova_latency
    % %                 pause;
    % %             end            
    %             
    %         end
    %         latency(k,n) = T + 19
    %     end
    % end
    
    
    % %------------------------------------------------------------------
    % %----------------------TEMP - LATENCY OUTPUT-----------------------
    % %------------------------------------------------------------------
    % 
    % buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', FILE, latency);
    % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\Dir1D_Az_VaryFix_Latency.dat'];
    % printflag = 0;
    % if (exist(outfile, 'file') == 0)    %file does not yet exist
    %     printflag = 1;
    % end
    % fid = fopen(outfile, 'a');
    % if (printflag)
    %     fprintf(fid, 'FILE\t Latency_minus\t Latency_zero\t Latency_plus\t Latency_minus2\t Latency_zero2\t Latency_plus2\t Latency_minus3\t Latency_zero3\t Latency_plus3\t');
    %     fprintf(fid, '\r\n');
    % end
    % fprintf(fid, '%s', buff);
    % fprintf(fid, '\r\n');
    % fclose(fid);
    % 
    % return;
    % %------------------------------------------------------------------
    % %----------------------TEMP - LATENCY OUTPUT-----------------------
    % %------------------------------------------------------------------
    
    
    %------------------------------------------------------------------
    % Define figure
    
    figure(2);
    orient landscape;
    set(2,'Position', [5,25 1200,900], 'Name', '1D Azimuth Tuning VaryFix');
    axis off;
    for n=1:length(unique_stim_type)
        for k=1: length(unique_condition)
            xoffset = (n-1) * 0.33;
            yoffset = (k-1) * -0.28;
            axes('position',[0.02+xoffset 0.62+yoffset 0.28 0.22]);
            x_azimuth = [1 2 3 4 4.5 5 5.5 6 7 8 9];
            errorbar(x_azimuth, resp_mat_tran(k,n,:), resp_mat_std_tran(k, n, :)/sqrt(repetitions), 'ko-');
            hold on;
            spon_line(1:length(x_azimuth)) = spon_resp(k);  % display spontaneous as a dotted line
            plot(x_azimuth,spon_line,'b:');
            xlim( [x_azimuth(1), x_azimuth(end)] );
            set(gca, 'xtick', x_azimuth );
            set(gca, 'xdir' , 'reverse');
            set(gca, 'xticklabel','-90|-45|0|45| |90| |135|180|225|270');
            if k == length(unique_condition)
                xlabel('Azimuth');
            end
            title( h_title{k,n} );
        end
    end
    
    % close;
    
    proceed = 0;
    proceed = input('proceed with DI, HDI, & Congruency?');
    if proceed == 1
    
        
    %----------------------------------------------------------
    % Cross-covariance method for finding shift ratio ('displacement index', Avillac et al. 2005)
    
    bin = 1; % interval (in degrees) between interpolated points -- 360/bin must be an even integer
    method = 'linear';
    showplots = 0;
    x_values = [0 45 67.5 90 112.5 135 180 225 270 315 360];
    DI{bin_num} = cross_covariance(resp_mat, unique_stim_type, unique_condition, bin, method, showplots, x_values);
    % DI = zeros(3,3);
    
    
    % %----------------------------------------------------------
    % % Congruency metric: simple correlation between vis+ves (interpolated) tuning curves
    % 
    % xi = 0 : bin : 359;
    % resp_mat_360 = resp_mat;
    % for n = 1:length(unique_stim_type)
    %     for k = 1:length(unique_condition)
    %         resp_mat_360(k,n,length(x_values)) = resp_mat(k,n,1);
    %         clear temp_y;
    %         for i = 1:length(x_values)
    %             temp_y(i) = resp_mat_360(k,n,i);
    %         end
    %         yi{k,n} = interp1(x_values, temp_y, xi, method);
    %     end
    % end
    % 
    % corr_temp = corrcoef(yi{2,1},yi{2,2});
    % congruency = corr_temp(1,2);
    
    
    %----------------------------------------------------------
    % Bootstrap DI's to classify significantly eye, head, intermed, or unclassif
    bootstraps = 1;
    if bootstraps > 1
    %tic
    for t = 1:bootstraps
        if t/10 == floor(t/10)
            FILE
            t
        end
        resp_mat_boot = [];  % response matrix first
        for k=1:length(unique_condition)
            for n=1:length(unique_stim_type)
                for i=1:length(unique_azimuth)
                    select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
                    if (sum(select) > 0)
                        spike_select = spike_rates(select);
                        for r = 1:repetitions
                            spike_select = spike_select( randperm(length(spike_select)));
                            spike_bootstrap(r) = spike_select(1);   %(sampling with replacement)
                        end
                        resp_mat_boot(k, n, i) = mean(spike_bootstrap);
                    else
                        resp_mat_boot(k, n, i) = 0;                
                    end
                end        
            end
        end
        
        bin = 2; % interval (in degrees) between interpolated points -- 360/bin must be an even integer
        method = 'linear';
        showplots = 0;
        x_values = [0 45 67.5 90 112.5 135 180 225 270 315 360];
        DI_boot{t} = cross_covariance(resp_mat_boot, unique_stim_type, unique_condition, bin, method, showplots, x_values);
    %    DI_avg_boot(t,:) = mean(DI_boot{t});
        DI_avg_boot(t,:) = DI_boot{t}(2,:);
    end
    %toc
    
    % For 95% confidence interval, clip off 2.5% of each side of the distribution
    clip = floor(bootstraps * .025);
    for n = 1:length(unique_stim_type)
    
        DI_avg_sort(:,n) = sort(DI_avg_boot(:,n));
        DI_avg_95(:,n) = DI_avg_sort(clip + 1 : end - clip, n);
        
        % now assign head-centered, eye-centered, or intermediate based on whether confidence interval...
        % includes 0 but not 1 (head):
        if (DI_avg_95(end,n) >= 0 & DI_avg_95(1,n) <= 0) & ~(DI_avg_95(end,n) >= 1.0 & DI_avg_95(1,n) <= 1.0)
            DI_frame(n) = 1;
        % includes 1 but not 0 (eye):
        elseif (DI_avg_95(end,n) >= 1.0 & DI_avg_95(1,n) <= 1.0) & ~(DI_avg_95(end,n) >= 0 & DI_avg_95(1,n) <= 0)
            DI_frame(n) = 2;
        % includes neither ('intermediate', including large/negative DI's):
        elseif ~(DI_avg_95(end,n) >= 1.0 & DI_avg_95(1,n) <= 1.0) & ~(DI_avg_95(end,n) >= 0 & DI_avg_95(1,n) <= 0)
            DI_frame(n) = 3;
        % or includes both (unclassifiable):
        else
            DI_frame(n) = 4;
        end
        
    end
    
    % Plot histograms
    % for n = 1: length(unique_stim_type)
    %     figure(n*11);
    %     hist(DI_avg_sort(:,n));
    % end
    
    % Save bootstrap distributions to .mat file
%     save(['Z:\Users\Chris2\bootstraps\DI\' FILE(1:end-4) '.mat'], 'DI_avg_sort');
    
    else
        DI_frame = [NaN NaN NaN];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMPORARY - uncomment for batch DI_boot (and comment everything below)
    
    % buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', FILE, DI, DI_frame);
    % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\Dir1D_Az_VaryFix_DIBoot.dat'];
    % printflag = 0;
    % if (exist(outfile, 'file') == 0)    %file does not yet exist
    %     printflag = 1;
    % end
    % fid = fopen(outfile, 'a');
    % if (printflag)
    %     fprintf(fid, 'FILE\t DI_minus\t DI_plusminus\t DI_plus\t DI_minus2\t DI_plusminus2\t DI_plus2\t DI_minus3\t DI_plusminus3\t DI_plus3\t Frame_ves\t Frame_vis\t Frame_comb\t');
    %     fprintf(fid, '\r\n');
    % end
    % fprintf(fid, '%s', buff);
    % fprintf(fid, '\r\n');
    % fclose(fid);
    % 
    % toc
    
    % TEMPORARY - uncomment for batch DI_boot (and comment everything below)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %****************************************************************
    % New code for Heading Discrimination Index, HDI
    % (identical to DDI and TDI; see Prince et al. 2002, Nguyenkim and DeAngelis, 2003)
    
    % close;
    
    spike_rates_sqrt = sqrt(spike_rates);
    
    for n=1:length(unique_stim_type)
        for k=1: length(unique_condition)
            for i=1:length(unique_azimuth)
                clear select;
                select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
                SSE(k,n,i) = sum( (spike_rates_sqrt(select) - mean(spike_rates_sqrt(select))).^2 );
            end
            SSE_total(k,n) = sum(SSE(k,n,:));
            max_min_term = sqrt(Max_resp(k,n)) - sqrt(Min_resp(k,n));
            num_conditions = length(unique_azimuth);
            num_trials = repetitions * num_conditions;
            HDI(k,n) = max_min_term / (max_min_term + 2*sqrt(SSE_total(k,n)/(num_trials - num_conditions)));
        end
    end
    
    
    end  % (from 'Proceed with DI, HDI, & Congruency?')
    
    
    % ---------------------------------------------------------------------------------------
    % Rayleigh test (Batschelet, 1981): third test for tuning significance
    % Also see Duffy 1998 (very confusing)
    
    % My take:
    % r = R/n, where R is the true magnitude of the vectorsum (as we compute
    % it), and n is the sample size (which Duffy calls 'total number of spikes', 
    % but should be equivalent to mean spike rate across all 8 directions, since my interval is 1 sec).
    % For grouped data, need to multiply r by c, 
    % where c = (lambda/2) / (sin(lambda/2)) = 1.0262
    % (lambda = 'class length' or group size = 45 deg = pi/4 rads for us)
    
    for k=1:length(unique_condition)
        for n=1:length(unique_stim_type)
            M = squeeze(resp_mat(k,n,:));
            M = M';
            M(3) = [];  % for now, remove the extra sampled points at front of sphere
            M(4) = [];  % (can't do vectorsum and still keep them in there)
            c = 1.0262;              % correction for grouped data at 45 deg intervals
            R = Vec_sum{k,n}(3) * c; % resultant length = vectorsum magnitude * correction
            num = sum(M);            % n = sum of all vectors (same as 'total number of spikes')
            Z(k,n) = R^2 / num;      % Z = n*r^2 = R^2/n (because r = R/n)
        end
    end
    
    %------------------------------------------------------------------
    % Now show vectorsum, HTI, p, spontaneous, and DI at the top of figure
    
    figure(2);
    for n = 1:length(unique_stim_type)
        xoffset = (n-1) * 0.33;
        axes('position',[0.02+xoffset,0.83 0.31,0.12] );
        xlim( [0,100] );
        ylim( [0,length(unique_condition)+0.5] );
        % removed p-perm for now
        text(0.3,length(unique_condition),'Azi          Std         HTI       p-anova        Z        DI');
        for k=1:length(unique_condition)
            K = [2 3 1];  % (for SR_true, which has shift ratios as plus, minus, plusminus)
            if p_anova(k,n) < 0.001
                p_anova(k,n) = 0; % round tiny p_anovas to zero, for display purposes
            end            
            h_text{k,n} = num2str( [Vec_sum{k,n}(1), resp_std(k,n), HTI_temp(k,n), p_anova(k,n), Z(k,n), DI(k,n)], 4);
            text(0.3,length(unique_condition)-k/2, h_text{k,n});
        end
        if n == 2
            text(0.3, length(unique_condition)+0.4, FILE);
        end
        axis off;
    end
    
    % print;
    % close;
     
    proceed1 = 0;
    proceed1 = input('proceed with eye vs head?');
    if proceed1 == 1
    %----------------------------------------------------------
    % Eye-centered model versus head-centered model comparison
    % (fitting all gaze angles simultaneously)
    
    % options = optimset('MaxFunEvals', 10000, 'MaxIter', 5000, 'LargeScale', 'off', 'LevenbergMarquardt', 'on', 'Display', 'off');
    options = optimset('MaxFunEvals', 5000, 'MaxIter', 1000, 'LargeScale', 'off', 'LevenbergMarquardt', 'on', 'Display', 'off');
    A = []; b = []; Aeq = []; beq = []; nonlcon = [];
    
    clear global xdata ydata ydata_merged; % necessary when using tempo_gui multiple times (?)
    global xdata ydata ydata_merged;
    
    xdata = [0 45 67.5 90 112.5 135 180 225 270 315] * pi/180;          % for fitting
    xdata_tran = [-90 -45 0 45 67.5 90 112.5 135 180 225 270] * pi/180; % for display
    for j = 1:repetitions
        xdata_mat(j,:) = xdata;
    end
    xdata = xdata_mat;
    
    for n=1:length(unique_stim_type)
        for k=1: length(unique_condition) 
            for i=1:length(unique_azimuth)
                select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
                ydata{k,n}(:,i) = spike_rates(select)';
            end
        end
        ydata_merged{n} = [ydata{1,n}; ydata{2,n}; ydata{3,n}];
    end
    
    for n=1:length(unique_stim_type)

    % model = ' Charlie Special ';
    % 16 params: 
    % X = [ A1  mu-s  sigma1  K1  K_sig1  DC1
    %       A2  mu    sigma2  K2  K_sig2  DC2
    %       A3  mu+s  sigma3  K3  K_sig3  DC3  ], where s = 20 for eye and 0 for head
    
        LB = [0.001 -2*pi pi/6 0 0.5 0 ; 
              0.001 -2*pi pi/6 0 0.5 0 ;
              0.001 -2*pi pi/6 0 0.5 0];   % lower bounds
    
        UB = [1.5*(Max_resp(1,n)-Min_resp(1,n)) 2*pi 2*pi 0.95 0.95 0.8*Max_resp(1,n) ;
              1.5*(Max_resp(2,n)-Min_resp(2,n)) 2*pi 2*pi 0.95 0.95 0.8*Max_resp(2,n) ;
              1.5*(Max_resp(3,n)-Min_resp(3,n)) 2*pi 2*pi 0.95 0.95 0.8*Max_resp(3,n)];   % upper bounds
            
        % initial parameter guesses
        X0 = [];
        for k = 1:length(unique_condition)
            
            % A = peak-trough modulation
            X0(k,1) = Max_resp(k,n) - Min_resp(k,n);
            if X0(k,1) < 1
                X0(k,1) = 1;
            end
            
            % mu = azimuth of max response
            whereis_max = logical(resp_mat == Max_resp(k,n));
            whereis_max = squeeze(whereis_max(k,n,:));
            max_azi{k,n} = unique_azimuth(find(whereis_max));
            if length(max_azi{k,n}) > 1   % need these to be unique, and believe it or not they sometimes aren't
                max_azi{k,n} = max_azi{k,n}(1);
            end
            X0(k,2) = max_azi{k,n} * pi/180;
    
            % DC = min response
            X0(k,6) = Min_resp(k,n);
            
            % search for best starting values of sigma, K, and K-sig
            N = 40;
            min_err = 9999999999999999.99;
    		x3range = LB(k,3) : (UB(k,3)-LB(k,3))/N : UB(k,3);
            x4range = LB(k,4) : (UB(k,4)-LB(k,4))/N : UB(k,4); 
            x5range = LB(k,5) : (UB(k,5)-LB(k,5))/N : UB(k,5);
    		for i = 1:N
                x3temp = x3range(i);
                if x3temp == 0
                    x3temp = 0.0001;
                end
                for j = 1:N
                    x4temp = x4range(j);
                    for h = 1:N
                        x5temp = x5range(h);
                        if x5temp == 0
                            x5temp = 0.0001;
                        end
                        x_temp = [X0(k,1) X0(k,2) x3temp x4temp x5temp X0(k,6)];
                        error = VF_1D_Curvefit_err_temp(x_temp,n,k,repetitions);
                        if (error < min_err)
                            x3min = x3temp; x4min = x4temp; x5min = x5temp;
                            min_err = error;    
                        end
                    end
                end
    		end
    		X0(k,3) = x3min; X0(k,4) = x4min; X0(k,5) = x5min;
    
        end
        
        % fit multiple times with some jitter in the initial params
        N_reps = 15;
        wiggle = 0.3;
        clear testpars_head testpars_eye;
        global stimtype;
        stimtype = n;
        min_err_eye = 9999999999999999.99; min_err_head = 9999999999999999.99;
        for j=1:N_reps
    
            FILE
            n        
            j
            
            rand_factor = rand(size(X0)) * wiggle + (1-wiggle/2); % ranges from 1-wiggle/2 -> 1 + wiggle/2
            temp_X0 = X0 .* rand_factor;
    
            testpars_head = fmincon('VF_1D_Curvefit_err_head', temp_X0, A, b, Aeq, beq, LB, UB, nonlcon, options);
            err_head = VF_1D_Curvefit_err_head(testpars_head);
            if (err_head < min_err_head)
                testpars_head_min = testpars_head;
                min_err_head = err_head;    
            end
    
            testpars_eye = fmincon('VF_1D_Curvefit_err_eye', temp_X0, A, b, Aeq, beq, LB, UB, nonlcon, options);
            err_eye = VF_1D_Curvefit_err_eye(testpars_eye);
            if (err_eye < min_err_eye)
                testpars_eye_min = testpars_eye;
                min_err_eye = err_eye;    
            end
    
        end
    
        X_0_head{n} = X0;
        X_0_eye{n} = X0;
        X_head{n} = testpars_head_min;
        X_eye{n} = testpars_eye_min;
        yfit_head{n} = (VF_1D_Curvefit_head(X_head{n},xdata));
        yfit_eye{n} = (VF_1D_Curvefit_eye(X_eye{n},xdata));
    
        % R^2's using means, but need to parse and avg the data first
        for k = 1:length(unique_condition)
            for i = 1:length(unique_azimuth)
                ydata_stacked{n}(i+length(unique_azimuth)*(k-1)) = mean(ydata_merged{n}((k-1)*repetitions+1:k*repetitions,i));
                yfit_head_stacked{n}(i+length(unique_azimuth)*(k-1)) = mean(yfit_head{n}((k-1)*repetitions+1:k*repetitions,i));
                yfit_eye_stacked{n}(i+length(unique_azimuth)*(k-1)) = mean(yfit_eye{n}((k-1)*repetitions+1:k*repetitions,i));
            end
        end
    
        clear coef P;
        [coef,P] = corrcoef(ydata_stacked{n},yfit_head_stacked{n});
        R_head(n) = coef(1,2);
        rsquared_head(n) = coef(1,2)^2;
        p_fit_head(n) = P(1,2);
    
        clear coef P;
        [coef,P] = corrcoef(ydata_stacked{n},yfit_eye_stacked{n});
        R_eye(n) = coef(1,2);
        rsquared_eye(n) = coef(1,2)^2;
        p_fit_eye(n) = P(1,2);
        
    %    for partial correlation, need corrcoef between eye and head themselves
        clear coef P;
        [coef,P] = corrcoef(yfit_head_stacked{n},yfit_eye_stacked{n});
        R_headeye(n) = coef(1,2);
        rsquared_headeye(n) = coef(1,2)^2;
        partialcorr_head(n) = (R_head(n) - R_eye(n) * R_headeye(n)) / sqrt( (1-rsquared_eye(n)) * (1-rsquared_headeye(n)) );
        partialcorr_eye(n) = (R_eye(n) - R_head(n) * R_headeye(n)) / sqrt( (1-rsquared_head(n)) * (1-rsquared_headeye(n)) );
        partialZ_head(n) = 0.5 * log((1+partialcorr_head(n))/(1-partialcorr_head(n))) / (1/sqrt(length(unique_azimuth)*length(unique_condition)-3));
        partialZ_eye(n) = 0.5 * log((1+partialcorr_eye(n))/(1-partialcorr_eye(n))) / (1/sqrt(length(unique_azimuth)*length(unique_condition)-3));
        
    end
    
    
    figure(5);
    orient landscape;
    set(5,'Position', [50,25 1200,900], 'Name', '1D Azimuth VaryFix CurveFits');
    axis off;
    
    X_head_temp = X_head;
    X_eye_temp = X_eye;
    for n=1:length(unique_stim_type)
        
        X_head_temp{n}(1,2) = X_head_temp{n}(2,2);
        X_head_temp{n}(3,2) = X_head_temp{n}(2,2);
        X_eye_temp{n}(1,2) = X_eye_temp{n}(2,2) + pi/8;
        X_eye_temp{n}(3,2) = X_eye_temp{n}(2,2) - pi/8;
        
        for k=1: length(unique_condition) 
            figure(5);
            xoffset = (n-1) * 0.33;
            yoffset = (k-1) * -0.28;
            axes('position',[0.02+xoffset 0.62+yoffset 0.28 0.20]);
            errorbar(xdata_tran, resp_mat_tran(k,n,:), resp_mat_std_tran(k, n, :)/sqrt(repetitions), 'ro');
            hold on;
            spon_line2(1:length(xdata_tran)) = spon_resp(k);  % display spontaneous as a dotted line
            plot(xdata_tran,spon_line2,'b:');
            
            x_smooth = 0:0.01:2*pi;
            y_smooth_head = VF_1D_Curvefit(X_head_temp{n}(k,:),x_smooth);
            y_smooth_eye = VF_1D_Curvefit(X_eye_temp{n}(k,:),x_smooth);
            x_smooth_tran = xdata_tran(1):0.01:xdata_tran(end);
            y_smooth_tran_head = VF_1D_Curvefit(X_head_temp{n}(k,:),x_smooth_tran);
            y_smooth_tran_eye = VF_1D_Curvefit(X_eye_temp{n}(k,:),x_smooth_tran);
    
            plot(x_smooth_tran,y_smooth_tran_head,'b',x_smooth_tran,y_smooth_tran_eye,'g');
            xlim( [xdata_tran(1) xdata_tran(end)] );
            set(gca, 'xtick', xdata_tran);
            set(gca, 'xdir' , 'reverse');
            set(gca, 'xticklabel','-90|-45|0|45| |90| |135|180|225|270');
            if k == length(unique_condition)
                xlabel('Azimuth');
            end
            
            % show peaks (mu) for eye and head
            param_text = [num2str(X_head_temp{n}(k,2)*180/pi,4) '       ' num2str(X_eye_temp{n}(k,2)*180/pi,4)];
            clear y_lim y_range;
            y_lim = ylim;
            y_range = y_lim(2)-y_lim(1);
            text(4, y_lim(2)+.12*y_range, 'azi-head    azi-eye');
            text(4, y_lim(2)+.04*y_range, param_text);
        end
        
    
        
        % and the R^2's at the top
        xoffset = (n-1) * 0.33;
        axes('position',[0.02+xoffset,0.83 0.31,0.12] );
        xlim( [0,100] );
        ylim( [0,length(unique_condition)] );
        text(0.3,length(unique_condition),'     Rsquared-head          Rsquared-eye     ');
        text(0.3,length(unique_condition)-0.5, ['     ' num2str(rsquared_head(n)) '                    ' num2str(rsquared_eye(n))]);
        if n == 2
            text(0.35, length(unique_condition)+0.4, FILE);
        end
        axis off;
        
    end
    
    % print;
    % close;
    
    
    % % %--------------------------------------------------------------------------------------
    % % % Bootstrap eye vs. head model fits  --- OMIT --- need to re-fit each bootstrap; will do later *****
    % % bootstraps = 1;
    % % if bootstraps > 1
    % %     
    % % for t = 1:bootstraps
    % %     resp_mat_boot = [];  % response matrix first
    % %     for k=1:length(unique_condition)
    % %         for n=1:length(unique_stim_type)
    % %             for i=1:length(unique_azimuth)
    % %                 select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
    % %                 if (sum(select) > 0)
    % %                     spike_select = spike_rates(select);
    % %                     for r = 1:repetitions
    % %                         spike_select = spike_select( randperm(length(spike_select)));
    % %                         spike_bootstrap(r) = spike_select(1);   %(sampling with replacement)
    % %                     end
    % %                     resp_mat_boot(k, n, i) = mean(spike_bootstrap);
    % %                 else
    % %                     resp_mat_boot(k, n, i) = 0;                
    % %                 end
    % %             end        
    % %         end
    % %     end
    % % 
    % %     % then re-compute R-squared for each model using previous params
    % %     clear ydata_stacked_boot;
    % %     for n=1:length(unique_stim_type)
    % %         for k=1: length(unique_condition) 
    % %             for i=1:length(unique_azimuth)
    % %                 ydata_stacked_boot{n}(i+length(unique_azimuth)*(k-1)) = resp_mat_boot(k,n,i);
    % %             end
    % %         end
    % %         clear coef P;
    % %         [coef,P] = corrcoef(ydata_stacked_boot{n},yfit_head_stacked{n});
    % %         rsquared_head_boot(n,t) = coef(1,2)^2;
    % %         clear coef P;
    % %         [coef,P] = corrcoef(ydata_stacked_boot{n},yfit_eye_stacked{n});
    % %         rsquared_eye_boot(n,t) = coef(1,2)^2;
    % %     end
    % %     
    % % end
    % %     
    % % % For 95% confidence interval, clip off 2.5% of each side of the distribution
    % % clip = floor(bootstraps * .025);
    % % for n = 1: length(unique_stim_type)
    % %     rsquared_head_sort(n,:) = sort(rsquared_head_boot(n,:));
    % %     rsquared_head_95(n,:) = rsquared_head_sort(n, clip + 1 : end - clip);
    % %     rsquared_eye_sort(n,:) = sort(rsquared_eye_boot(n,:));
    % %     rsquared_eye_95(n,:) = rsquared_eye_sort(n, clip + 1 : end - clip);
    % %     
    % %     % then assign significance based on whether confidence intervals overlap
    % %     if rsquared_head(n) > rsquared_eye(n)
    % %         if rsquared_head_95(n,1) > rsquared_eye_95(n,end)
    % %             models_sig_diff(n) = 1;
    % %         else
    % %             models_sig_diff(n) = 0;
    % %         end
    % %     else
    % %         if rsquared_eye_95(n,1) > rsquared_head_95(n,end)
    % %             models_sig_diff(n) = 1;
    % %         else
    % %             models_sig_diff(n) = 0;
    % %         end
    % %     end
    % % 
    % % end
    % % 
    % % % % Plot histograms
    % % % for n = 1: length(unique_stim_type)
    % % %     figure(n*10);
    % % %     hist(rsquared_head_boot(n,:));
    % % %     title('Rsquared - Head');
    % % %     figure(n*10+1);
    % % %     hist(rsquared_eye_boot(n,:));
    % % %     title('Rsquared - Eye');
    % % % end
    % % 
    % % % % Save bootstrap distributions to .mat file
    % % % if bootstraps > 1
    % % %     save(['Z:\Users\Chris2\bootstraps\' FILE '.mat'], 'rsquared_eye_sort' 'rsquared_head_sort');
    % % % end
    % % 
    % % else
    % %     models_sig_diff = [NaN NaN NaN];
    % % end
    % % 
    % % %--------------------------------------------------------------------------------------
     
    else  % to skip entire section above, following PROCEED1 input
        R_head = [0 0 0];
        R_eye = [0 0 0];
        rsquared_head = [0 0 0];
        rsquared_eye = [0 0 0];
    end   % END Proceed1 loop
    
    
    
    % -------------------------------------------------------------------------
    % Curve fitting
    % -------------------------------------------------------------------------
    %
    %
    
    proceed2 = 1;
    % proceed2 = input('proceed with fitting?');
    if proceed2 == 1
    
    options = optimset('MaxFunEvals', 10000, 'MaxIter', 5000, 'LargeScale', 'off', 'LevenbergMarquardt', 'on', 'Display', 'off');
    A = []; b = []; Aeq = []; beq = []; nonlcon = [];
    
    clear global xdata; % necessary when using tempo_gui multiple times (?)
    global xdata;
    
    xdata = [0 45 67.5 90 112.5 135 180 225 270 315] * pi/180;          % for fitting
    xdata_tran = [-90 -45 0 45 67.5 90 112.5 135 180 225 270] * pi/180; % for display
    for j = 1:repetitions
        xdata_mat(j,:) = xdata;
    end
    xdata = xdata_mat;
    
    figure(3);
    orient landscape;
    set(3,'Position', [50,25 1200,900], 'Name', '1D Azimuth VaryFix CurveFits');
    axis off;
    for n=1:length(unique_stim_type)
        for k=1: length(unique_condition) 
            figure(3);
            xoffset = (n-1) * 0.33;
            yoffset = (k-1) * -0.28;
            axes('position',[0.02+xoffset 0.62+yoffset 0.28 0.20]);
            errorbar(xdata_tran, resp_mat_tran(k,n,:), resp_mat_std_tran(k, n, :)/sqrt(repetitions), 'ro');
            hold on;
            spon_line2(1:length(xdata_tran)) = spon_resp(k);  % display spontaneous as a dotted line
            plot(xdata_tran,spon_line2,'b:');
    
          % fitting individual trial data, so:
            clear global ydata
            global ydata
            for i=1:length(unique_azimuth)
                i
                select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
                ydata(:,i) = spike_rates(select)';
            end
    
    %**********************************************
    % Function 4 (6 params, see VF_1D_Curvefit.m):
    
            model = ' Charlie Special ';
            param_label = 'A         mu     sigma     K     K-sig      DC';
            lb = [0.001 -2*pi pi/6 0 0.5 0];   % lower bounds
            ub = [1.5*(Max_resp(k,n)-Min_resp(k,n)) 2*pi 2*pi 0.95 0.95 0.8*Max_resp(k,n)];   % upper bounds
            
          % initial parameter guesses
            x0 = [];
            
            x0(1) = Max_resp(k,n) - Min_resp(k,n);   % A = peak-trough modulation
            if x0(1) < 1
                x0(1) = 1;
            end
    
    %     % mu = azimuth of vector sum:
    %         x0(2) = Vec_sum{k,n}(1) * pi/180;
    %         if x0(2) < 0
    %             x0(2) = x0(2) + 2*pi;
    %         end
    
          % OR, mu = azimuth of max response:
            whereis_max = logical(resp_mat == Max_resp(k,n));
            whereis_max = squeeze(whereis_max(k,n,:));
            max_azi{k,n} = unique_azimuth(find(whereis_max));
            if length(max_azi{k,n}) > 1   % need these to be unique, and believe it or not they sometimes aren't
                max_azi{k,n} = max_azi{k,n}(1);
            end
            x0(2) = max_azi{k,n} * pi/180;
    
            x0(6) = Min_resp(k,n);   % DC = min response
            
          % search for best starting values of x0(3), x0(4), and x0(5)
            N = 30;
            min_err = 9999999999999999.99;
    		x3range = lb(3) : (ub(3)-lb(3))/N : ub(3);
            x4range = lb(4) : (ub(4)-lb(4))/N : ub(4); 
            x5range = lb(5) : (ub(5)-lb(5))/N : ub(5); 
    		for i = 1:N
                x3temp = x3range(i);
                if x3temp == 0
                    x3temp = 0.0001;
                end
                for j = 1:N
                    x4temp = x4range(j);
                    for h = 1:N
                        x5temp = x5range(h);
                        if x5temp == 0
                            x5temp = 0.0001;
                        end
                        x_temp = [x0(1) x0(2) x3temp x4temp x5temp x0(6)];
                        error = VF_1D_Curvefit_err(x_temp);
                        if (error < min_err)
                            x3min = x3temp; x4min = x4temp; x5min = x5temp;
                            min_err = error;    
                        end
                    end
                end
    		end
    		x0(3) = x3min; x0(4) = x4min; x0(5) = x5min;
    
    %**********************************************
    
            % fit multiple times with some jitter in the initial params
            N_reps = 30;
            wiggle = 0.3;
            clear testpars fval_temp exitflag_temp;
            for j=1:N_reps
                rand_factor = rand(length(x0),1) * wiggle + (1-wiggle/2); % ranges from 1-wiggle/2 -> 1 + wiggle/2
                temp_x0 = x0' .* rand_factor;
                [testpars{j}, fval_temp(j), exitflag_temp(j)] = fmincon('VF_1D_Curvefit_err', temp_x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    %            testpars{j} = fmincon('VF_1D_Curvefit_err', temp_x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
                err{k,n}(j) = VF_1D_Curvefit_err(testpars{j});
                yfit_test{j} = (VF_1D_Curvefit(testpars{j},xdata));
                [coef_test,P_test] = corrcoef(ydata',yfit_test{j}');   % R^2's using means
                rsquared_test(j) = coef_test(1,2)^2;
            end
    
          % so the best-fit values of x are:
            best_err = find(err{k,n} == min(err{k,n}));
            best_r2 = find(rsquared_test(best_err) == max(rsquared_test(best_err)));  % if multiple occurences of lowest 
            x{k,n} = testpars{best_err(best_r2)};                                     % error, use the one with best R^2
          % and the rest...
            final_error(k,n) = min(err{k,n});
            fval(k,n) = fval_temp(best_err(best_r2));
            exitflag(k,n) = exitflag_temp(best_err(best_r2));
            yfit{k,n} = (VF_1D_Curvefit(x{k,n},xdata));
            [coef,P] = corrcoef(mean(ydata)',mean(yfit{k,n})');   % R^2's using means
            rsquared(k,n) = coef(1,2)^2;
            p_fit(k,n) = P(1,2);
            x_0{k,n} = x0;
    
          % finish plotting
            x_smooth = 0:0.01:2*pi;
            y_smooth = (VF_1D_Curvefit(x{k,n},x_smooth));
            x_smooth_tran = xdata_tran(1):0.01:xdata_tran(end);
            y_smooth_tran = (VF_1D_Curvefit(x{k,n},x_smooth_tran));
            plot(x_smooth_tran,y_smooth_tran,'b');
            xlim( [xdata_tran(1) xdata_tran(end)] );
            set(gca, 'xtick', xdata_tran);
            set(gca, 'xdir' , 'reverse');
            set(gca, 'xticklabel','-90|-45|0|45| |90| |135|180|225|270');
            if k == length(unique_condition)
                xlabel('Azimuth');
            end
    
            % show params for each fit
            param_text = [num2str(x{k,n}(1),4) '    ' num2str(x{k,n}(2)*180/pi,4) '    ' num2str(x{k,n}(3)*180/pi,4) '    ' num2str(x{k,n}(4),4) '    ' num2str(x{k,n}(5),4) '    ' num2str(x{k,n}(6),4)];            
            x0_text = [num2str(x0(1),4) '    ' num2str(x0(2)*180/pi,4) '    ' num2str(x0(3)*180/pi,4) '    ' num2str(x0(4),4) '    ' num2str(x0(5),4) '    ' num2str(x0(6),4)];
            y_lim = ylim;
            y_range = y_lim(2)-y_lim(1);
            text(4, y_lim(2)+.20*y_range, param_label);
            text(4, y_lim(2)+.12*y_range, x0_text);
            text(4, y_lim(2)+.04*y_range, param_text);
    
          % take peak of fitted function as preferred heading 
            peak(k,n) = x_smooth(find(y_smooth == max(y_smooth))) * 180/pi;
    %*******peak(k,n) = x{k,n}(2); (or should it be phi/mu parameter? Why aren't they identical? *******
    
        end
    end
    
    
    % plot error for each jittered fit, to make sure fits are converging properly
    figure(4);
    orient landscape;
    set(4,'Position', [95,25 1200,900], 'Name', '1D Az VF CurveFit Error');
    axis off;
    for n=1:length(unique_stim_type)
        for k=1: length(unique_condition) 
            xoffset = (n-1) * 0.33;
            yoffset = (k-1) * -0.28;
            axes('position',[0.02+xoffset 0.62+yoffset 0.28 0.22]);
            plot(err{k,n});
            x_text = ['min err = ' num2str(min(err{k,n})) '   fval = ' num2str(fval(k,n))];
            xlabel(x_text);
            if n == 2 & k == 1
                title([FILE '   ' model]);
            end
        end
    end
    
    % print;
    % close;
    
    
    % -------------------------------------------------------------------------
    % Compute shift ratios from the curve-fitted peaks
    for n = 1:length(unique_stim_type)
        
        if abs(peak(1,n) - peak(2,n)) > 180   % (to deal with wrap-around problem)
            if peak(1,n) > peak(2,n)
                SR_fit(1,n) = (peak(1,n) - (peak(2,n)+360)) / unique_condition(3);
            else
                SR_fit(1,n) = ((peak(1,n)+360) - peak(2,n)) / unique_condition(3);
            end
        else
            SR_fit(1,n) = (peak(1,n) - peak(2,n)) / unique_condition(3);      % minus
        end
        
        if abs(peak(1,n) - peak(3,n)) > 180
            if peak(1,n) > peak(3,n)
                SR_fit(2,n) = (peak(1,n) - (peak(3,n)+360)) / (2*unique_condition(3));
            else
                SR_fit(2,n) = ((peak(1,n)+360) - peak(3,n)) / (2*unique_condition(3));
            end
        else
            SR_fit(2,n) = (peak(1,n) - peak(3,n)) / (2*unique_condition(3));  % plusminus
        end
        
        if abs(peak(2,n) - peak(3,n)) > 180
            if peak(2,n) > peak(3,n)
                SR_fit(3,n) = (peak(2,n) - (peak(3,n)+360)) / unique_condition(3);      
            else
                SR_fit(3,n) = ((peak(2,n)+360) - peak(3,n)) / unique_condition(3);
            end
        else
            SR_fit(3,n) = (peak(2,n) - peak(3,n)) / unique_condition(3);      % plus            
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Lastly, display fit params, etc. at the top of the figure
    
    figure(3);
    for n = 1:length(unique_stim_type)
        xoffset = (n-1) * 0.33;
        axes('position',[0.02+xoffset,0.83 0.31,0.12] );
        xlim( [0,100] );
        ylim( [0,length(unique_condition)] );
        text(0.3,length(unique_condition),'Azi-fit    exitflag   error        r^2          p-fit      SR-fit');
        for k=1:length(unique_condition)
            if p_fit(k,n) < 0.001
                p_fit(k,n) = 0; % round tiny p_fits to zero, for display purposes
            end            
            h_text{k,n} = num2str( [peak(k,n), exitflag(k,n), final_error(k,n), rsquared(k,n), p_fit(k,n), SR_fit(k,n)], 4);
            text(0.3,length(unique_condition)-k/2, h_text{k,n});
        end
        if n == 2
            text(0.35, length(unique_condition)+0.4, FILE);
        end
        if n == 3
            text(0.35, length(unique_condition)+0.4, model);
        end
        axis off;
    end
    
    % print;
    % close;
    
    
    %-------------------------------------------------------------------------
    % Fisher information analysis: compute estimate of variance at each point on 
    % fitted function, using the average fano factor for this cell
    %-------------------------------------------------------------------------
    
    for n=1:length(unique_stim_type)
        for k=1: length(unique_condition)
            
            fano(k,n) = mean(squeeze(resp_mat_std(k,n,:).^2) ./ squeeze(resp_mat(k,n,:)) );
    
            x_smooth_tran = xdata_tran(1):0.01:xdata_tran(end);
    		y_smooth_tran = VF_1D_Curvefit(x{k,n},x_smooth_tran);
            fisher_variance = y_smooth_tran * fano(k,n);
            fisher_slope = diff(y_smooth_tran);
            fisher{k,n} = (fisher_slope.^2) ./ fisher_variance(1:end-1);
            
            %%%%% TEMP - for testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         figure; errorbar(x_smooth_tran,y_smooth_tran,fisher_variance);
    %         xlim( [xdata_tran(1) xdata_tran(end)] );
    % 		set(gca, 'xtick', xdata_tran);
    % 		set(gca, 'xdir' , 'reverse');
    %         set(gca, 'xticklabel','-90|-45|0|45| |90| |135|180|225|270');
    %         legend('variance');
            
            figure; plot(x_smooth_tran,y_smooth_tran/max(y_smooth_tran));
            hold on;
            plot(x_smooth_tran(1:end-1), -(fisher_slope/max(fisher_slope)), 'g'); % sign-reversed for plotting; actual sign doesn't matter
            plot(x_smooth_tran(1:end-1), fisher{k,n}/max(fisher{k,n}), 'r');
            xlim( [xdata_tran(1) xdata_tran(end)] );
    		set(gca, 'xtick', xdata_tran);
    		set(gca, 'xdir' , 'reverse');
            set(gca, 'xticklabel','-90|-45|0|45| |90| |135|180|225|270');
            legend('normalized response','normalized slope','normalized fisher info');
            %%%%% TEMP - for testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
    end
    
    
    else   % to skip entire curve-fitting section, following PROCEED2 input
        peak = zeros(3,3);
        rsquared = zeros(3,3);
        p_fit = zeros(3,3);
        SR_fit = zeros(3,3);
    end 
    
    
    % %--------------------------------------------------------------------------------------
    % % Bootstrap vectorsum and shift ratios   [OBSOLETE]
    % bootstraps = 1;
    % if bootstraps > 1
    %     
    % for t = 1:bootstraps
    %     resp_mat_boot = [];  % response matrix first
    %     for k=1:length(unique_condition)
    %         for n=1:length(unique_stim_type)
    %             for i=1:length(unique_azimuth)
    %                 select = logical( (azimuth==unique_azimuth(i)) & (condition==unique_condition(k)) & (stim_type==unique_stim_type(n)) );
    %                 if (sum(select) > 0)
    %                     spike_select = spike_rates(select);
    %                     for r = 1:repetitions
    %                         spike_select = spike_select( randperm(length(spike_select)));
    %                         spike_bootstrap(r) = spike_select(1);   %(sampling with replacement)
    %                     end
    %                     resp_mat_boot(k, n, i) = mean(spike_bootstrap);
    %                 else
    %                     resp_mat_boot(k, n, i) = 0;                
    %                 end
    %             end        
    %         end
    %     end
    % 
    %     % then bootstrap vectorsum for each gaze condition...
    %     for n = 1: length(unique_stim_type)
    %         for k = 1: length(unique_condition)
    %             N = squeeze(resp_mat_boot(k,n,:));
    %             N = N';
    %             N(3) = [];  % for now, remove the extra sampled points at front of sphere
    %             N(4) = [];  % (can't do vectorsum and still keep them in there)
    %             [Azi_boot, Ele_boot, Amp_boot] = vectorsum(N);
    %             Vec_sum_boot{k,n} = [Azi_boot, Ele_boot, Amp_boot];
    %         end
    %     end
    % 
    %     % which is used to compute SR_true_boots
    %     for n = 1: length(unique_stim_type)
    %         if abs(Vec_sum_boot{3*n-1}(1) - Vec_sum_boot{3*n}(1)) > 180
    %             SR_true_boot{n}(1) = ((Vec_sum_boot{3*n-1}(1)+360) - Vec_sum_boot{3*n}(1)) / unique_condition(3);      
    %         else
    %             SR_true_boot{n}(1) = (Vec_sum_boot{3*n-1}(1) - Vec_sum_boot{3*n}(1)) / unique_condition(3);     % plus            
    %         end
    %         if abs(Vec_sum_boot{3*n-2}(1) - Vec_sum_boot{3*n-1}(1)) > 180
    %             SR_true_boot{n}(2) = ((Vec_sum_boot{3*n-2}(1)+360) - Vec_sum_boot{3*n-1}(1)) / unique_condition(3);
    %         else
    %             SR_true_boot{n}(2) = (Vec_sum_boot{3*n-2}(1) - Vec_sum_boot{3*n-1}(1)) / unique_condition(3);   % minus
    %         end
    %         if abs(Vec_sum_boot{3*n-2}(1) - Vec_sum_boot{3*n}(1)) > 180
    %             SR_true_boot{n}(3) = ((Vec_sum_boot{3*n-2}(1)+360) - Vec_sum_boot{3*n}(1)) / (2*unique_condition(3));
    %         else
    %             SR_true_boot{n}(3) = (Vec_sum_boot{3*n-2}(1) - Vec_sum_boot{3*n}(1)) / (2*unique_condition(3)); % plusminus
    %         end
    % 
    % % formerly used rotation.m, but not necessary for horiz gaze shift
    % %         [E_plus, E_minus, SR_plus, SR_minus, SR_plusminus, NR_plus, NR_minus, NR_plusminus] = rotation_no_plots(Vec_sum_boot{3*n-1}, Vec_sum_boot{3*n}, Vec_sum_boot{3*n-2}, unique_condition(3), gaze_dir);
    % %         SR_plus_boot(n,t) = SR_plus;
    % %         SR_minus_boot(n,t) = SR_minus;
    % %         SR_plusminus_boot(n,t) = SR_plusminus;
    % %         NR_plus_boot(n,t) = NR_plus;
    % %         NR_minus_boot(n,t) = NR_minus;
    % %         NR_plusminus_boot(n,t) = NR_plusminus;
    %     end
    %     
    % end
    % 
    % % For 95% confidence interval, clip off 2.5% of each side of the distribution
    % clip = floor(bootstraps * .025);
    % for n = 1: length(unique_stim_type)
    %     SR_plus_sort(n,:) = sort(SR_plus_boot(n,:));
    %     SR_plus_95(n,:) = SR_plus_sort(n, clip + 1 : end - clip);
    %     SR_minus_sort(n,:) = sort(SR_minus_boot(n,:));
    %     SR_minus_95(n,:) = SR_minus_sort(n, clip + 1 : end - clip);
    %     SR_plusminus_sort(n,:) = sort(SR_plusminus_boot(n,:));
    %     SR_plusminus_95(n,:) = SR_plusminus_sort(n, clip + 1 : end - clip);
    %     
    %     % now assign head-centered, eye-centered, or intermediate based on whether confidence interval...
    %     % includes 1 but not 0 (eye):
    %     if (SR_plus_95(n,end) >= 1.0 & SR_plus_95(n,1) <= 1.0) & ~(SR_plus_95(n,end) >= 0 & SR_plus_95(n,1) <= 0)
    %         plus_frame(n) = 1;
    %     % includes 0 but not 1 (head):
    %     elseif (SR_plus_95(n,end) >= 0 & SR_plus_95(n,1) <= 0) & ~(SR_plus_95(n,end) >= 1.0 & SR_plus_95(n,1) <= 1.0)
    %         plus_frame(n) = 2;
    %     % includes neither ('intermediate', including large/negative SRs):
    %     elseif ~(SR_plus_95(n,end) >= 1.0 & SR_plus_95(n,1) <= 1.0) & ~(SR_plus_95(n,end) >= 0 & SR_plus_95(n,1) <= 0)
    %         plus_frame(n) = 3;
    %     % or includes both (unclassifiable):
    %     else
    %         plus_frame(n) = 4;
    %     end
    %     
    %     if (SR_minus_95(n,end) >= 1.0 & SR_minus_95(n,1) <= 1.0) & ~(SR_minus_95(n,end) >= 0 & SR_minus_95(n,1) <= 0)
    %         minus_frame(n) = 1;
    %     elseif (SR_minus_95(n,end) >= 0 & SR_minus_95(n,1) <= 0) & ~(SR_minus_95(n,end) >= 1.0 & SR_minus_95(n,1) <= 1.0)
    %         minus_frame(n) = 2;
    %     elseif ~(SR_minus_95(n,end) >= 1.0 & SR_minus_95(n,1) <= 1.0) & ~(SR_minus_95(n,end) >= 0 & SR_minus_95(n,1) <= 0)
    %         minus_frame(n) = 3;
    %     else
    %         minus_frame(n) = 4;
    %     end
    % 
    %     if (SR_plusminus_95(n,end) >= 1.0 & SR_plusminus_95(n,1) <= 1.0) & ~(SR_plusminus_95(n,end) >= 0 & SR_plusminus_95(n,1) <= 0)
    %         plusminus_frame(n) = 1;
    %     elseif (SR_plusminus_95(n,end) >= 0 & SR_plusminus_95(n,1) <= 0) & ~(SR_plusminus_95(n,end) >= 1.0 & SR_plusminus_95(n,1) <= 1.0)
    %         plusminus_frame(n) = 2;
    %     elseif ~(SR_plusminus_95(n,end) >= 1.0 & SR_plusminus_95(n,1) <= 1.0) & ~(SR_plusminus_95(n,end) >= 0 & SR_plusminus_95(n,1) <= 0)
    %         plusminus_frame(n) = 3;
    %     else
    %         plusminus_frame(n) = 4;
    %     end
    %     
    %     boot_data{n} = [plus_frame(n) minus_frame(n) plusminus_frame(n)];
    %     
    % end
    % 
    % % Plot histograms
    % % for n = 1: length(unique_stim_type)
    % %     figure(n*10);
    % %     hist(SR_plus_boot(n,:));
    % %     title('Shift Ratio Plus');
    % %     figure(n*10+1);
    % %     hist(SR_minus_boot(n,:));
    % %     title('Shift Ratio Minus');
    % %     figure(n*10+2);
    % %     hist(SR_plusminus_boot(n,:));
    % %     title('Shift Ratio Plusminus');
    % % end
    % 
    % % Save bootstrap distributions to .mat file
    % save(['Z:\Users\Chris2\bootstraps\' FILE '.mat'], 'SR_plus_boot', 'SR_minus_boot', 'SR_plusminus_boot');
    % 
    % else
    %     boot_data{1} = [999 999 999];
    %     boot_data{2} = [999 999 999];
    %     boot_data{3} = [999 999 999];
    % end
    
    
    % if (proceed1 == 1 | proceed2 == 1)
    % %--------------------------------------------------------------------------------------
    % % Write out all data to a cumulative summary file
    % 
    % % if proceed2 == 1    ***** OMIT *****   (keeps giving 'file may be corrupt' error msg)
    % %     % save fit params to .mat file
    % %     save(['Z:\Users\Chris2\fitparams\' 'fit_' FILE '.mat'], 'xdata', 'x', 'x_0', 'peak', 'rsquared', 'p_fit', 'SR_fit');
    % % end
    % 
    % buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %3.0f\t %s\t', ...
    %        FILE, spon_resp, Min_resp(:), Max_resp(:), Vec_sum{:}, HTI_temp(:), Z(:), p_anova(:), resp_std(:), gf_anovas{:}, SR_true{:}, boot_data{:}, peak(:), rsquared(:), p_fit(:), SR_fit(:), EndTrial-(BegTrial-1), dir_text);
    %        % currently 153 fields (11/05)
    % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\Dir1D_Az_VaryFix_Sum.dat'];
    % printflag = 0;
    % if (exist(outfile, 'file') == 0)    %file does not yet exist
    %     printflag = 1;
    % end
    % fid = fopen(outfile, 'a');
    % if (printflag)
    %     fprintf(fid, 'FILE\t Spon_minus\t Spon_zero\t Spon_plus\t Min_minus\t Min_zero\t Min_plus\t Min_minus2\t Min_zero2\t Min_plus2\t Min_minus3\t Min_zero3\t Min_plus3\t Max_minus\t Max_zero\t Max_plus\t Max_minus2\t Max_zero2\t Max_plus2\t Max_minus3\t Max_zero3\t Max_plus3\t Azi_minus\t Ele_minus\t Amp_minus\t Azi_zero\t Ele_zero\t Amp_zero\t Azi_plus\t Ele_plus\t Amp_plus\t Azi_minus2\t Ele_minus2\t Amp_minus2\t Azi_zero2\t Ele_zero2\t Amp_zero2\t Azi_plus2\t Ele_plus2\t Amp_plus2\t Azi_minus3\t Ele_minus3\t Amp_minus3\t Azi_zero3\t Ele_zero3\t Amp_zero3\t Azi_plus3\t Ele_plus3\t Amp_plus3\t HTI_minus\t HTI_zero\t HTI_plus\t HTI_minus2\t HTI_zero2\t HTI_plus2\t HTI_minus3\t HTI_zero3\t HTI_plus3\t Z_minus\t Z_zero\t Z_plus\t Z_minus2\t Z_zero2\t Z_plus2\t Z_minus3\t Z_zero3\t Z_plus3\t P_anova_minus\t P_anova_zero\t P_anova_plus\t P_anova_minus2\t P_anova_zero2\t P_anova_plus2\t P_anova_minus3\t P_anova_zero3\t P_anova_plus3\t Std_minus\t Std_zero\t Std_plus\t Std_minus2\t Std_zero2\t Std_plus2\t Std_minus3\t Std_zero3\t Std_plus3\t ANOVA_max\t ANOVA_spon\t ANOVA_maxmin\t ANOVA_maxspon\t ANOVA_max2\t ANOVA_spon2\t ANOVA_maxmin2\t ANOVA_maxspon2\t ANOVA_max3\t ANOVA_spon3\t ANOVA_maxmin3\t ANOVA_maxspon3\t Shift_ratio_plus\t Shift_ratio_minus\t Shift_ratio_plusminus\t Shift_ratio_plus2\t Shift_ratio_minus2\t Shift_ratio_plusminus2\t Shift_ratio_plus3\t Shift_ratio_minus3\t Shift_ratio_plusminus3\t Frame_plus\t Frame_minus\t Frame_plusminus\t Frame_plus2\t Frame_minus2\t Frame_plusminus2\t Frame_plus3\t Frame_minus3\t Frame_plusminus3\t Peak_minus\t Peak_zero\t Peak_plus\t Peak_minus2\t Peak_zero2\t Peak_plus2\t Peak_minus3\t Peak_zero3\t Peak_plus3\t Rsquared_minus\t Rsquared_zero\t Rsquared_plus\t Rsquared_minus2\t Rsquared_zero2\t Rsquared_plus2\t Rsquared_minus3\t Rsquared_zero3\t Rsquared_plus3\t P_fit_minus\t P_fit_zero\t P_fit_plus\t P_fit_minus2\t P_fit_zero2\t P_fit_plus2\t P_fit_minus3\t P_fit_zero3\t P_fit_plus3\t SR_fit_minus\t SR_fit_plusminus\t SR_fit_plus\t SR_fit_minus2\t SR_fit_plusminus2\t SR_fit_plus2\t SR_fit_minus3\t SR_fit_plusminus3\t SR_fit_plus3\t Num_trials\t Gaze_dir');
    %     fprintf(fid, '\r\n');
    % end
    % fprintf(fid, '%s', buff);
    % fprintf(fid, '\r\n');
    % fclose(fid);
    % 
    % % -------------------------------------------------------------------------
    % % Write out data relevant for model fitting analysis to separate file
    % 
    % buff2 = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %3.0f\t %s\t', ...
    %        FILE, peak(:), p_anova(:), Z(:), SR_true{:}, boot_data{:}, SR_fit(:), rsquared_head, rsquared_eye, partialZ_head, partialZ_eye, EndTrial-(BegTrial-1), dir_text);
    %        % currently 69 fields (8/19/05)
    % outfile2 = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\Dir1D_Az_VaryFix_Models.dat'];
    % printflag2 = 0;
    % if (exist(outfile2, 'file') == 0)    %file does not yet exist
    %     printflag2 = 1;
    % end
    % fid2 = fopen(outfile2, 'a');
    % if (printflag2)
    %     fprintf(fid2, 'FILE\t Peak_minus\t Peak_zero\t Peak_plus\t Peak_minus2\t Peak_zero2\t Peak_plus2\t Peak_minus3\t Peak_zero3\t Peak_plus3\t P_anova_minus\t P_anova_zero\t P_anova_plus\t P_anova_minus2\t P_anova_zero2\t P_anova_plus2\t P_anova_minus3\t P_anova_zero3\t P_anova_plus3\t Z_minus\t Z_zero\t Z_plus\t Z_minus2\t Z_zero2\t Z_plus2\t Z_minus3\t Z_zero3\t Z_plus3\t Shift_ratio_plus\t Shift_ratio_minus\t Shift_ratio_plusminus\t Shift_ratio_plus2\t Shift_ratio_minus2\t Shift_ratio_plusminus2\t Shift_ratio_plus3\t Shift_ratio_minus3\t Shift_ratio_plusminus3\t Frame_plus\t Frame_minus\t Frame_plusminus\t Frame_plus2\t Frame_minus2\t Frame_plusminus2\t Frame_plus3\t Frame_minus3\t Frame_plusminus3\t SR_fit_minus\t SR_fit_plusminus\t SR_fit_plus\t SR_fit_minus2\t SR_fit_plusminus2\t SR_fit_plus2\t SR_fit_minus3\t SR_fit_plusminus3\t SR_fit_plus3\t Rsquared_head1\t Rsquared_head2\t Rsquared_head3\t Rsquared_eye1\t Rsquared_eye2\t Rsquared_eye3\t partialZ_head1\t partialZ_head2\t partialZ_head3\t partialZ_eye1\t partialZ_eye2\t partialZ_eye3\t Num_trials\t Gaze_dir');
    %     fprintf(fid2, '\r\n');
    % end
    % fprintf(fid2, '%s', buff2);
    % fprintf(fid2, '\r\n');
    % fclose(fid2);
    % 
    % end % (proceed1 or proceed2)
    
    if proceed
    % -------------------------------------------------------------------------
    % DI's only
                    
    buff3 = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', ...
           FILE, DI);
    outfile3 = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\Dir1D_Az_VaryFix_DI.dat'];
    printflag3 = 0;
    if (exist(outfile3, 'file') == 0)    %file does not yet exist
        printflag3 = 1;
    end
    fid3 = fopen(outfile3, 'a');
    if (printflag3)
        fprintf(fid3, 'FILE\t DI_minus\t DI_plusminus\t DI_plus\t DI_minus2\t DI_plusminus2\t DI_plus2\t DI_minus3\t DI_plusminus3\t DI_plus3\t');
        fprintf(fid3, '\r\n');
    end
    fprintf(fid3, '%s', buff3);
    fprintf(fid3, '\r\n');
    fclose(fid3);
    
    
    % -------------------------------------------------------------------------
    % HDI's and F values (tuning strength), and congruency
    buff4 = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', ... 
                    FILE, congruency, F_val, HDI, SSE_total);
    outfile4 = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation\Dir1D_Az_VaryFix_HDI.dat'];
    printflag4 = 0;
    if (exist(outfile4, 'file') == 0)    %file does not yet exist
        printflag4 = 1;
    end
    fid4 = fopen(outfile4, 'a');
    if (printflag4)
        fprintf(fid4, 'FILE\t Congruency\t F_minus\t F_zero\t F_plus\t F_minus2\t F_zero2\t F_plus2\t F_minus3\t F_zero3\t F_plus3\t HDI_minus\t HDI_zero\t HDI_plus\t HDI_minus2\t HDI_zero2\t HDI_plus2\t HDI_minus3\t HDI_zero3\t HDI_plus3\t SSE_minus\t SSE_zero\t SSE_plus\t SSE_minus2\t SSE_zero2\t SSE_plus2\t SSE_minus3\t SSE_zero3\t SSE_plus3\t');
        fprintf(fid4, '\r\n');
    end
    fprintf(fid4, '%s', buff4);
    fprintf(fid4, '\r\n');
    fclose(fid4);
    
    end % (proceed)
    
    time = clock;
    disp('END TIME = ');
    disp(time(4:6));
    toc
    
    return;

end
end