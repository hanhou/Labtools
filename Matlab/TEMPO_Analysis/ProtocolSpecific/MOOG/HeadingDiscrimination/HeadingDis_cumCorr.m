%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cumCorr(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; 

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);
temp_spike_data = data.spike_data(1,:);   % spike rasters
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
spike_rates = temp_spike_rates( select_trials);
total_trials = temp_total_trials( select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
disc_heading = unique_heading( floor(length(unique_heading)/2)+1 : end );

%%%%%%%% PSTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );
spike_data( find(spike_data>100) ) = 1; % something is absolutely wrong

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
for k=1: length(unique_stim_type)
    count_con(1:100) = 0;
    for i=1:length(unique_heading)
        select = logical( (heading==unique_heading(i)) & (stim_type==unique_stim_type(k)) );            
        act_found = find( select==1 );
        % count spikes per timebin on every same condition trials
        temp_count = [];
        for repeat=1:length(act_found) 
            for n=1:(x_length)
                temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                time_step=time_step+timebin;
            end
            time_step=1;  
        end
        count_y_trial{i,k}(:,:) = temp_count;  % each trial's PSTH 
        dim{i,k} = size(temp_count);
        count_y{i,k} = mean(temp_count);  % take the mean
        max_count_y(i,k) = max(count_y{i,k});
        resp_mean(i,k) = sum(count_y{i,k}(26:55)); % mean firing rate  
      
    end 
%     % normalize PSTH to 1 for each stimulus condition
%     for i=1:length(unique_heading) 
%         count_y{i,k} = count_y{i,k} / max(max_count_y(:,k));
%     end  
end

% now calculate corr based on trials
t = 0:0.05:2;
ampl = 0.30;
num_sigs = 4;
pos = ampl*0.5*(erf(2*num_sigs/3*(t-1)) + 1);
veloc = diff(pos)/0.01;  % velocity has 40 points
veloc = veloc/max(veloc);  % normalize
accel = diff(veloc)/0.01; % accel has 39 points
accel = 2*(accel - min(accel))/max(accel - min(accel)); % make it positive and normalized, so that the middle line is 1, max is 2, min is 0

% decide shift_win for each cell
%for k = 1 : length(unique_stim_type)
for k = 1 : 1 % for vestibular temporarily
    pref = find( resp_mean(:,k) == max(resp_mean(:,k)) );
    nullf = find( resp_mean(:,k) == min(resp_mean(:,k)) );
    pre = pref(1); % in case two same peaks
    null = nullf(1); 
    %arbitrarily add 3 more bins after 67 to kill the saccade response,
    %which usually starts from 68 bin
    count_y{pre,k}(68:70) = count_y{pre,k}(67); 
    count_y{null,k}(68:70) = count_y{null,k}(67);
    count_y{pre,k}(21:66) = count_y{pre,k}(21:66)/max(count_y{pre,k}(21:66)); % normalize
    
    % for acceleration, stimulus profile need to be clipped to avoid
    % negative values
    max_peak = 0;
    for i = 21:60
        temp_peak = mean(count_y{pre,k}(i-1:i+1));% find peak of response (highest 3 bins instead of highest bin)
        if temp_peak > max_peak
            max_peak = temp_peak;
            peak_locate = i;
        end
    end
    if max_peak > 2*mean(count_y{pre,k}(1:20)) 
        ratio = mean(count_y{pre,k}(1:20)) / ( max_peak - mean(count_y{pre,k}(1:20)) ) ;
        accel(accel < (1-ratio)) = 1-ratio;
    end
    accel(40) = accel(39); % make accel 40 points also

    for s = 1 : 7 % move 6 bin back, which are 300ms each direction 
        cc = corrcoef( veloc, count_y{pre,k}(21+s-1 : 60+s-1) ); % vel
        aa_pre(s) = cc(1,2);
        cc = corrcoef( veloc, count_y{null,k}(21+s-1 : 60+s-1) ); %vel
        aa_null(s) = cc(1,2);
        acce = corrcoef( accel, count_y{pre,k}(21+s-1 : 60+s-1) ); % acceleration
        acc(s) = acce(1,2);
    end
    % if more negatively corelated than positively correlated, use negative
    % correlation then
    pp = mean(count_y{pre,k}(26:55)) - mean(count_y{pre,k}(1:20));
    nn = mean(count_y{null,k}(26:55)) - mean(count_y{null,k}(1:20));
%     if abs(nn) > abs(pp)  
%         shift_win_vel = find(aa_null==min(aa_null));
%         Rv = max( aa_null );     
%     else
        shift_win_vel = find(aa_pre==max(aa_pre));
        Rv = max( aa_pre);
%    end
    % for acceleration
    shift_win_acc = find( abs(acc)==max(abs(acc)) ); % could be negatively correlated   
    % correlation between models
    vel_model = zeros(1,50);
    acc_model = ones(1,50);
    vel_model(shift_win_vel : 39+shift_win_vel) = veloc;
    acc_model(shift_win_acc : 39+shift_win_acc) = accel;
    coef_mod = corrcoef(vel_model, acc_model);
    Rva = coef_mod(1,2);
    Ra = max(acc);
    % now calculate partial correlation
    Rv_p = (Rv-Ra*Rva) / sqrt( (1-Ra^2)*(1-Rva^2) ); % partial correlation for velocity
    Ra_p = (Ra-Rv*Rva) / sqrt( (1-Rv^2)*(1-Rva^2) ); % partial correlation for accel
    % Z-transform 
    Zv = 0.5*log( (1+Rv_p)/(1-Rv_p) ) / sqrt(1/37);
    Za = 0.5*log( (1+Ra_p)/(1-Ra_p) ) / sqrt(1/37);
%     if max_allow == 11;
%         for h = 1: 6 %shift delay from 0 to 500ms with 100ms step
%             Rv_p_fix = (Rv_fix(2*h-1)-Ra_fix(2*h-1)*Rva) / sqrt( (1-Ra_fix(2*h-1)^2)*(1-Rva^2) ); % partial correlation for velocity
%             Ra_p_fix = (Ra_fix(2*h-1)-Rv_fix(2*h-1)*Rva) / sqrt( (1-Rv_fix(2*h-1)^2)*(1-Rva^2) ); % partial correlation for accel
%             Zv_fix(h) = 0.5*log( (1+Rv_p_fix)/(1-Rv_p_fix) ) / sqrt(1/37);
%             Za_fix(h) = 0.5*log( (1+Ra_p_fix)/(1-Ra_p_fix) ) / sqrt(1/37);    
%             Rv_fixx(h) = Rv_fix(2*h-1);
%             Ra_fixx(h) = Ra_fix(2*h-1);
%         end
%     end
end    

%%%%%%% replace spike_rates with corr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save a copy for spike_rates first
spike_rates_copy = spike_rates;
for nn = 1 : 2 % either with velocity or acceleration
    
    for k = 1 : length(unique_stim_type)
        for i = 1 : length(unique_heading)
            select = logical(heading==unique_heading(i) & stim_type==unique_stim_type(k));
            act_find = find(select==1);
            for j = 1 : length(act_find)
                count_y_trial{i,k}(j,68:70) = count_y_trial{i,k}(j,67);
                if nn == 1
                   count_corr = count_y_trial{i,k}(j,21-1+shift_win_vel : 60-1+shift_win_vel );
                   cccc = sum(count_corr.*veloc) - sum(count_corr)*sum(veloc)/40; % this is the numerator of the Corrcoef
                else
                   count_corr = count_y_trial{i,k}(j,21-1+shift_win_acc : 60-1+shift_win_acc );
                   cccc = sum(count_corr.*accel) - sum(count_corr)*sum(accel)/40;
                end
                 corr_trial{i,k}(j) = cccc;
                 bb(j) = cccc;
                 spike_rates(act_find(j)) = corr_trial{i,k}(j); 
            end
            corr_mean(i,k) = mean(bb); %take the mean for each heading
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %neurometric dataset and calculate ROC, Choice Probability(CP)
    %determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
    LEFT = 1;
    RIGHT = 2;
    for i= 1 : length(spike_rates) 
        temp = data.event_data(1,:,i + BegTrial-1);
        events = temp(temp>0);  % all non-zero entries
        if (sum(events == IN_T1_WIN_CD) > 0)
            choice(i) = RIGHT;
        elseif (sum(events == IN_T2_WIN_CD) > 0)
            choice(i) = LEFT;
        else
            disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
        end
    end

    one_repetition = length(unique_heading)*length(unique_stim_type);
    repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition
    resp_heading = [];
    Z_Spikes = spike_rates;
    % z-score data for later cp analysis across headings
    for k = 1:length(unique_stim_type)
        for i = 1:length(unique_heading)
            select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
            z_dist = spike_rates(select);
            z_dist = (z_dist - mean(z_dist))/std(z_dist);
            Z_Spikes(select) = z_dist;
        end
    end
    Z_Spikes_Ori = Z_Spikes; % keep a Z_Spikes unchanged for later use
    % now group neuronal data into two groups according to monkey's choice
    for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading    
        for i = 1:length(unique_heading)
            select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
            resp{k,i} = spike_rates(select);

            resp_mat_copy{k}(i) = mean(spike_rates_copy(select)); % metric based on mean firing rate
            resp_mat{k}(i) = mean(resp{k,i});     % use the DFT for each trial     
            resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);

            % calculate CP, group data based on monkey's choice 
            resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );
            resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );
            if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
                Z_Spikes(select) = 9999;   % similar to NaN, just make a mark            
            end 
        end  
        % now cross all data 
        resp_left_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
        resp_right_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
        resp_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) ); 
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [rr_DFT_trial,pp_DFT_trial] = corrcoef(unique_heading, resp_mat{1}(:));  % slope based on DFT trials
    line_re_DFT_trial = rr_DFT_trial(1,2);
    line_p_DFT_trial = pp_DFT_trial(1,2);
    line_r(nn) = line_re_DFT_trial;
    line_p(nn) = line_p_DFT_trial;
  
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
    % neurothreshold 
    fit_data_neuro = [];
    fit_data_neuro_cut = [];
    for i = 1 : length(unique_heading)-1   % subtract the 0 heading
        trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(1)) ) ;
        fit_data_neuro_cum{nn}(i,3) = sum(trials_n);  % for later function fit use
        if i < (1+length(unique_heading))/2
             % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction            
             Neuro_correct(i) =  rocN( resp{1,length(unique_heading)-i+1},resp{1,i},100 ); 
         else
             Neuro_correct(i) =  rocN( resp{1,length(unique_heading)-i}, resp{1,(i+1)},100 );  
         end
         if line_re_DFT_trial > 0 % if left heading is not the larger responses, then linear regression will be positive 
             % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
             % here we asume if left response larger than right response then asign the left to be preferred direction
             Neuro_correct(i) = 1 - Neuro_correct(i);            
         end  
     end
     % choice probability
    for i = 1 : length(unique_heading)  
        if (length(resp_left_all{1}) > 3) | (length(resp_right_all{1}) > 3)
            CP_all{nn} = rocN( resp_left_all{1},resp_right_all{1},100 );
        else
            CP_all{nn} = NaN; 
        end
        if  line_re_DFT_trial > 0
            CP_all{nn} = 1 - CP_all{nn};
        end
    end
 
    %%%%%% use Wichman's MLE method to estimate threshold and bias  
    fit_data_neuro_cum{nn}(:,1) = unique_heading(unique_heading~=0);
    fit_data_neuro_cum{nn}(:,2) = Neuro_correct(:);
    wichman_neu = pfit(fit_data_neuro_cum{nn},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
    Thresh_neu{nn} = wichman_neu.params.est(2);
    % negative and positive infinite value means flat tuning
    if Thresh_neu{nn}<0 | Thresh_neu{nn}> 300
        Thresh_neu{nn} = 300;
    end  

      % permutation test for CP
      %do permutation to test the significance of CP_all{k}, re-calculate CP 2000 times
    perm_num = 1000;
    Z_Spikes_perm = Z_Spikes;
    bin = 0.005;
    x_bin = 0 : bin : 1;
    for k = 1:1    
        for n = 1 : perm_num
            % temperarilly only use near-threshold heading angles where monkey make a guess mainly
            select = logical( (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
            Z_Spikes_con{k} = Z_Spikes_perm( select );
            Z_Spikes_con{k} = Z_Spikes_con{k}(randperm(length(Z_Spikes_con{k})));   % permute spike_rates
            Z_Spikes_perm(select) = Z_Spikes_con{k};    % now in spike_rates, the corresponding data were permuted already

            resp_left_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
            resp_right_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 

            if  (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
                CP_all_perm{k}(n) = rocN( resp_left_all_perm{k}, resp_right_all_perm{k},100 );
            else
                CP_all_perm{k}(n) = NaN; 
            end
            if  line_re_DFT_trial > 0  
                CP_all_perm{k}(n) = 1 - CP_all_perm{k}(n);             
            end  

            resp_left_choose_perm = Z_Spikes_perm( (heading == unique_heading(5)) & (stim_type == unique_stim_type(1)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
            resp_right_choose_perm = Z_Spikes_perm( (heading == unique_heading(5)) & (stim_type == unique_stim_type(1)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 

            if  (length(resp_left_choose{1,5}) > 3) & (length(resp_right_choose{1,5}) > 3)
                CP_perm(n) = rocN( resp_left_choose_perm, resp_right_choose_perm,100 );
            else
                CP_perm(n) = NaN; 
            end
            if  line_re_DFT_trial > 0  
                CP_perm(n) = 1 - CP_perm(n);             
            end  

        end
        % now calculate p value or significant test
        if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
            hist_perm(k,:) = hist( CP_all_perm{k}(:), x_bin );  % for permutation
            bin_sum = 0;
            n = 0;
            while ( n < (CP_all{k}/bin) )
                 n = n+1;
                 bin_sum = bin_sum + hist_perm(k, n);
                 if CP_all{k} > 0.5                  % note it's two tail test
                    p{k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
                 else
                    p{k} = 2* bin_sum / perm_num;
                 end
            end
        else
            p{k} = NaN;
        end 
    end 
    pp(nn) = p{k};
end
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
buff = sprintf(sprint_txt, FILE, Thresh_neu{1}, Thresh_neu{2}, CP_all{1}, CP_all{2}, pp(1),pp(2), line_r(:), line_p(:) );
%buff = sprintf(sprint_txt, FILE, Zv, Za, Thresh_neu{1}, CP_all{1},Thresh_neu{2}, CP_all{2} );
%buff = sprintf(sprint_txt, FILE, shift_win_vel, shift_win_acc );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cumCorr.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         theshold\t CP\t Rfiring\t Rdfttrial\t Rdftmean\t phase\t ratio\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
return;