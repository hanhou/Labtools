%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
% fluctuation of behavior and neuronal sensitivity over time within one
% session (not individual trial!). We want to see whether there is a
% correlation of the two fluctuation (Zohary et al., 1994) (not choice probability!) and
% more importantly, the cue integration effect in our own task
%--	03/17/08 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cum_session(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);
temp_spike_data = data.spike_data(SpikeChan,:);   % spike rasters
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part deals with for whatever reason, some trials are absolutely
% wrong, like >1000 spikes/s, makes no sense. tick this bad trial out,
% replace with one of the other data within same condition, this is similar
% to bootstrap. But in order to keep every running the same, use the mean
% of the rest instead. this happen only very rarely!!!
bad_tri = find(temp_spike_rates > 1000 );
if length(bad_tri) > 0
    for i=1:length(bad_tri)
		same_tri = find( (stim_type == stim_type(bad_tri(i))) & (heading == heading(bad_tri(i))) ); % find other trials within same condition
		rest_tri = same_tri(same_tri ~= bad_tri(i));
		spike_rates(bad_tri(i)) = mean( spike_rates(rest_tri) );
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
        select = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
        z_dist = spike_rates(select);
        z_dist = (z_dist - mean(z_dist))/std(z_dist);
        Z_Spikes(select) = z_dist;
    end
end
%choice(889) =2; % for cell m2c384r2 % for some reason the choice is 0 for this trial
% psychometric dataset
psycho_correct = [];
fit_data_psycho = [];
N_obs = [];
for i = 1:length(unique_heading)
    for k = 1:length(unique_stim_type)
         trials_p =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
         % make 'S' curve by using the rightward choice for y-axis
         correct_trials = (trials_p & (choice == RIGHT) );
         psycho_correct(k,i) = 1*sum(correct_trials) / sum(trials_p); 
         fit_data_psycho_cum{k}(i, 1) = unique_heading( i );  
         fit_data_psycho_cum{k}(i, 2) = psycho_correct(k,i);
         fit_data_psycho_cum{k}(i, 3) = sum(trials_p); 
     end
end

one_repetition = length(unique_heading)*length(unique_stim_type);
repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition
resp_heading = [];

for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading    
    for i = 1:length(unique_heading)
        select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
        resp{k,i} = spike_rates(select);   
        % calculate firing rate of each trial
%         for j = 1 : repetition; 
%             spike_temp = spike_rates(select);   
%             resp_heading{k}(i, j) = spike_temp(j);           
%         end
        resp_mat{k}(i) = mean(resp{k,i});  % the mean firing rate for each heading 
        resp_mat_std{k}(i)= std(resp{k,i});
        resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);
        % calculate CP, group data based on monkey's choice 
        resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );
        resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );
        if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
      %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)  
            Z_Spikes(select) = 9999;   % similar to NaN, just make a mark            
        end 
    end  
    % now across all data 
    resp_left_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
    resp_right_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
    resp_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) ); 
    
    [rr,pp] = corrcoef(unique_heading, resp_mat{k}(:));
    line_re{k} = rr(1,2);
    line_p{k} = pp(1,2);
end

% now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
% neurothreshold 
fit_data_neuro = [];
fit_data_neuro_cut = [];
for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_heading)-1   % subtract the 0 heading   
        trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
        fit_data_neuro_cum{k}(i,3) = sum(trials_n);  % for later function fit use
        if i < (1+length(unique_heading))/2
             Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i+1},resp{k,i},100 );
        else
             Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i}, resp{k,(i+1)},100 );
        end
         if line_re{k} > 0 % if left heading is not the larger responses, then linear regression will be positive 
             % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
             % here we asume if left response larger than right response then asign the left to be preferred direction
             Neuro_correct{k}(i) = 1 - Neuro_correct{k}(i);            
         end  
     end
end

% % choice probability
% for k = 1 : length(unique_stim_type)
%     if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3)
%         CP_all{k} = rocN( resp_left_all{k},resp_right_all{k},100 );
%     else
%         CP_all{k} = NaN; 
%     end
%     if  line_re{k} > 0
%         CP_all{k} = 1 - CP_all{k};
%     end
%  
%     % now bootstrap to detect 95% confidence interval of CPs    
%     for n=1:5000
%         for ll = 1:length(resp_left_all{k})
%             randselect = floor( length(resp_left_all{k})*rand+1 );
%             resp_left_bootstp(ll) = resp_left_all{k}(randselect);
%         end
%         for rr = 1:length(resp_right_all{k})
%             randselect = floor( length(resp_right_all{k})*rand+1 );
%             resp_right_bootstp(rr) = resp_right_all{k}(randselect);
%         end
%         CP_bootstp(k,n) = rocN( resp_left_bootstp(1:ll),resp_right_bootstp(1:rr),100 );
%         if  line_re{k} > 0
%             CP_bootstp(k,n) = 1 - CP_bootstp(k,n);
%         end        
%     end
%     % now calculate upper and lower bound for 95% confidence interval
%     CP_bootstp_sort = sort(CP_bootstp(k,:)); % sort CP 
%     CP_lowbound(k) = CP_bootstp_sort(5000*0.025); % 2.5% confidence on the left side
%     CP_upbound(k) = CP_bootstp_sort(5000-5000*0.025); % 2.5% confidence on the right side
%     CP_mean(k) = mean( CP_bootstp(k,:) );
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot psychometric, neurometric, CP over time
% % run the slide threshold over time, see whether performance fluctuate across time
% span = floor(repetition/2);  % calculate threshod every 5 repeats;
% slide = floor(repetition/2);  % slide threshod with increment of 1 repeats;
% BegTrial_shift = BegTrial;
% EndTrial_shift = BegTrial_shift + span*one_repetition-1;
% n=0;
% while EndTrial_shift <= EndTrial
%     n = n + 1;
%     select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) );
%     stim_type_shift = temp_stim_type( select_trials_shift );
%     heading_shift = temp_heading( select_trials_shift );
%     unique_stim_type_shift = munique(stim_type_shift');
%     unique_heading_shift = munique(heading_shift');
%     total_trials_shift = temp_total_trials( select_trials_shift);
%     spike_rates_shift = temp_spike_rates( select_trials_shift );
% 
%     for k = 1:length(unique_stim_type)
%         for i = 1:length(unique_heading)
%              trials_shift =logical( (heading_shift == unique_heading(i)) & (stim_type_shift == unique_stim_type(k)) ) ;
%              correct_trials_shift = (trials_shift & (total_trials_shift == RIGHT) );
%              
%              resp_heading_shift{k,i} = spike_rates_shift(trials_shift );
%              % make 'S' curve by using the rightward choice for y-axis
%              if ( unique_heading(i) < 0 )
%                  correct_rate_shift(i) = 1 - 1*sum(correct_trials_shift) / sum(trials_shift); 
%              else
%                  correct_rate_shift(i) = 1*sum(correct_trials_shift) / sum(trials_shift); 
%              end         
%          end
%          fit_data_psycho_cum_shift{k}(:, 1) = fit_data_psycho_cum{k}(:, 1);  
%          fit_data_psycho_cum_shift{k}(:, 2) = correct_rate_shift(:);
%          fit_data_psycho_cum_shift{k}(:, 3) = span;
%          [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum_shift{k});
%          psy_thresh_shift(k,n) = tt;
%          % for neuronal performence over time
%          for i = 1 : length(unique_heading)-1   % subtract the 0 heading
%              if i < (1+length(unique_heading))/2
%                  Neuro_correct_shift{k}(i) =  rocN( resp_heading_shift{k,(1+length(unique_heading))/2},resp_heading_shift{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
%              else
%                  Neuro_correct_shift{k}(i) =  rocN( resp_heading_shift{k,(1+length(unique_heading))/2},resp_heading_shift{k,i+1},100 ); % compare to the 0 heading condition, which is straght ahead
%              end
%              if  resp_mat{k}(1) < resp_mat{k}(end)  
%                  Neuro_correct_shift{k}(i) = 1 - Neuro_correct_shift{k}(i);            
%              end  
%          end
%          fit_data_neu_cum_shift{k}(:, 1) = unique_heading(unique_heading~=0);  
%          fit_data_neu_cum_shift{k}(:, 2) = Neuro_correct_shift{k}(:);
%          fit_data_neu_cum_shift{k}(:, 3) = span;
%          [bbb,ttt] = cum_gaussfit_max1(fit_data_neu_cum_shift{k});
%          if ttt>300
%              ttt=300;
%          end
%          neu_thresh_shift(k,n) = ttt;           
%     end   
%     BegTrial_shift = BegTrial_shift + slide*one_repetition;
%     EndTrial_shift = BegTrial_shift + span*one_repetition-1;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot psychometric, neurometric every single repetition!!
% span = 1;  
% slide = 1;  
% BegTrial_shift = BegTrial;
% EndTrial_shift = BegTrial_shift + span*one_repetition-1;
% n=0;
% while EndTrial_shift <= EndTrial
%     n = n + 1;
%     select_trials_shift = ( (trials >= BegTrial_shift) & (trials <= EndTrial_shift) );
%     stim_type_shift = temp_stim_type( select_trials_shift );
%     heading_shift = temp_heading( select_trials_shift );
%     unique_stim_type_shift = munique(stim_type_shift');
%     unique_heading_shift = munique(heading_shift');
%     total_trials_shift = temp_total_trials( select_trials_shift);
%     spike_rates_shift = temp_spike_rates( select_trials_shift );
% 
%     for k = 1:length(unique_stim_type)
%         % for psychometric
%         for i = 1:length(unique_heading)
%              trials_shift =find( (heading_shift == unique_heading(i)) & (stim_type_shift == unique_stim_type(k)) ) ;
%                           
%              resp_heading_shift{k,i} = spike_rates_shift(trials_shift );
%              % make 'S' curve by using the rightward choice for y-axis
%              if ( unique_heading(i) < 0 )
%                  correct_rate_shift(i) = 1 - 1*sum(correct_trials_shift) / sum(trials_shift); 
%              else
%                  correct_rate_shift(i) = 1*sum(correct_trials_shift) / sum(trials_shift); 
%              end   
%              
%          end
%                   
%          % for neuronal performence over time
%          for i = 1 : length(unique_heading)-1   % subtract the 0 heading
%              if i < (1+length(unique_heading))/2
%                  Neuro_correct_shift{k}(i) =  rocN( resp_heading_shift{k,(1+length(unique_heading))/2},resp_heading_shift{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
%              else
%                  Neuro_correct_shift{k}(i) =  rocN( resp_heading_shift{k,(1+length(unique_heading))/2},resp_heading_shift{k,i+1},100 ); % compare to the 0 heading condition, which is straght ahead
%              end
%              if  resp_mat{k}(1) < resp_mat{k}(end)  
%                  Neuro_correct_shift{k}(i) = 1 - Neuro_correct_shift{k}(i);            
%              end  
%          end
%          fit_data_neu_cum_shift{k}(:, 1) = unique_heading(unique_heading~=0);  
%          fit_data_neu_cum_shift{k}(:, 2) = Neuro_correct_shift{k}(:);
%          fit_data_neu_cum_shift{k}(:, 3) = span;
%          [bbb,ttt] = cum_gaussfit_max1(fit_data_neu_cum_shift{k});
%          if ttt>300
%              ttt=300;
%          end
%          neu_thresh_shift(k,n) = ttt;           
%     end   
%     BegTrial_shift = BegTrial_shift + slide*one_repetition;
%     EndTrial_shift = BegTrial_shift + span*one_repetition-1;
% end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 1000 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
%buff = sprintf(sprint_txt, FILE, repetition, psy_thresh_shift(1,:),psy_thresh_shift(2,:),psy_thresh_shift(3,:), neu_thresh_shift(1,:),neu_thresh_shift(2,:),neu_thresh_shift(3,:) );
buff = sprintf(sprint_txt, FILE, CP_lowbound(1),CP_all{1},CP_upbound(1), CP_mean(1),CP_lowbound(2),CP_all{2},CP_upbound(2),CP_mean(2),CP_lowbound(3),CP_all{3},CP_upbound(3),CP_mean(3) );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum_CPbootstp.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         Coher\t repet\t Ve_P_u\t Vi_P_u\t co_P_u\t Ve_P_thr\t vi_P_th\t Co_P_thr\t Ve_N_u\t Vi_N_u\t co_N_u\t Ve_N_thr\t vi_N_th\t Co_N_thr\t Ves_CP\t Vis_CP\t Com_CP\t Ves_p\t Vis_p\t Com_p\t sign\t vesMax\t visMax\t coMax\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
return;