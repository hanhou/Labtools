function [Thresh_neu,CP_all]=ChoiceProbability_simp(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OptimalWindow);

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
temp_spike_rates = data.spike_rates(SpikeChan, :); 

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
disc_heading = unique_heading( floor(length(unique_heading)/2)+1 : end );%???

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );
spike_data( find(spike_data>100) ) = 1; % something is absolutely wrong 
spike_rates_copy = spike_rates; % save a copy for spike_rates

%Calculate the choice probability during the optimal window
for ss =  1 : length(spike_rates) % ss marks the index of trial
    spike_rates(ss) = sum( spike_data(1,min(OptimalWindow)+5000*(ss-1) : max(OptimalWindow)+5000*(ss-1)) ) / 1.01; % 996~3006 every 1000ms
end

% monkey's choice:neurometric dataset and calculate ROC, Choice
% Probability(CP) determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2) 
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

% neural database
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
        resp_mat{k}(i) = mean(resp{k,i});  % the mean firing rate for each heading 
        resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);
        % calculate CP, group data based on monkey's choice 
        resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );
        resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );
        if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
            %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)  
            Z_Spikes(select) = 9999;   % similar to NaN, just make a mark            
            %     Z_Spikes( (heading == unique_heading(length(unique_heading)+1-i)) & (stim_type == unique_stim_type(k)) ) = 9999; % the corresponding heading is also excluded
            CP{k}(i) = 9999;
            %      CP{k}(length(unique_heading)+1-i) = 9999;
        else
            CP{k}(i) = 0;
        end
    end
    
    % now across all data
    resp_left_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
    resp_right_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
    resp_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );     
end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % decide sign of slope, decide preferred direction 
% % % %--------- plot neuronal response aross time, see whether it is stable--------------
% figure(9);
% kk{1}='b.-'; kk{2}='r.-'; kk{3}='k.-'; kk{4}='g.-'; kk{5}='y.'; kk{6}='y.-'; 
% nu{1}='b.--'; nu{2}='r.--'; nu{3}='k.--'; nu{4}='g.--'; nu{5}='y.--'; nu{6}='y.--';  
% plot(resp{k,1},kk{1});hold on;
% plot(resp{k,5},kk{4});hold on;plot(resp{k,9}, kk{2});hold on;ylim([0 100]);
% set(gca, 'ytick',[0:20:100]);set(gca, 'xtick',[0 repetition]);
% 
for k = 1 : length(unique_stim_type)
    % decide whether ves and vis is congruent tuning. Fit line by linear
    % regression first and compare the sign of each condition to decide whether
    % congruent or opposite, this is used to check whether congruent cells lead
    % to better neuronal performance in combined condition, and vice versa
    [rr,pp] = corrcoef(unique_heading, resp_mat{k}(:));
    line_re{k} = rr(1,2);
    line_p{k} = pp(1,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
% neurothreshold 
fit_data_neuro = [];
fit_data_neuro_cut = [];       
for k = 1 : length(unique_stim_type)
    fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
    for i = 1 : length(unique_heading)-1   % subtract the 0 heading
        trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;            
        if i < (1+length(unique_heading))/2
            % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction            
            Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i+1},resp{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
        else
            Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i}, resp{k,(i+1)},100 ); % compare to the 0 heading condition, which is straght ahead
        end
        %        if  resp_mat{k}(1) < resp_mat{k}(end)
        if line_re{k} > 0 % if left heading is not the larger responses, then linear regression will be positive 
            % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
            % here we asume if left response larger than right response then asign the left to be preferred direction
            Neuro_correct{k}(i) = 1 - Neuro_correct{k}(i);            
        end
        fit_data_neuro_cum{k}(i,2) = Neuro_correct{k}(i);
        fit_data_neuro_cum{k}(i,3) = sum(trials_n);  % for later function fit use
    end
end

% choice probability
for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_heading)  
        if CP{k}(i)~=9999
            CP{k}(i) = rocN( resp_left_choose{k,i},resp_right_choose{k,i},100 );
        else
            CP{k}(i) = NaN;
        end
        if (length(resp_left_all{k}) > 3) | (length(resp_right_all{k}) > 3)
            %if  (length(resp_left_all{k}) / length(resp_all{k}) >0.25) |  (length(resp_left_all{k}) / length(resp_all{k}) <0.75)
            CP_all{k} = rocN( resp_left_all{k},resp_right_all{k},100 );
        else
            CP_all{k} = NaN;
        end
        if  line_re{k} > 0
            CP{k}(i) = 1 - CP{k}(i);
            CP_all{k} = 1 - CP_all{k};
        end
    end
end

%%%%%% use Wichman's MLE method to estimate threshold and bias
for k = 1:length(unique_stim_type)        
    [bb,tt] = cum_gaussfit_max1(fit_data_neuro_cum{k});
    Thresh_neu{k} = tt;
    % negative and positive infinite value means flat tuning
    if Thresh_neu{k}<0 | Thresh_neu{k}> 300
        Thresh_neu{k} = 300;
        %                wichman_neu.params.est(2) = 300;
    end
end

for k = 1:length(unique_stim_type)    
    wichman_psy = pfit(fit_data_psycho_cum{k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
    Thresh_psy{k} = wichman_psy.params.est(2);
    Bias_psy{k} = wichman_psy.params.est(1);
    psy_perf{k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
    fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
    fit_data_neuro_cum{k}(:,2) = Neuro_correct{k}(:);
    wichman_neu = pfit(fit_data_neuro_cum{k}(:,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
    Thresh_neu{k} = wichman_neu.params.est(2);
    % negative and positive infinite value means flat tuning
    if Thresh_neu{k}<0 | Thresh_neu{k}> 300
        Thresh_neu{k} = 300;
        wichman_neu.params.est(2) = 300;
    end
    Bias_neu{k} = wichman_neu.params.est(1);
    neu_perf{k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
end

% %--------------------------------------------------------------------------
% %do permutation to test the significance of CP_all{k}, re-calculate CP 2000 times
% perm_num = 1000;
% Z_Spikes_perm = Z_Spikes;
% bin = 0.005;
% x_bin = 0 : bin : 1;
% for k = 1:length(unique_stim_type)    
%     for n = 1 : perm_num
%         % temperarilly only use near-threshold heading angles where monkey make a guess mainly
%         select = logical( (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
%         Z_Spikes_con{k} = Z_Spikes_perm( select );
%         Z_Spikes_con{k} = Z_Spikes_con{k}(randperm(length(Z_Spikes_con{k})));   % permute spike_rates
%         Z_Spikes_perm(select) = Z_Spikes_con{k};    % now in spike_rates, the corresponding data were permuted already
%         
%         resp_left_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
%         resp_right_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
%    
%         if  (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
%             CP_all_perm{k}(n) = rocN( resp_left_all_perm{k}, resp_right_all_perm{k},100 );
%         else
%             CP_all_perm{k}(n) = NaN; 
%         end
%         if  line_re{k} > 0  
%             CP_all_perm{k}(n) = 1 - CP_all_perm{k}(n);             
%         end  
%         
%         resp_left_choose_perm = Z_Spikes_perm( (heading == unique_heading(5)) & (stim_type == unique_stim_type(1)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
%         resp_right_choose_perm = Z_Spikes_perm( (heading == unique_heading(5)) & (stim_type == unique_stim_type(1)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
%    
%         if  (length(resp_left_choose{1,5}) > 3) & (length(resp_right_choose{1,5}) > 3)
%             CP_perm(n) = rocN( resp_left_choose_perm, resp_right_choose_perm,100 );
%         else
%             CP_perm(n) = NaN; 
%         end
%         if  line_re{1} > 0  
%             CP_perm(n) = 1 - CP_perm(n);             
%         end  
%         
%     end
%     % now calculate p value or significant test
%     if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
%   
%         hist_perm(k,:) = hist( CP_all_perm{k}(:), x_bin );  % for permutation
%         bin_sum = 0;
%         n = 0;
%         while ( n < (CP_all{k}/bin) )
%              n = n+1;
%              bin_sum = bin_sum + hist_perm(k, n);
%              if CP_all{k} > 0.5                  % note it's two tail test
%                 p{k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
%              else
%                 p{k} = 2* bin_sum / perm_num;
%              end
%         end
%     else
%         p{k} = NaN;
%     end 
%     
%     % calculate p value for CP during straight ahead motion
%     if (length(resp_left_choose{1,round(length(unique_heading)/2)}) > 3) & (length(resp_right_choose{1,round(length(unique_heading)/2)}) > 3)  
%         hist_perm(k,:) = hist( CP_perm(:), x_bin );  % for permutation
%         bin_sum = 0;
%         n = 0;
%         while ( n < (CP{1}(round(length(unique_heading)/2))/bin) )
%              n = n+1;
%              bin_sum = bin_sum + hist_perm(k, n);
%              if CP{k}(round(length(unique_heading)/2)) > 0.5                  % note it's two tail test
%                 pp = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
%              else
%                 pp = 2* bin_sum / perm_num;
%              end
%         end
%     else
%         pp = NaN;
%     end  
%     p_a(k) = pp;
% end

% h_title{1}='Vestibular';
% h_title{2}='Visual';
% h_title{3}='Combined';

% %plot psychometric and neurometric function here
% h{1} = 'bo';  f{1} = 'b-';  g{1} = 'bo-';
% h{2} = 'rd';  f{2} = 'r-';  g{2} = 'rd-';
% h{3} = 'gs';  f{3} = 'g-';  g{3} = 'gs-';
% figure(10);
% set(10,'Position', [5,25, 980,650], 'Name', 'psycho_neurometic function');
% orient landscape;
% %plot psychometric function
% axes('position',[0.05 0.3, 0.26 0.45]);
% title('psychometric');
% legend_txt = [];
% for k = 1:length(unique_stim_type)
%     xi = min(unique_heading) : 0.1 : max(unique_heading);   
%     beta = [0, 1.0];
%  %   plot data in logarithmical space instead of linspace
%     plot(unique_heading, psycho_correct(k,:), h{k}, xi, cum_gaussfit(psy_perf{k}, xi),  f{k} );
%     xlabel('Heading Angles');   
%     ylim([0,1]);
%     ylabel('Rightward Choices');
%     hold on;
%     legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
%     legend_txt{k*2} = [''];
% %    also fit data with weibull function
% %    [psycho_alpha(k) psycho_beta(k)]= weibull_fit(fit_data_psycho{k});
% end
% 
% %plot neurometric function
% axes('position',[0.36 0.3, 0.26 0.45]);
% title('neuroometric');
% for k = 1:length(unique_stim_type)
%     neu_heading = unique_heading(unique_heading~=0);
%     xi = min(unique_heading) : 0.1 : max(unique_heading); 
%     plot(neu_heading, Neuro_correct{k}, h{k},  xi, cum_gaussfit(neu_perf{k}, xi),  f{k} );
%     xlabel('Heading Angles');   
%     ylim([0,1]);
%     hold on;
% end
% 
% %%%%%%  neurological raw data based on firing rate instead of ROC
% axes('position',[0.7 0.3, 0.26 0.45]);
% title('firing rate');
% for k = 1:length(unique_stim_type)
%     errorbar(unique_heading, resp_mat{k}(:), resp_mat_err{k}(:),g{k} );
%     xlabel('Heading Angle (deg)');
%     ylabel('Firing rate(spikes/s)');   
%     xlim([min(unique_heading),max(unique_heading)]);
%     hold on;
% end
% 
% % output some text of basic parameters in the figure
% axes('position',[0.05,0.8, 0.9,0.15] );
% xlim( [0,100] );
% ylim( [2,10] );
% text(0, 10, FILE);
% text(20,10,'coherence =');
% text(30,10,num2str(unique_motion_coherence) );
% text(45,10,'repetition =');
% text(55,10,num2str(repetition) ); 
% text(5,8, 'Psy: u      threshold        err           Neu:u         threshold         err              CP        p');
% text(0,8, 'stim');
% for k = 1:length(unique_stim_type)
%     text(0,8-k, num2str(unique_stim_type(k)));
%     text(5,8-k,num2str(Bias_psy{k} ));
%     text(12,8-k,num2str(Thresh_psy{k} ));
%     text(30,8-k,num2str(Bias_neu{k} ));
%     text(40,8-k,num2str(Thresh_neu{k} ));
%     text(50,8-k,num2str(CP_all{k}') ); 
%     text(60,8-k,num2str(p{k}') );   
% end
% axis off;