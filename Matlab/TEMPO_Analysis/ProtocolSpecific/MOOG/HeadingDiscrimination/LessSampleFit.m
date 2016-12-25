%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
% remove data at largest heading, compare resultant threshold both behavior
% and neuronal
%--	09/20/06 GY
%-----------------------------------------------------------------------------------------------------------------------

function LessSampleFit(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    

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
% if FILE=='m2c384r2.htb' %choice(889) =2; % for cell m2c384r2 % for some reason the choice is 0 for this trial
%      choice(889) = 2;
% end

% psychometric dataset
psycho_correct = [];
fit_data_psycho = [];
N_obs = [];
for i = 1:length(unique_heading)
     trials_p =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(1)) ) ;
     % make 'S' curve by using the rightward choice for y-axis
     correct_trials = (trials_p & (choice == RIGHT) );
     psycho_correct(1,i) = 1*sum(correct_trials) / sum(trials_p); 
     fit_data_psycho_cum(i, 1) = unique_heading( i );  
     fit_data_psycho_cum(i, 2) = psycho_correct(1,i);
     fit_data_psycho_cum(i, 3) = sum(trials_p);          
end

% remove 2 data points (1 pair at edge)
fit_data_psycho_2(1:7, 1) = fit_data_psycho_cum(2:8, 1);
fit_data_psycho_2(1:7, 2) = fit_data_psycho_cum(2:8, 2);
fit_data_psycho_2(1:7, 3) = fit_data_psycho_cum(2:8, 3);

% remove 4 data points (2 pair at edge)
fit_data_psycho_4(1:5, 1) = fit_data_psycho_cum(3:7, 1);
fit_data_psycho_4(1:5, 2) = fit_data_psycho_cum(3:7, 2);
fit_data_psycho_4(1:5, 3) = fit_data_psycho_cum(3:7, 3);         
         
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
        % calculate firing rate of each trial
%         for j = 1 : repetition; 
%             spike_temp = spike_rates(select);   
%             resp_heading{k}(i, j) = spike_temp(j);           
%         end
        resp_mat{k}(i) = mean(resp{k,i});  % the mean firing rate for each heading 
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
end

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %----------- plot neuronal response aross time, see whether it is stable-------------------------------------------------------------------
% figure(9);
% kk{1}='b.-'; kk{2}='r.-'; kk{3}='k.-'; kk{4}='g.-'; kk{5}='y.'; kk{6}='y.-'; 
% nu{1}='b.--'; nu{2}='r.--'; nu{3}='k.--'; nu{4}='g.--'; nu{5}='y.--'; nu{6}='y.--';  
for k = 1 : length(unique_stim_type)
%     subplot(1,3,k);
%     plot(resp{k,1},kk{1});
%     hold on;
%     plot(resp{k,5},kk{4});
%     hold on;
%     plot(resp{k,9}, nu{2});
%     hold on;
    % decide whether ves and vis is congruent tuning. Fit line by linear
	% regression first and compare the sign of each condition to decide whether
	% congruent or opposite, this is used to check whether congruent cells lead
	% to better neuronal performance in combined condition, and vice versa
    [rr,pp] = corrcoef(unique_heading, resp_mat{k}(:));
    line_re{k} = rr(1,2);
    line_p{k} = pp(1,2);
end

if length(unique_stim_type)>=2 % make sure there are ves and vis conditions
    if sign(line_re{1}) == sign(line_re{2})
        tuning_sign_vis = 0; % congruent
    else
        tuning_sign_vis = 180; % opposite
    end
    if sign(line_re{1}) == sign(line_re{3})
        tuning_sign_com = 0; % congruent
    else
        tuning_sign_com = 180; % opposite
    end
    tuning_sign_p(1) = line_p{1};
    tuning_sign_p(2) = line_p{2};
    % run correlation between ves and vis
    [rrr,ppp] = corrcoef(resp_mat{1}(:), resp_mat{2}(:));
    line2_re = rrr(1,2);
    line2_p = ppp(1,2);
else
    tuning_sign_vis = NaN;
    tuning_sign_com = NaN;
    tuning_sign_p(1) = NaN;
    tuning_sign_p(2) = NaN;
    line2_re = NaN;
    line2_p = NaN;
end    
hold off;
% %------------------------------------------------------------------------

% now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
% neurothreshold 
fit_data_neuro = [];
for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_heading)-1   % subtract the 0 heading   
        trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
        fit_data_neuro_cum{k}(i,3) = sum(trials_n);  % for later function fit use
        if i < (1+length(unique_heading))/2
             % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction            
         %    Neuro_correct{k}(i) =  rocN( resp{k,(1+length(unique_heading))/2},resp{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
        % use anti-neuron model instead, so compare each plus and minus headings 
             Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i+1},resp{k,i},100 );
        else
         %    Neuro_correct{k}(i) =  rocN( resp{k,(1+length(unique_heading))/2}, resp{k,(i+1)},100 ); % compare to the 0 heading condition, which is straght ahead
             Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i}, resp{k,(i+1)},100 );
        end
 %        if  resp_mat{k}(1) < resp_mat{k}(end)
         if line_re{k} > 0 % if left heading is not the larger responses, then linear regression will be positive 
             % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
             % here we asume if left response larger than right response then asign the left to be preferred direction
             Neuro_correct{k}(i) = 1 - Neuro_correct{k}(i);            
         end  
    end
    fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
    fit_data_neuro_cum{k}(:,2) = Neuro_correct{k}(:);
end

% remove 2 data points (1 pair at edge)
fit_data_neuro_2(1:6, 1) = fit_data_neuro_cum{1}(2:7, 1);
fit_data_neuro_2(1:6, 2) = fit_data_neuro_cum{1}(2:7, 2);
fit_data_neuro_2(1:6, 3) = fit_data_neuro_cum{1}(2:7, 3);

% remove 4 data points (2 pair at edge)
fit_data_neuro_4(1:4, 1) = fit_data_neuro_cum{1}(3:6, 1);
fit_data_neuro_4(1:4, 2) = fit_data_neuro_cum{1}(3:6, 2);
fit_data_neuro_4(1:4, 3) = fit_data_neuro_cum{1}(3:6, 3);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% use Wichman's MLE method to estimate threshold and bias
wichman_psy = pfit(fit_data_psycho_cum,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
Thresh_psy = wichman_psy.params.est(2);
psy_perf = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
xi = min(unique_heading) : 0.01 : max(unique_heading); 
yi_psy = cum_gaussfit(psy_perf, xi);   
for j = 1: length(fit_data_psycho_2(:, 1)) % catch the corresponding value
    diff = abs(xi - fit_data_psycho_2(j, 1));
    index_psy(j) = find(diff==min(diff));
end
err_psy = 100*sum( (yi_psy(index_psy)-fit_data_psycho_2(:, 2)').^2 ); % magnify 100 times

wichman_psy2 = pfit(fit_data_psycho_2,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
Thresh_psy2 = wichman_psy2.params.est(2);
psy_perf2 = [wichman_psy2.params.est(1),wichman_psy2.params.est(2)];
xi2 = min(fit_data_psycho_2(:, 1)) : 0.01 : max(fit_data_psycho_2(:, 1)); 
yi_psy2 = cum_gaussfit(psy_perf2, xi2); 
for j = 1: length(fit_data_psycho_2(:, 1))
    diff2 = abs(xi2 - fit_data_psycho_2(j, 1));
    index_psy2(j) = find(diff2==min(diff2));
end

err_psy2 = 100*sum( (yi_psy2(index_psy2)-fit_data_psycho_2(:, 2)').^2 );

wichman_psy4 = pfit(fit_data_psycho_4,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
Thresh_psy4 = wichman_psy4.params.est(2);
% psy_perf4 = [wichman_psy4.params.est(1),wichman_psy4.params.est(2)];
% xi4 = min(fit_data_psycho_4(:, 1)) : 0.01 : max(fit_data_psycho_4(:, 1)); 
% yi_psy4 = cum_gaussfit(psy_perf4, xi4); 
% for j = 1: length(fit_data_psycho_4(:, 1))
%     diff4 = abs(xi4 - fit_data_psycho_4(j, 1));
%     index_psy4(j) = find(diff4==min(diff4));
% end
% err_psy4 = 100*sum( (yi_psy4(index_psy4)-fit_data_psycho_4(:, 2)').^2 );

% neuron
wichman_neu = pfit(fit_data_neuro_cum{1},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
Thresh_neu = wichman_neu.params.est(2);
if Thresh_neu<0 | Thresh_neu>300
    Thresh_neu = 300;
end
neu_perf = [wichman_neu.params.est(1),Thresh_neu];
yi_neu = cum_gaussfit(neu_perf, xi);   
for j = 1: length(fit_data_neuro_2(:, 1)) % catch the corresponding value
    diff = abs(xi - fit_data_neuro_2(j,1));
    index_neu(j) = find(diff==min(diff));
end
err_neu = 100*sum( (yi_neu(index_neu)-fit_data_neuro_2(:, 2)').^2 ); % magnify 100 times

wichman_neu2 = pfit(fit_data_neuro_2,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
Thresh_neu2 = wichman_neu2.params.est(2);
if Thresh_neu2<0 | Thresh_neu2>300
    Thresh_neu2 = 300;
end
neu_perf2 = [wichman_neu2.params.est(1),Thresh_neu2];
yi_neu2 = cum_gaussfit(neu_perf2, xi2); 
for j = 1: length(fit_data_neuro_2(:, 1))
    diff2 = abs(xi2 - fit_data_neuro_2(j, 1));
    index_neu2(j) = find(diff2==min(diff2));
end
err_neu2 = 100*sum( (yi_neu2(index_neu2)-fit_data_neuro_2(:, 2)').^2 );

wichman_neu4 = pfit(fit_data_neuro_4,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
Thresh_neu4 = wichman_neu4.params.est(2);
if Thresh_neu4<0 | Thresh_neu4>300
    Thresh_neu4 = 300;
end
% neu_perf4 = [wichman_neu4.params.est(1),Thresh_neu4];
% yi_neu4 = cum_gaussfit(neu_perf4, xi4); 
% for j = 1: length(fit_data_neuro_4(:, 1))
%     diff4 = abs(xi4 - fit_data_neuro_4(j, 1));
%     index_neu4(j) = find(diff4==min(diff4));
% end
% err_neu4 = 100*sum( (yi_neu4(index_neu4)-fit_data_neuro_4(:, 2)').^2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure(2);
plot(fit_data_neuro_cum{1}(:,1), fit_data_neuro_cum{1}(:,2), 'o', xi, cum_gaussfit(neu_perf, xi),  'b-' );
hold on;
plot(xi2, cum_gaussfit(neu_perf2, xi2),  'r-' );

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
buff = sprintf(sprint_txt, FILE, Thresh_psy, Thresh_psy2, Thresh_psy4, Thresh_neu, Thresh_neu2, Thresh_neu4, err_psy, err_psy2, err_neu, err_neu2 );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum_VesNarrowRange.dat'];
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