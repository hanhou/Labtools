%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cum_shiftwindow_HH(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);

TEMPO_Defs;
Path_Defs;
% ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

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
    disp('Bad trial removed...');
    for i=1:length(bad_tri)
        same_tri = find( (stim_type == stim_type(bad_tri(i))) & (heading == heading(bad_tri(i))) & (motion_coherence == motion_coherence(bad_tri(i))) ); % find other trials within same condition
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
%choice(889) =2; % for cell m2c384r2 % for some reason the choice is 0 for this trial
% psychometric dataset
psycho_correct = [];
fit_data_psycho = [];
N_obs = [];

for k = 1:length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        for i = 1:length(unique_heading)
            trials_p =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) ) ;
            % make 'S' curve by using the rightward choice for y-axis
            correct_trials = (trials_p & (choice == RIGHT) );
            psycho_correct{c}(k,i) = 1*sum(correct_trials) / sum(trials_p);
            fit_data_psycho_cum{c,k}(i, 1) = unique_heading( i );
            fit_data_psycho_cum{c,k}(i, 2) = psycho_correct{c}(k,i);
            fit_data_psycho_cum{c,k}(i, 3) = sum(trials_p);
        end
    end
end

% this part needs work later, for vestibular condition, coherence does not duplicate
one_repetition = length(unique_heading)*length(unique_stim_type)*length(unique_motion_coherence);
repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition

resp_heading = [];
Z_Spikes = spike_rates;
% z-score data for later cp analysis across headings
for k = 1:length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        for i = 1:length(unique_heading)
            select = logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) ) ;
            z_dist = spike_rates(select);
            z_dist = (z_dist - mean(z_dist))/std(z_dist);
            Z_Spikes(select) = z_dist;
        end
    end
end
Z_Spikes_Ori = Z_Spikes; % keep a Z_Spikes unchanged for later use

% now group neuronal data into two groups according to monkey's choice
for k = 1:length(unique_stim_type)    % notice that the stim_type is double than disc_heading
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        for i = 1:length(unique_heading)
            select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) ) ;
            resp{c,k,i} = spike_rates(select);
            
            resp_mat{c,k}(i) = mean(resp{c,k,i});  % the mean firing rate for each heading
            resp_mat_std{c,k}(i)= std(resp{c,k,i});
            resp_mat_err{c,k}(i) = std(resp{c,k,i}) / sqrt(repetition);
            
            % calculate CP, group data based on monkey's choice
            resp_left_choose{c,k,i} = spike_rates(select & (choice == LEFT) );
            resp_right_choose{c,k,i} = spike_rates(select & (choice == RIGHT) );
            
            if (length(resp_left_choose{c,k,i}) <= 3) | (length(resp_right_choose{c,k,i}) <= 3)   % make sure each stim_type has at least 3 data values
                %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)
                Z_Spikes(select) = 9999;   % similar to NaN, just make a mark
            end
        end
        
        % now across all data (z-score for grand CP? HH)
        resp_left_all{c,k} = Z_Spikes( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) );
        resp_right_all{c,k} = Z_Spikes( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) );
        resp_all{c,k} = Z_Spikes( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
    end
end

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %----------- plot neuronal response aross time, see whether it is stable-------------------------------------------------------------------
% figure(9);
% kk{1}='b.-'; kk{2}='r.-'; kk{3}='k.-'; kk{4}='g.-'; kk{5}='y.'; kk{6}='y.-';
% nu{1}='b.--'; nu{2}='r.--'; nu{3}='k.--'; nu{4}='g.--'; nu{5}='y.--'; nu{6}='y.--';
% plot(resp{1,2,1},kk{1});
% hold on;
% plot(resp{1,2,5},kk{4});
% hold on;
% plot(resp{1,2,9}, kk{2});
% hold on;
% ylim([0 100]);
% set(gca, 'ytick',[0:20:100]);
% set(gca, 'xtick',[0 repetition]);

for k=1:length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        % decide whether ves and vis is congruent tuning. Fit line by linear
        % regression first and compare the sign of each stim_type to decide whether
        % congruent or opposite, this is used to check whether congruent cells lead
        % to better neuronal performance in combined stim_type, and vice versa
        [rr,pp] = corrcoef(unique_heading, resp_mat{c,k}(:));
        line_re{c,k} = rr(1,2);
        line_p{c,k} = pp(1,2);
    end
end

% hold off;
% %------------------------------------------------------------------------

% now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
% neurothreshold
fit_data_neuro = [];
fit_data_neuro_cut = [];
for k = 1 : length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        
        for i = 1 : sum(unique_heading ~=0)   % subtract the 0 heading
            trials_n =logical( (motion_coherence == unique_motion_coherence(c)) & (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
            fit_data_neuro_cum{c,k}(i,3) = sum(trials_n);  % for later function fit use
            fit_data_neuro_cum_anti{c,k}(i,3) = sum(trials_n);  % for later function fit use (anti neuron model)
            
            
            if i < (1+length(unique_heading))/2
                % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction
                
                % Comparing with 0 heading. HH20140510
                if sum(unique_heading ==0) ==0      % If we don't have 0 heading
                    Neuro_correct{c,k}(i) = NaN;
                else
                    Neuro_correct{c,k}(i) =  rocN( resp{c,k,unique_heading == 0},resp{c,k,i},100 ); % compare to the 0 heading stim_type, which is straght ahead
                end
                
                % Anti-neuron model, so compare each plus and minus headings
                Neuro_correct_anti{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i+1},resp{c,k,i},100 );
                
            else
                
                % Comparing with 0 heading. HH20140510
                if sum(unique_heading ==0) ==0      % If we don't have 0 heading
                    Neuro_correct{c,k}(i) = NaN;
                else
                    Neuro_correct{c,k}(i) =  rocN( resp{c,k,(1+length(unique_heading))/2}, resp{c,k,(i+1)},100 ); % compare to the 0 heading stim_type, which is straght ahead
                end
                
                % Anti-neuron model, so compare each plus and minus headings
                Neuro_correct_anti{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i+ 1-sum(unique_heading == 0)}, resp{c,k,(i+sum(unique_heading==0))},100 );
                
            end
            
            %        if  resp_mat{k}(1) < resp_mat{k}(end)
            if line_re{c,k} > 0 % if left heading is not the larger responses, then linear regression will be positive
                % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
                % here we asume if left response larger than right response then asign the left to be preferred direction
                Neuro_correct{c,k}(i) = 1 - Neuro_correct{c,k}(i);
                Neuro_correct_anti{c,k}(i) = 1- Neuro_correct_anti{c,k}(i);
            end
        end
        
        %          for i = 1 : length(unique_heading)-2   % for symetric headings/no 0 heading
        %             trials_n =logical( (motion_coherence == unique_motion_coherence(c)) & (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
        %             fit_data_neuro_cum{c,k}(i,3) = sum(trials_n);  % for later function fit use
        %             resp_comparison = [resp{c,k,length(unique_heading)/2},resp{c,k,1+length(unique_heading)/2} ];
        %             if i < length(unique_heading)/2
        %                  % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction
        %                  Neuro_correct{c,k}(i) =  rocN( resp_comparison,resp{c,k,i},100 ); % compare to the 0 heading stim_type, which is straght ahead
        %             % use anti-neuron model instead, so compare each plus and minus headings
        %            %      Neuro_correct{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i+1},resp{c,k,i},100 );
        %             else
        %                  Neuro_correct{c,k}(i) =  rocN( resp_comparison, resp{c,k,(i+2)},100 ); % compare to the 0 heading stim_type, which is straght ahead
        %         %         Neuro_correct{c,k}(i) =  rocN( resp{c,k,length(unique_heading)-i}, resp{c,k,(i+1)},100 );
        %             end
        %      %        if  resp_mat{k}(1) < resp_mat{k}(end)
        %              if line_re{c,k} > 0 % if left heading is not the larger responses, then linear regression will be positive
        %                  % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
        %                  % here we asume if left response larger than right response then asign the left to be preferred direction
        %                  Neuro_correct{c,k}(i) = 1 - Neuro_correct{c,k}(i);
        %              end
        %         end
    end
end

%choice probability
for k = 1 : length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        if unique_stim_type(k) ==1
            c=1;
        end
        
        % Each unique heading
        for i = 1 : length(unique_heading)
            if (length(resp_left_choose{c,k,i}) > 3) & (length(resp_right_choose{c,k,i}) > 3)
                CP{c,k}(i) = rocN( resp_left_choose{c,k,i},resp_right_choose{c,k,i},100 );  % raw firing rate
            else
                CP{c,k}(i) = NaN;
                
            end
            if  line_re{c,k} > 0
                CP{c,k}(i) = 1 - CP{c,k}(i);
            end
        end
        
        % Grand CP
        if (length(resp_left_all{c,k}) > 3) & (length(resp_right_all{c,k}) > 3)
            CP_all{c,k} = rocN( resp_left_all{c,k},resp_right_all{c,k},100 );   % z-scored firing rate
        else
            CP_all{c,k} = NaN;
        end
        if  line_re{c,k} > 0
            CP_all{c,k} = 1 - CP_all{c,k};
        end
        
        
    end
end

% HH20130824 Commentted
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% use Wichman's MLE method to estimate threshold and bias

% for k = 1:length(unique_stim_type)
%      for c = 1:length(unique_motion_coherence)
%        wichman_psy = pfit(fit_data_psycho_cum{c,k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
%         Thresh_psy{c,k} = wichman_psy.params.est(2);
%         Bias_psy{c,k} = wichman_psy.params.est(1);
%         psy_perf{c,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
%          fit_data_neuro_cum{c,k}(:,1) = unique_heading(unique_heading~=0);
% %        fit_data_neuro_cum{c,k}(:,1) = unique_heading( unique_heading~=unique_heading(length(unique_heading)/2) & unique_heading~=unique_heading(1+length(unique_heading)/2) );
%          fit_data_neuro_cum{c,k}(:,2) = Neuro_correct{c,k}(:);
%        wichman_neu = pfit(fit_data_neuro_cum{c,k}(:,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
%        Thresh_neu{c,k} = wichman_neu.params.est(2);
% %         [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_cum{c,k});
% %         if ttt>300
% %             ttt=300;
% %         end
% %         Thresh_neu{c,k} = ttt;
%         % negative and positive infinite value means flat tuning
%         if Thresh_neu{c,k}<0 | Thresh_neu{c,k}> 300
%             Thresh_neu{c,k} = 300;
%             wichman_neu.params.est(2) = 300;
%         end
%         Bias_neu{c,k} = wichman_neu.params.est(1);
%         neu_perf{c,k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
%     end
% end



% Uncommentted HH20130824
% % use self-written maximum liklihood method, this should run faster than
% above
for c = 1:length(unique_motion_coherence)
    for k = 1:length(unique_stim_type)
        [bb,tt] = cum_gaussfit_max1(fit_data_psycho_cum{c,k});
        Bias_psy{c,k} = bb;
        Thresh_psy{c,k} = tt;
        psy_perf{c,k}=[bb, tt];
        
        fit_data_neuro_cum{c,k}(:,1) = unique_heading(unique_heading~=0);
        fit_data_neuro_cum{c,k}(:,2) = Neuro_correct{c,k}(:);
        
        fit_data_neuro_cum_anti{c,k}(:,1) = unique_heading(unique_heading~=0);
        fit_data_neuro_cum_anti{c,k}(:,2) = Neuro_correct_anti{c,k}(:);
        
        if sum(isnan(fit_data_neuro_cum{c,k})) == 0 % OK
            
            [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_cum{c,k});
            Bias_neu{c,k} = bbb;
            Thresh_neu{c,k} = ttt;
            neu_perf{c,k}=[bbb, ttt];
        else
            Bias_neu{c,k} = NaN;
            Thresh_neu{c,k} = NaN;
            neu_perf{c,k}=[NaN, NaN];
        end
        
        [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_cum_anti{c,k});
        Bias_neu_anti{c,k} = bbb;
        Thresh_neu_anti{c,k} = ttt;
        neu_perf_anti{c,k}=[bbb, ttt];
        
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%do permutation to test the significance of CP_all{k}, re-calculate CP 2000 times
perm_num = 1000;
Z_Spikes_perm = Z_Spikes;
bin = 0.005;
x_bin = 0 : bin : 1;



if isempty(batch_flag) ;progressbar('Stim type','Perm'); end

for k = 1:length(unique_stim_type)
    for c = 1:length(unique_motion_coherence)
        
    % Automatic change the "dead-ahead"
    if abs(Bias_psy{c,k}) > 0.5
        %      [~,i_0]=min(abs(Bias_psy{c,k}-unique_heading)); % HH20130905
        [~,i_0]=min(abs(psycho_correct{c}(k,:)-0.5)); % HH20140510
    else
        i_0 = find(unique_heading >= 0,1);
    end        

        if unique_stim_type(k) ==1
            c=1;
        end
        for n = 1 : perm_num
            
            % temperarilly only use near-threshold heading angles where monkey make a guess mainly
            select = logical( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
            Z_Spikes_con{c,k} = Z_Spikes_perm( select );
            Z_Spikes_con{c,k} = Z_Spikes_con{c,k}(randperm(length(Z_Spikes_con{c,k})));   % permute spike_rates
            Z_Spikes_perm(select) = Z_Spikes_con{c,k};    % now in spike_rates, the corresponding data were permuted already
            
            resp_left_all_perm{c,k} = Z_Spikes_perm( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) );
            resp_right_all_perm{c,k} = Z_Spikes_perm( (motion_coherence == unique_motion_coherence(c)) & (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) );
            
            if  (length(resp_left_all{c,k}) > 3) && (length(resp_right_all{c,k}) > 3)
                CP_all_perm{c,k}(n) = rocN( resp_left_all_perm{c,k}, resp_right_all_perm{c,k},100 );
            else
                CP_all_perm{c,k}(n) = NaN;
            end
            if  line_re{c,k} > 0
                CP_all_perm{c,k}(n) = 1 - CP_all_perm{c,k}(n);
            end
            
            resp_left_choose_perm = Z_Spikes_perm((heading == unique_heading(i_0)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) & (choice == LEFT) & (Z_Spikes~=9999) );
            resp_right_choose_perm = Z_Spikes_perm((heading == unique_heading(i_0)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c)) & (choice == RIGHT) & (Z_Spikes~=9999) );
            
            if  (length(resp_left_choose{c,k,i_0}) > 3) && (length(resp_right_choose{c,k,i_0}) > 3)
              try
                  CP_perm(n) = rocN( resp_left_choose_perm, resp_right_choose_perm,100 );
              catch
                  keyboard
              end
                  else
                CP_perm(n) = NaN;
            end
            if  line_re{c,k} > 0
                CP_perm(n) = 1 - CP_perm(n);
            end
            
            if isempty(batch_flag); progressbar([],[n/perm_num]);end
            
        end
        % now calculate p value or significant test
        if (length(resp_left_all{c,k}) > 3) & (length(resp_right_all{c,k}) > 3)
            hist_perm = [];
            hist_perm = hist( CP_all_perm{c,k}(:), x_bin );  % for permutation
            bin_sum = 0;
            n = 0;
            while ( n < (CP_all{c,k}/bin) )
                n = n+1;
                bin_sum = bin_sum + hist_perm(n);
                if CP_all{c,k} >= 0.5                  % note it's two tail test
                    p{c,k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
                else
                    p{c,k} = 2* bin_sum / perm_num;
                end
            end
        else
            p{c,k} = NaN;
        end
        
        % calculate p value for CP during straight ahead motion
        
        if (length(resp_left_choose{c,k,i_0}) > 3) && (length(resp_right_choose{c,k,i_0}) > 3)
            hist_perm = [];
            hist_perm = hist( CP_perm(:), x_bin );  % for permutation
            bin_sum = 0;
            n = 0;
           
            while ( n <= (CP{c,k}(i_0)/bin) )
                n = n+1;
                bin_sum = bin_sum + hist_perm(n);
                if CP{c,k}(i_0) > 0.5    % note it's two tail test
                    pp = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
                else
                    pp = 2* bin_sum / perm_num;
                end
            end
        else
            pp = NaN;
        end
        p_a(c,k) = pp;
    end
    if isempty(batch_flag) ;progressbar(k/length(unique_stim_type),[]);end
    
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot psychometric and neurometric function here
% modified by HH20130901

h{1} = 'ko';  f{1} = 'k-';  g{1} = 'ko-';
h{2} = 'rd';  f{2} = 'r-';  g{2} = 'rd-';
h{3} = 'gs';  f{3} = 'g-';  g{3} = 'gs-';

figure(11); clf;

if ~isempty(batch_flag); set(11,'visible','off'); end

set(11,'Position', [20,40, 750,750], 'Name', 'Psycho-neurometic function & CP', 'color','w');
orient landscape;


% plot psychometric and neurometric function
% axes('position',[0.32 0.1, 0.25 0.6],'xtickmode','auto');

gaps = 0.1;
margins = [0.18 0.1 0.1 0.1];
subplot_tight(2,2,2,gaps,margins);

legend_txt = [];
for k = 1:length(unique_stim_type)
    xi = min(unique_heading) : 0.1 : max(unique_heading);
    beta = [0, 1.0];
    
    % Psychometric
    plot(unique_heading, psycho_correct{1}(k,:), h{unique_stim_type(k)},'markerfacecolor',h{unique_stim_type(k)}(1),'markersize',10);
    hold on;
    plot(xi, cum_gaussfit(psy_perf{1,k}, xi),  f{unique_stim_type(k)},'linewidth',2.5);
    
    % Neurometric
    neu_heading = unique_heading(unique_heading~=0);
    plot(neu_heading, Neuro_correct{1,k}, '^','linewidth',2,'markersize',10,'color',h{unique_stim_type(k)}(1));
    plot(xi, cum_gaussfit(neu_perf{1,k}, xi),  '--','linewidth',2.5,'color',h{unique_stim_type(k)}(1));
    plot(xi, cum_gaussfit(neu_perf_anti{1,k}, xi),  ':','linewidth',1.5,'color',h{unique_stim_type(k)}(1));
    
    title('Psychometric & Neurometric');
    set(gca,'xtickmode','auto');
    xlabel('Heading Angle (deg)');
    ylim([0,1]);
    ylabel('Rightward Choices');
    legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
    legend_txt{k*2} = [''];
    
    % annotation
    text(max(unique_heading)/6,0.25-k*0.07,sprintf('%6.3g %6.3g(%6.3g)',psy_perf{1,k}(2),neu_perf{1,k}(2),neu_perf_anti{1,k}(2)),'color',h{unique_stim_type(k)}(1));
end

text(max(unique_heading)/6,0.25,sprintf('    \\itpsy  neu'),'color','b');

xlim([min(unique_heading)*1.1 max(unique_heading)*1.1]);



%%%%%%  neurological raw data based on firing rate instead of ROC
% axes('position',[0.03 0.1, 0.25 0.6]);
subplot_tight(2,2,1,gaps,margins);

for k = 1:length(unique_stim_type)
    set(errorbar(unique_heading, resp_mat{1,k}(:), resp_mat_err{1,k}(:),g{unique_stim_type(k)} ),'linewidth',2);
    xlabel('Heading Angle (deg)');
    ylabel('Firing rate(spikes/s)');
    xlim([min(unique_heading)*1.1,max(unique_heading)]*1.1);
    hold on;
    text(min(unique_heading),resp_mat{1,k}(1)*0.9,sprintf('%4.2g',line_re{k}),'color',h{unique_stim_type(k)}(1));
end
title(['Tuning curve']);
set(gca,'xtickmode','auto');


%{

%% CP distribution for all trials (grand CP). HH20130901
axes('position',[0.7, 0.52, 0.25, 0.2]);

% Temporarily for only one condition
c=1;
k=1;

if ~isempty(resp_all{c,k})
    xbins = linspace(min(resp_all{c,k}(:)),max(resp_all{c,k}(:)),15);
    if line_re{1} < 0  % Prefers left
        prefN = hist(resp_left_all{c,k},xbins);
        nullN = hist(resp_right_all{c,k},xbins);
    else
        prefN = hist(resp_right_all{c,k},xbins);
        nullN = hist(resp_left_all{c,k},xbins);
    end
    
    hold on;
    h1_ = bar(xbins,[nullN' prefN'],1.5,'grouped','faceColor','k');
    set(h1_(1),'facecolor','w');
    xlabel('Z-scored firing rate');
    ylabel('Frequency');
    legend('null choice', 'pref choice');
    title(sprintf('Grand \\itCP = \\rm%6.3g, \\itp = \\rm%4.3g, rep = %g',CP_all{c,k},p{1,k},repetition));
    set(findall(gcf,'FontSize',10),'FontSize',13);
end

%% CP for 0 heading. HH20130901
axes('position',[0.7, 0.18, 0.25, 0.2]);

% Temporarily for only one condition
c=1;
k=1;

% i_0 = find(unique_heading == 0);

xbins = linspace(min([resp_left_choose{c,k,i_0}(:);resp_right_choose{c,k,i_0}(:)]),max([resp_left_choose{c,k,i_0}(:);resp_right_choose{c,k,i_0}(:)]),15);
if line_re{1} < 0  % Prefers left
    
    prefN = hist(resp_left_choose{c,k,i_0},xbins);
    nullN = hist(resp_right_choose{c,k,i_0},xbins);
else
    prefN = hist(resp_right_choose{c,k,i_0},xbins);
    nullN = hist(resp_left_choose{c,k,i_0},xbins);
end

hold on;
h1_ = bar(xbins,[nullN' prefN'],1.5,'grouped','faceColor','k');
set(h1_(1),'facecolor','w');
xlabel('Firing rate (Hz)');
ylabel('Frequency');
xlims = xlim();
title(sprintf('\\itCP \\rmat %g\\circ = \\rm%6.3g, \\itp = \\rm%4.3g, rep = %g ', unique_heading(i_0), CP{c,k}(i_0), p_a(1,k),repetition));
set(findall(gcf,'FontSize',10),'FontSize',13);

%}

%% Shift Window (from HeadingDis_cum_shiftwindow_HH.   HH20140601)

spike_data = squeeze(data.spike_data(SpikeChan,:,select_trials))';   % spike rasters  % Trial Nums * 5000
spike_data( spike_data>100 ) = 1; % something is absolutely wrong

% shift window, currently using 100ms shift with 1000ms width
window_interval = 1000;
step=100;
endloop = floor((3000-window_interval/2)/100);

if isempty(batch_flag) ;progressbar('ShiftWindows'); end

clear spike_rates;


for s = 1 : endloop
   
    % Each trial in this shifting window
    spike_rates = sum(spike_data(:,StartEventBin(1)-0.5*window_interval+step*(s-1): StartEventBin(1)+ 0.5*window_interval+step*(s-1) ),2)';
    
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
    %     if FILE=='m2c384r2.htb'
    %        choice(889) =2; % for cell m2c384r2 % for some reason the choice is 0 for
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    % now group neuronal data into two groups according to monkey's choice
    for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading
        for i = 1:length(unique_heading)
            select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
            resp{k,i} = spike_rates(select);
            resp_mat{k}(i) = mean(resp{k,i});  % the mean firing rate for each heading
            resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);            
             resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );
             resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );
            if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
          %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)  
                Z_Spikes(select) = 9999;   % similar to NaN, just make a mark            
           %     Z_Spikes( (heading == unique_heading(length(unique_heading)+1-i)) & (stim_type == unique_stim_type(k)) ) = 9999; % the corresponding heading is also excluded
          %      CP{k}(length(unique_heading)+1-i) = 9999;
            else
            end 

        end
        % now across all data
        resp_left_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) );
        resp_right_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) );
        resp_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
        
    end
    
    % Note that here we fit a local tuning curve for each shifting window, so the "preferred" heading may change its side
    % during the shifting. Therefore, it would be neccessary to mark the preferred heading across time.  HH20140602
    
    for k = 1 : length(unique_stim_type)
        % decide whether ves and vis is congruent tuning. Fit line by linear
        % regression first and compare the sign of each condition to decide whether
        % congruent or opposite, this is used to check whether congruent cells lead
        % to better neuronal performance in combined condition, and vice versa
        [rr,pp] = corrcoef(unique_heading, resp_mat{k}(:));
        line_re_shift{k}(s) = rr(1,2);  % Resotre the preferred directions for each window
        line_p{k} = pp(1,2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
    % neurothreshold
    fit_data_neuro = [];
    fit_data_neuro_cut = [];
    for k = 1 : length(unique_stim_type)
        fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
        for i = 1 : sum(unique_heading ~=0)   % subtract the 0 heading, if we have
            trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
                        
            if i < (1+length(unique_heading))/2
                % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction
                Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i+1},resp{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
            else
                Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i+ 1-sum(unique_heading == 0)}, resp{k,(i+ sum(unique_heading == 0))},100 ); % compare to the 0 heading condition, which is straght ahead
            end
            
            %        if  resp_mat{k}(1) < resp_mat{k}(end)
            if line_re_shift{k}(s) > 0 % if left heading is not the larger responses, then linear regression will be positive
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
            if (length(resp_left_all{k}) > 3) | (length(resp_right_all{k}) > 3)
                %if  (length(resp_left_all{k}) / length(resp_all{k}) >0.25) |  (length(resp_left_all{k}) / length(resp_all{k}) <0.75)
                CP_all_thisbin{k} = rocN( resp_left_all{k},resp_right_all{k},100 );
            else
                CP_all_thisbin{k} = NaN;
            end
            if  line_re_shift{k}(s) > 0
                CP_all_thisbin{k} = 1 - CP_all_thisbin{k};
            end
end
    % %--------------------------------------------------------------------------
    %do permutation to test the significance of CP_all{k}, re-calculate CP 1000 times
    %     perm_num = 1000;
    %     Z_Spikes_perm = Z_Spikes;
    %     bin = 0.005;
    %     x_bin = 0 : bin : 1;
    %     Z_Spike_con = cell(1,perm_num);
    %
    %     for k = 1:1
    %  %   for k = 1:length(unique_stim_type)
    %         for n = 1 : perm_num
    %             % temperarilly only use near-threshold heading angles where monkey make a guess mainly
    %             select = logical( (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
    %             Z_Spikes_con{k} = Z_Spikes_perm( select );
    %             Z_Spikes_con{k} = Z_Spikes_con{k}(randperm(length(Z_Spikes_con{k})));   % permute spike_rates
    %             Z_Spikes_perm(select) = Z_Spikes_con{k};    % now in spike_rates, the corresponding data were permuted already
    %
    %             resp_left_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) );
    %             resp_right_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) );
    %
    %             if  (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3)
    %                 CP_all_perm{k}(n) = rocN( resp_left_all_perm{k}, resp_right_all_perm{k},100 );
    %             else
    %                 CP_all_perm{k}(n) = NaN;
    %             end
    %             if  line_re_shift{k}(s) > 0
    %                 CP_all_perm{k}(n) = 1 - CP_all_perm{k}(n);
    %             end
    %         end
    %         % now calculate p value or significant test
    %         if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3)
    %
    %             hist_perm(k,:) = hist( CP_all_perm{k}(:), x_bin );  % for permutation
    %             bin_sum = 0;
    %             n = 0;
    %             while ( n < (CP_all{k}/bin) )
    %                  n = n+1;
    %                  bin_sum = bin_sum + hist_perm(k, n);
    %                  if CP_all{k} >= 0.5                  % note it's two tail test
    %                     p{k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
    %                  else
    %                     p{k} = 2* bin_sum / perm_num;
    %                  end
    %             end
    %         else
    %             p{k} = NaN;
    %         end
    %     end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% use Wichman's MLE method to estimate threshold and bias
    
    for k = 1:length(unique_stim_type)
        
        [bb,tt] = cum_gaussfit_max1(fit_data_neuro_cum{k});
        Thresh_neu_thisbin{k} = tt;
        % negative and positive infinite value means flat tuning
        if Thresh_neu_thisbin{k}<0 | Thresh_neu_thisbin{k}> 300
            Thresh_neu_thisbin{k} = 300;
        end
    end
    for k = 1:length(unique_stim_type)
        % calculate cp and threshold based on different analyze window
        CP_all_shift(k,s) = CP_all_thisbin{k};
        %         p_shift(1,s) = p{1};
        [~,p_shift(k,s)] = ttest2(resp_left_all{k},resp_right_all{k}); % Use ttest2 instead of permutation for speed (the outcomes are quite similar). HH20140601
        Thresh_neu_shift(k,s) = Thresh_neu_thisbin{k};
    end
    
    if isempty(batch_flag)  ;progressbar(s/endloop); end
end


% HH20140511 plot CP time course

stim_type_names = {'Vest','Vis','Comb'}; % stim_type = 0, 1, 2, 3
colors = {'k','r','g'};
xx = 0:step:step*(endloop-1);

% set(figure(777),'color','w','position',[5 50 1400 700]);
% if ~isempty(batch_flag); set(777,'visible','off'); end
ax = subplot_tight(2,2,[3 4],gaps,margins);

clear h;
ylimCP = [0.3 1.06];

for k = 1:length(unique_stim_type)
    % Plot neurometric and CP
    if k==1
        [ax,h(k,1),h(k,2)] = plotyy(ax,xx,CP_all_shift(k,:),xx,Thresh_neu_shift(k,:));
    else
        axes(ax(1)); hold on; h(k,1) = plot(xx,CP_all_shift(k,:));
        axes(ax(2)); hold on; h(k,2) = plot(xx,Thresh_neu_shift(k,:));
    end
    
    set(h(k,1),'linewidth',1,'marker','o','markersize',10,'color',colors{unique_stim_type(k)})
    set(h(k,2),'linewidth',1,'linestyle','--','marker','^','markersize',10,'color',colors{unique_stim_type(k)})
    
    axes(ax(1)); hold on;
    plot(xx(p_shift(k,:)<0.05),CP_all_shift(k,p_shift(k,:)<0.05),['o' colors{unique_stim_type(k)}],'markerfacecolor',colors{unique_stim_type(k)},'markersize',10);   % Sig. CP

    % Plot the "preferred direction" for each window. HH20140602
    plot(xx,ylimCP(1)+(line_re_shift{k}>0)*(range(ylimCP)-0.06)+(3-k)*0.02,['s' colors{unique_stim_type(k)}],'markerfacecolor',colors{unique_stim_type(k)},'markersize',5);

    % Plot Psychometric threshold
    axes(ax(2)); hold on; plot([xx(1) xx(end)],[Thresh_psy{1,k} Thresh_psy{1,k}],['-' colors{unique_stim_type(k)}]);
end

set(ax,'xlim',[0 xx(end)]);

axes(ax(1)); hold on;
lims = axis; hold on;
plot([lims(1) lims(2)],[0.5 0.5],'k--');
ylim(ylimCP);
set(ax(1),'xtick',[],'ytick',[ylimCP(1):0.1:ylimCP(2)],'ycolor','k');
ylabel(ax(1),sprintf('Grand CP (O)'))

visDur = round(mean(StopEventBin-StartEventBin)/100)*100;
plot(ax(1),[visDur visDur],ylimCP,'k');

axes(ax(2));
ylim([1 200]);
xlabel(ax(2),['Center of ' num2str(window_interval) ' ms time window (ms)']);
ylabel(ax(2),sprintf('Neuronal Threshold (\\Delta anti-model)'))
set(ax(2),'yscale','log','yticklabel',[1 10 100],'ytick',[1 10 100],'ycolor','k');

% title([FILE '_' num2str(SpikeChan)]);
SetFigure(10); 

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output some text of basic parameters in the figure
axes('position',[0.05,0.88, 0.9,0.08] );

xlim( [0,100] );
ylim( [2,10] );
text(0,10, [FILE 'unit' num2str(SpikeChan)]);
text(20,10,sprintf('coherence = %g   rep = %g',unique_motion_coherence,repetition)  );
text(0,8, 'stim  Psy: u   threshold    Neu:u     threshold(anti-model)    CP\_all      p      CP\_0     p');

for k = 1:length(unique_stim_type)
    text(0,8-k*2,sprintf('%g    %5.3f    %5.3f    %5.3f    %5.3f(%5.3f)    %5.3f(p= %5.3f)    %5.3f(p= %5.3f)', unique_stim_type(k),...
        Bias_psy{1,k},Thresh_psy{1,k},Bias_neu{1,k},Thresh_neu{1,k},Thresh_neu_anti{1,k},...
                     CP_all{1,k}',p{1,k}',CP{c,k}(i_0),p_a(1,k)));
end
axis off;
% tightfig;

% Output to command line for Excel
fprintf(' rep,	Psy_u,	Psy_thres,	Neu_u,	Neu_thres, Neu_thres(anti-model),CP_all,	p_CP_all,	CP_0,	p_CP_0\n');
fprintf('%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n', repetition, Bias_psy{1,k},Thresh_psy{1,k},Bias_neu{1,k},Thresh_neu{1,k},Thresh_neu_anti{1,k},CP_all{1,k}',p{1,k}',CP{c,k}(i_0),p_a(1,k));


%%%%%%%%%%%%%%%%%%%%%%%  Batch Output.   HH20140510  %%%%%%%%%%%%%%%%%

if ~isempty(batch_flag)
    
    outpath = ['Z:\Data\Tempo\Batch\' batch_flag(1:end-2) '\'];
   
    if ~exist(outpath,'dir')
        mkdir(outpath);    
    end
    
    % Save figures
    orient landscape;
    savefilename = [outpath [FILE '_' num2str(SpikeChan)] '_ShiftCP.png'];
    if exist(savefilename)
        delete(savefilename);
    end
    saveas(11,savefilename,'png');
        
    % Print results
    sprint_txt_temp = 'gggggggggssss';
    sprint_txt = [];
    for i = 1:length(sprint_txt_temp)
        sprint_txt = [sprint_txt '%' sprint_txt_temp(i) '\t '];
    end
    
    outfile = [outpath 'ShiftCP.dat'];
    printHead = 0;
    if (exist(outfile, 'file') == 0)   % file does not yet exist
        printHead = 1;
    end
    
    fid = fopen(outfile, 'a');
    if (printHead)
        fprintf(fid, 'FILE\t   rep,	Psy_u,	Psy_thres,	Neu_u,	Neu_thres, Neu_thres, (anti-model),CP_all,	p_CP_all,	CP_0,	p_CP_0,  xx, Thresh_neu_shift(k,:), CP_all_shift(k,:), p_shift(k,:)  ');
        fprintf(fid, '\r\n');
    end
    
    fprintf(fid,'%s\t %g\t ',[FILE '_' num2str(SpikeChan)], repetition);

    for conditions = 1:3 % Always output 3 conditions (if not exist, fill with NaNs)
        if sum(unique_stim_type == conditions)==0
            buff = sprintf(sprint_txt, ones(1,length(sprint_txt_temp))*NaN);
        else
            k = find(unique_stim_type == conditions);
            buff = sprintf(sprint_txt,  Bias_psy{1,k},Thresh_psy{1,k},Bias_neu{1,k},Thresh_neu{1,k},Thresh_neu_anti{1,k},...
                CP_all{1,k}',p{1,k}',CP{c,k}(i_0),p_a(1,k), ...
                num2str(xx), num2str(Thresh_neu_shift(k,:)), num2str(CP_all_shift(k,:)), num2str(p_shift(k,:)));
        end
        fprintf(fid, '%s', buff);
    end
    
    fprintf(fid, '\r\n');
    fclose(fid);
    
end;

%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
%{
sprint_txt = ['%s'];
for i = 1 : 100 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];
end
%buff = sprintf(sprint_txt, FILE, unique_motion_coherence, repetition, Bias_psy{:}, Thresh_psy{:}, Bias_neu{:}, Thresh_neu{:}, CP_all{:}, p{:}, line_re{:}, line_p{:}, line2_re, line2_p );
if length(unique_stim_type)==1
   buff = sprintf(sprint_txt, FILE, 2, Thresh_neu{length(unique_motion_coherence),1} );   % visual 100% coherence
elseif length(unique_stim_type)>1
    buff = sprintf(sprint_txt, FILE, 1,CP{1,1}(5),CP{1,2}(5),CP{1,3}(5) );   % vestibular
end
%buff = sprintf(sprint_txt, resp_mat{1}, resp_mat_std{1}.^2, resp_mat{2}, resp_mat_std{2}.^2, resp_mat{3}, resp_mat_std{3}.^2  );
%buff = sprintf(sprint_txt, FILE, Thresh_neu{:}, CP_all{:}, p{:}, line_re{:} );
%buff = sprintf(sprint_txt, FILE, Thresh_neu{:},line_re{1}, line_re{2}, line_p{1}, line_p{2} );
%buff = sprintf(sprint_txt, FILE, resp_mat{1}(:),resp_mat{2}(:),resp_mat{3}(:) );
%buff = sprintf(sprint_txt, FILE, fit_data_psycho_cum{1}(:, 2),fit_data_psycho_cum{2}(:, 2),fit_data_psycho_cum{3}(:, 2),Neuro_correct{1}(:),Neuro_correct{2}(:),Neuro_correct{3}(:) );
outfile = ['Z:\Users\Yong\CP.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cum.dat'];
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

% figure(10);saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '.png'],'png')
% saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '.fig'])
% close (9);close (10);
% [StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
% HeadingDis_cum_PSTH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
%
% figure(2);orient landscape
% saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '_PSTH.png'],'png')
% saveas(gcf,['C:\Aihua\Chaos\' FILE(1:end-4) '_PSTH.fig'])
% close(2);
%}

if ~isempty(batch_flag)
    % Do this together
    HeadingDis_cum_PSTH_HH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);
end

return;