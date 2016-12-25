% % isolate more than 2 units from single unit recording data by offline
% spikesorting, analyze clustering structure and noise correlation among
% units --YG, 03/08
% %-----------------------------------------------------------------------------------------------------------------------
function HeadingDis_cum_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_conflict_angle = data.moog_params(CONFLICT_ANGLE,:,MOOG); %for chris' 2I conflict study data, unecessary if not analyzing his files
temp_spike_data = data.spike_data(:,:);

%remove trials that do not fall between BegTrial and EndTrial
%also parse out null trials. use for measuring spontaneous activity later
trials = 1:length(temp_heading);
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to remove trials from the middle of a block
% take_out = 514:756;
% select_trials(take_out)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

null_trials = logical(temp_heading == data.one_time_params(NULL_VALUE));
stim_type = temp_stim_type( ~null_trials & select_trials );
heading = temp_heading( ~null_trials & select_trials );
coherence = temp_motion_coherence( ~null_trials & select_trials );
conflict_angle = temp_conflict_angle( ~null_trials & select_trials );

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_coherence = munique(coherence');
% unique_conflict_angle = munique(conflict_angle');
% disc_heading = unique_heading( floor(length(unique_heading)/2)+1 : end ); %non mirrored headings

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_heading);
Discard_trials = find(null_trials==1 | ~select_trials);
for i = 1 : length(Discard_trials)
    temp_spike_data( : , ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 99;
end

repetition = floor( length(heading)/(length(unique_heading)*length(unique_coherence)) ); %to find the number of complete repetitions
% ^ABOVE repetition WILL NOT WORK IF IMBALANCED NUMBER OF TRIALS LIKE WITH VESTIBULAR ONLY -SF 7.09
% s = length(unique_stim_type);
% h = length(unique_heading);
% c = length(unique_coherence);
% a = length(unique_conflict_angle);
% repetition = floor( length(heading) / (h+ h*c + h*c*a));
% ^ ABOVE repetition FOR CHRIS'S CONFLICT STUDY ONLY.  COMMENT IT OUT IF NOT IN USE-SF 7.09

% extract channel information
% channelnum_temp = size(data.spike_data);
% channelnum = channelnum_temp(1,1);
% channelcount = 0;
% channel_analyze = [];
% for c=1:channelnum
% %   if (sum(data.spike_data(c,:))>100) && c~=2 % dont use sympulse 2nd channel.  
%     if ((sum(data.spike_data(c,:))>100) && c~=2) % dont use sympulse 2nd channel nor original 1st channel if spike sorted.  
%         channelcount = channelcount+1;
%         channel_analyze(channelcount) = c; %get indeces of channels to analyze
%     end
% end
channelnum_temp = size(temp_spike_data);
channelnum = channelnum_temp(1,1); % how many channels
channelcount = 0;
SpikeChan = 1; % define the first channel you want to start here
for c = 1: channelnum
    temp(1,:) = temp_spike_data( c, find(temp_spike_data(1,:)~=99) );
    spikesum(c) = sum(temp(1,:));
    if c>=SpikeChan & spikesum(c)>20 & c~=2 % all those channels later than the first channel, but exclude the second synpulse channel
        channelcount = channelcount+1;
        channel_analyze(channelcount) = c; % the final channels that need to analyze
    end    
end
% channelcount=3;
% channel_analyze=[4 5 6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%begin noise correlation analysis
StartEventBin(1)=996; %determines when during trial to begin analysis
offset = 615;
duration = 1000;
for c=1:channelcount
    spike_data(c,:)=temp_spike_data( channel_analyze(c), (temp_spike_data(channel_analyze(c),:)~=9999) );
    spike_data(c , (spike_data(c,:)>10) )=1; % something is absolutely wrong
    spike_rates = zeros(1,length(heading));
    for ss=1:length(heading)%ss is index of the trials
        spike_rates(ss) = sum( spike_data(c,StartEventBin(1)+offset+5000*(ss-1) : StartEventBin(1)+offset+duration+5000*(ss-1)) ) ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this part deals with for whatever reason, some trials are absolutely
    % wrong, like >1000 spikes/s, makes no sense. tick this bad trial out,
    % replace with one of the other data within same condition, this is similar
    % to bootstrap. But in order to keep every running the same, use the mean
    % of the rest instead. this happen only very rarely!!!
    bad_tri = find(spike_rates > 1000 );
    if ~isempty(bad_tri)
        for i=1:length(bad_tri)
            same_tri = find( (stim_type == stim_type(bad_tri(i))) & (heading == heading(bad_tri(i))) ); % find other trials within same condition
            rest_tri = same_tri(same_tri ~= bad_tri(i));
            spike_rates(bad_tri(i)) = mean( spike_rates(rest_tri) );
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spike_rates_channel{c} = spike_rates;
    spike_rates_z = spike_rates; % no high-pass filter

    % create z spikes
%     for h=1:length(unique_stim_type)
%         if unique_stim_type(h)~=1
            for k=1:length(unique_coherence)
                for i=1:length(unique_heading)
%                     select = logical( (coherence==unique_coherence(k)) & (heading==unique_heading(i)) & stim_type==unique_stim_type(h)); %select for Chris conflict study data
                    select = logical( (coherence==unique_coherence(k)) & (heading==unique_heading(i)) ); %select for coherence study
                    % z-score data for spike count correlation analysis
                    z_dist = spike_rates_z(select);
                    warning off;
                    if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                        z_dist = (z_dist - mean(z_dist))/std(z_dist);
                    else
                        z_dist = 0;
                    end
                    warning on;
                    Z_Spikes(select) = z_dist;
                end
            end
%         end
%     end
    Z_Spikes_channel{c} = Z_Spikes; % contains z-scores of each and every selected trial sorted by channel
end
clear temp_spike_data;



% now analyze noise correlation between pairs
% first compute all possible pairs
channelcount_temp = channelcount-1;
num_pairs = 0;
while channelcount_temp>=1
    num_pairs = num_pairs + channelcount_temp;
    channelcount_temp = channelcount_temp - 1;
end

pair_member1 = 1;
pair_member2 = pair_member1+1;
for i = 1:num_pairs % all possible pairs
    % remove slow fluctuations at every 20 trials
    zdivisions = floor(length(spike_rates)/20);
    zz1 = Z_Spikes_channel{pair_member1};
    zz2 = Z_Spikes_channel{pair_member2};
    % ^ FOR WHEN INCLUDING VESTIBULAR CONDITION, NOT IN COHERENCE STUDY
%     zdivisions = floor(length(spike_rates(stim_type~=1 & conflict_angle==0))/20);
%     zz1 = Z_Spikes_channel{pair_member1}(stim_type~=1 & conflict_angle==0);
%     zz2 = Z_Spikes_channel{pair_member2}(stim_type~=1 & conflict_angle==0);
    % ^ FOR WHEN EXCLUDING VESTIBULAR CONDITION, AS IN COHERENCE STUDY
    
    z1all=[];
    z2all=[];
    for zz=1:zdivisions
        if zz<zdivisions
            z1=( zz1(1+(zz-1)*20:20+(zz-1)*20) - mean(zz1(1+(zz-1)*20:20+(zz-1)*20)) ) / std( zz1(1+(zz-1)*20:20 + (zz-1)*20) );
            z2=( zz2(1+(zz-1)*20:20+(zz-1)*20) - mean(zz2(1+(zz-1)*20:20+(zz-1)*20)) ) / std( zz2(1+(zz-1)*20:20 + (zz-1)*20) );
        else
            z1=( zz1(1+(zz-1)*20:end) - mean(zz1(1+(zz-1)*20:end)) ) / std( zz1(1+(zz-1)*20:end) );
            z2=( zz2(1+(zz-1)*20:end) - mean(zz2(1+(zz-1)*20:end)) ) / std( zz2(1+(zz-1)*20:end) );
        end
        z1all=[z1all z1];
        z2all=[z2all z2];
        z1all(z1all>3)=3;  % cutoff between -3 and 3
        z1all(z1all<-3)=-3;  % cutoff between -3 and 3
        z2all(z2all<-3)=-3;  % cutoff between -3 and 3
        z2all(z2all>3)=3;  % cutoff between -3 and 3
    end
    % noise correlation with all stimuli conditions included
%     [rr,pp] = corrcoef(z1all,z2all);
%     noise_r(i) = rr(1,2);
%     noise_p(i) = pp(1,2);

    % separated between stimuli conditions
%     stim_type_temp = stim_type((stim_type~=1 & conflict_angle==0));%for excluding vestibular condition as in Chris conflict study
%     for k=1:length(unique_stim_type) % ananlyze noise correlation in different conditions, if find no difference, combine later
%         select_stim = logical( stim_type_temp==unique_stim_type(k) );
%         % noise correlation with stim type separated
%         Z_Spikes1 = z1all(select_stim);
%         Z_Spikes2 = z2all(select_stim);
% 
%         [rr,pp] = corrcoef(Z_Spikes1,Z_Spikes2);
%         noise_r_stim(i,k) = rr(1,2);
%         noise_p_stim(i,k) = pp(1,2);
%     end

    % separated between coherence conditions
%     coherence_temp = coherence((stim_type~=1 & conflict_angle==0));%for excluding vestibular condition as in Chris conflict study
    for k=1:length(unique_coherence) % ananlyze noise correlation in different conditions, if find no difference, combine later
%         select = logical( coherence_temp==unique_coherence(k));
        select = logical( coherence==unique_coherence(k));
        % noise correlation with coherence separated
        Z_Spikes1 = z1all(select);
        Z_Spikes2 = z2all(select);

        [rr,pp] = corrcoef(Z_Spikes1,Z_Spikes2);
        noise_r_coh(i,k) = rr(1,2);
        noise_p_coh(i,k) = pp(1,2);
    end

    if pair_member2 < num_pairs % more to run
        pair_member2 = pair_member2+1;
    else
        pair_member1 = pair_member1+1;
        pair_member2 = pair_member1+1;
    end
end

%end noise correlation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin CP analysis
LEFT = 1;
RIGHT = 2;
choice = zeros(1,length(spike_rates));
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
    for k = 1:length(unique_coherence)
         trials_p =logical( (heading == unique_heading(i)) & (coherence == unique_coherence(k)) ) ;
         % make 'S' curve by using the rightward choice for y-axis
         correct_trials = (trials_p & (choice == RIGHT) );
         psycho_correct(k,i) = 1*sum(correct_trials) / sum(trials_p); 
         fit_data_psycho_cum{k}(i, 1) = unique_heading( i );  
         fit_data_psycho_cum{k}(i, 2) = psycho_correct(k,i);
         fit_data_psycho_cum{k}(i, 3) = sum(trials_p); 
     end
end

% now group neuronal data into two groups according to monkey's choice
for h = 1:channelcount
    for k = 1:length(unique_coherence)    % notice that the condition is double than disc_heading
        for i = 1:length(unique_heading)
            select =logical( (heading == unique_heading(i)) & (coherence == unique_coherence(k)) ) ;
            resp{h,k,i} = spike_rates_channel{h}(select);
            resp_mat{h,k}(i) = mean(resp{h,k,i});  % the mean firing rate for each heading
%             resp_mat_std{h,k}(i)= std(resp{h,k,i});
            resp_mat_err{h,k}(i) = std(resp{h,k,i}) / sqrt(repetition);
            % calculate CP, group data based on monkey's choice
            resp_left_choose{h,k,i} = spike_rates_channel{h}(select & (choice == LEFT) );
            resp_right_choose{h,k,i} = spike_rates_channel{h}(select & (choice == RIGHT) );
            if (length(resp_left_choose{h,k,i}) <= 3) || (length(resp_right_choose{h,k,i}) <= 3)   % make sure each condition has at least 3 data values
                Z_Spikes_channel{h}(select) = 9999;   % similar to NaN, just make a mark
            end
        end
        % now across all headings
        resp_left_all{h,k} = Z_Spikes_channel{h}(  (coherence == unique_coherence(k)) & (choice == LEFT) & (Z_Spikes_channel{h}~=9999) );
        resp_right_all{h,k} = Z_Spikes_channel{h}(  (coherence == unique_coherence(k)) & (choice == RIGHT) & (Z_Spikes_channel{h}~=9999) );
%         resp_all{h,k} = Z_Spikes_channel{h}(  (coherence == unique_coherence(k)) & (Z_Spikes_channel{h}~=9999) );
       
        % decide whether ves and vis is congruent tuning. Fit line by linear
        % regression first and compare the sign of each condition to decide whether
        % congruent or opposite, this is used to check whether congruent cells lead
        % to better neuronal performance in combined condition, and vice versa
        [rr,pp] = corrcoef(unique_heading, resp_mat{h,k}(:));
        line_re{h,k} = rr(1,2);
        line_p{h,k} = pp(1,2);
    end
end

% now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
% neurothreshold 
for h = 1:channelcount
    for k = 1 : length(unique_coherence)
        for i = 1 : length(unique_heading)-1   % subtract the 0 heading
            trials_n =logical( (heading == unique_heading(i)) & (coherence == unique_coherence(k)) ) ;
            fit_data_neuro_cum{h,k}(i,3) = sum(trials_n);  % for later function fit use
            if i < (1+length(unique_heading))/2
                % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction
                %      Neuro_correct{k}(i) =  rocN( resp{k,(1+length(unique_heading))/2},resp{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
                % use anti-neuron model instead, so compare each plus and minus headings
                Neuro_correct{h,k}(i) =  rocN( resp{h,k,length(unique_heading)-i+1},resp{h,k,i},100 );
            else
                Neuro_correct{h,k}(i) =  rocN( resp{h,k,length(unique_heading)-i}, resp{h,k,(i+1)},100 );
            end
            %        if  resp_mat{k}(1) < resp_mat{k}(end)
            if line_re{h,k} > 0 % if left heading is not the larger responses, then linear regression will be positive
                % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
                % here we asume if left response larger than right response then asign the left to be preferred direction
                Neuro_correct{h,k}(i) = 1 - Neuro_correct{h,k}(i);
            end
        end
    end
end

%choice probability
for h = 1 : channelcount
    for k = 1 : length(unique_coherence)
        for i = 1 : length(unique_heading)
            if (length(resp_left_choose{h,k,i}) > 3) && (length(resp_right_choose{h,k,i}) > 3)
                CP{h,k}(i) = rocN( resp_left_choose{h,k,i},resp_right_choose{h,k,i},100 );
            else
                CP{h,k}(i) = NaN;
            end
            if  line_re{h,k} > 0
                CP{h,k}(i) = 1 - CP{h,k}(i);
            end
        end
        if (length(resp_left_all{h,k}) > 3) && (length(resp_right_all{h,k}) > 3)
            CP_all(h,k) = rocN( resp_left_all{h,k},resp_right_all{h,k},100 );
        else
            CP_all(h,k) = NaN;
        end
        if  line_re{h,k} > 0
            CP_all(h,k) = 1 - CP_all(h,k);
        end
    end
end

%%%%%% use Wichman's MLE method to estimate threshold and bias
for k = 1:length(unique_coherence) %psychometric
    wichman_psy = pfit(fit_data_psycho_cum{k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
    Thresh_psy{k} = wichman_psy.params.est(2);
    Bias_psy{k} = wichman_psy.params.est(1);
    psy_perf{k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
end
for h = 1:channelcount %neurometric
    for k = 1:length(unique_coherence)
        fit_data_neuro_cum{h,k}(:,1) = unique_heading(unique_heading~=0);
        fit_data_neuro_cum{h,k}(:,2) = Neuro_correct{h,k}(:);
        wichman_neu = pfit(fit_data_neuro_cum{h,k}(:,:),'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
        Thresh_neu{h,k} = wichman_neu.params.est(2);
        % negative and positive infinite value means flat tuning
        if Thresh_neu{h,k}<0 || Thresh_neu{h,k}> 300
            Thresh_neu{h,k} = 300;
            wichman_neu.params.est(2) = 300;
        end
        Bias_neu{h,k} = wichman_neu.params.est(1);
        neu_perf{h,k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%do permutation to test the significance of CP_all{h,k}, re-calculate CP 2000 times
perm_num = 1000;
Z_Spikes_perm_channel = Z_Spikes_channel;
bin = 0.005;
x_bin = 0 : bin : 1;
for h = 1:channelcount
    for k = 1:length(unique_coherence)
        for n = 1 : perm_num
            % temporarily only use near-threshold heading angles where monkey make a guess mainly
            select = logical( (coherence == unique_coherence(k)) & (Z_Spikes_channel{h}~=9999) );
            Z_Spikes_con{h,k} = Z_Spikes_perm_channel{h}( select );
            Z_Spikes_con{h,k} = Z_Spikes_con{h,k}(randperm(length(Z_Spikes_con{h,k})));   % permute spike_rates
            Z_Spikes_perm{h}(select) = Z_Spikes_con{h,k};    % now in spike_rates, the corresponding data were permuted already

            resp_left_all_perm{h,k} = Z_Spikes_perm{h}(  (coherence == unique_coherence(k)) & (choice == LEFT) & (Z_Spikes_channel{h}~=9999) );
            resp_right_all_perm{h,k} = Z_Spikes_perm{h}(  (coherence == unique_coherence(k)) & (choice == RIGHT) & (Z_Spikes_channel{h}~=9999) );

            if  (length(resp_left_all{h,k}) > 3) && (length(resp_right_all{h,k}) > 3)
                CP_all_perm{h,k}(n) = rocN( resp_left_all_perm{h,k}, resp_right_all_perm{h,k},100 );
            else
                CP_all_perm{h,k}(n) = NaN;
            end
            if  line_re{h,k} > 0
                CP_all_perm{h,k}(n) = 1 - CP_all_perm{h,k}(n);
            end

            resp_left_choose_perm = Z_Spikes_perm{h}( (heading == unique_heading(4)) & (coherence == unique_coherence(k)) & (choice == LEFT) & (Z_Spikes_channel{h}~=9999) );
            resp_right_choose_perm = Z_Spikes_perm{h}( (heading == unique_heading(4)) & (coherence == unique_coherence(k)) & (choice == RIGHT) & (Z_Spikes_channel{h}~=9999) );

            if  (length(resp_left_choose{h,k,4}) > 3) && (length(resp_right_choose{h,k,4}) > 3)
                s = warning('off', 'all');
                CP_perm{h}(n) = rocN( resp_left_choose_perm, resp_right_choose_perm,100 );
                warning(s);
            else
                CP_perm{h}(n) = NaN;
            end
            if  line_re{h,1} > 0
                CP_perm{h}(n) = 1 - CP_perm{h}(n);
            end

        end
        % now calculate p value or significant test
        if (length(resp_left_all{h,k}) > 3) && (length(resp_right_all{h,k}) > 3)

            hist_perm{h}(k,:) = hist( CP_all_perm{h,k}(:), x_bin );  % for permutation
            bin_sum = 0;
            n = 0;
            while ( n < (CP_all(h,k)/bin) )
                n = n+1;
                bin_sum = bin_sum + hist_perm{h}(k, n);
                if CP_all(h,k) > 0.5                  % note it's two tail test
                    p(h,k) = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
                else
                    p(h,k) = 2* bin_sum / perm_num;
                end
            end
        else
            p(h,k) = NaN;
        end
        % calculate p value for CP during straight ahead motion
%         if (length(resp_left_choose{h,1,round(length(unique_heading)/2)}) > 3) && (length(resp_right_choose{1,round(length(unique_heading)/2)}) > 3)
%             hist_perm(k,:) = hist( CP_perm{h}(:), x_bin );  % for permutation
%             bin_sum = 0;
%             n = 0;
%             while ( n < (CP{1}(round(length(unique_heading)/2))/bin) )
%                 n = n+1;
%                 bin_sum = bin_sum + hist_perm(k, n);
%                 if CP{k}(round(length(unique_heading)/2)) > 0.5                  % note it's two tail test
%                     pp = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
%                 else
%                     pp = 2* bin_sum / perm_num;
%                 end
%             end
%         else
%             pp = NaN;
%         end
%         p_a{h}(k) = pp;
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot psychometric function
e{1} = 'bo';  f{1} = 'b-';  g{1} = 'bo-';
e{2} = 'rd';  f{2} = 'r-';  g{2} = 'rd-';
e{3} = 'gs';  f{3} = 'g-';  g{3} = 'gs-';
e{4} = 'mx';  f{4} = 'm-';  g{4} = 'mx-';
e{5} = 'kv';  f{5} = 'k-';  g{5} = 'kv-';
e{6} = 'cp';  f{6} = 'c-';  g{6} = 'cp-';
figure(10);
set(10,'Position', [5,25, 980,650], 'Name', 'psycho_neurometic function');
orient landscape;
axes('position',[0.05 0.3, 0.26 0.45]);
legend_txt = [];
for k = 1:length(unique_coherence)
    xi = min(unique_heading) : 0.1 : max(unique_heading);   
 %  plot data in logarithmical space instead of linspace
    plot(unique_heading, psycho_correct(k,:), e{k}, xi, cum_gaussfit(psy_perf{k}, xi),  f{k} );
    xlabel('Heading Angles');   
    ylim([0,1]);
    ylabel('Rightward Choices');
    hold on;
    legend_txt{k*2-1} = num2str(unique_coherence(k));
    legend_txt{k*2} = '';
end
legend(legend_txt,'Location','Best');
title('Psychometric');
clear legend_txt;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot neurometric function

axes('position',[0.36 0.3, 0.26 0.45]);
c = 0;
for h = 1:channelcount
    for k = 1:length(unique_coherence)
        neu_heading = unique_heading(unique_heading~=0);
        xi = min(unique_heading) : 0.1 : max(unique_heading);
        plot(neu_heading, Neuro_correct{h,k}, e{k+c},  xi, cum_gaussfit(neu_perf{h,k}, xi),  f{k+c} );
        xlabel('Heading Angles');
        ylim([0,1]);
        hold on;
        legend_txt{(h-1)*4+k*2-1}= ['ch' num2str(channel_analyze(h)) ' co' num2str(unique_coherence(k))];
        legend_txt{(h-1)*4+k*2}= ' ';
    end
    c=c+2;
end
legend(legend_txt, 'Location','BestOutside');
title('Neurometric');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%neurological raw firing rate
axes('position',[0.7 0.3, 0.26 0.45]);
c = 0;
for h = 1:channelcount
    for k = 1:length(unique_coherence)
        errorbar(unique_heading, resp_mat{h,k}(:), resp_mat_err{h,k}(:),g{k+c} );
        xlabel('Heading Angle (deg)');
        ylabel('Firing rate(spikes/s)');
        xlim([min(unique_heading),max(unique_heading)]);
        hold on;
    end
    c=c+2;
end
title('Firing Rate');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output some text of basic parameters in the figure
axes('position',[0.05,0.8, 0.9,0.15] );
xlim( [0,100] );
ylim( [2,10] );
text(0, 10, FILE);
text(15,10,'repetition =');
text(25,10,num2str(repetition) ); 
text(0,8,'chan');
text(4,8,'coherence');
text(12,8,'Bias:Psy');
text(20,8,'threshold');
text(30,8,'Bias:Neu');
text(40,8,'threshold');
text(52,8,'CP');
text(61,8,'p');
text(72,8,'noisecorr');
% text(68,8,'noise');
% text(76,8,'p');
for k = 1:length(unique_coherence)
        text(6,8-2*k+1, num2str(unique_coherence(k)));
        text(12,8-2*k+1,num2str(Bias_psy{k} ));
        text(20,8-2*k+1,num2str(Thresh_psy{k} ));
        text(67,8-2*k+1,num2str(noise_r_coh(1,k)));
        text(75,8-2*k+1,num2str(noise_p_coh(1,k)));
end
for k = 1:length(unique_coherence)
    for h = 1:channelcount
        text(0, 8-2*(k-1)-h,num2str(channel_analyze(h)));
        text(30,8-2*(k-1)-h,num2str(Bias_neu{h,k} ));
        text(40,8-2*(k-1)-h,num2str(Thresh_neu{h,k} ));
        text(50,8-2*(k-1)-h,num2str(CP_all(h,k)') );
        text(60,8-2*(k-1)-h,num2str(p(h,k)') );
    end
end
axis off;

% -------------------------------------------------------------------------
% Also, write out some summary data to a cumulative summary file
sprint_txt = '%s\t';
for i = 1 : 1000
    sprint_txt = [sprint_txt, ' %1.3f\t'];
end

% buff = sprintf(sprint_txt, FILE, noise_r_coh(1,1),noise_p_coh(1,1),Bias_neu{1,1},Thresh_neu{1,1},Bias_neu{2,1},Thresh_neu{2,1},CP_all{1,1},p{1,1},CP_all{2,1},p{2,1},noise_r_coh(1,2),noise_p_coh(1,2),Bias_neu{1,2},Thresh_neu{1,2},Bias_neu{2,2},Thresh_neu{2,2},CP_all{1,2},p{1,2},CP_all{2,2},p{2,2});
% buff = sprintf(sprint_txt, FILE, noise_r_coh(1,1),noise_p_coh(1,1),CP_all{1,1},p{1,1},CP_all{2,1},p{2,1},noise_r_coh(1,2),noise_p_coh(1,2),CP_all{1,2},p{1,2},CP_all{2,2},p{2,2});
% buff = sprintf(sprint_txt, FILE, Bias_neu{1,1},Thresh_neu{1,1},Bias_neu{2,1},Thresh_neu{2,1},CP_all{1,1},p{1,1},CP_all{2,1},p{2,1},noise_r_coh(1,2),noise_p_coh(1,2),Bias_neu{1,2},Thresh_neu{1,2},Bias_neu{2,2},Thresh_neu{2,2},CP_all{1,2},p{1,2},CP_all{2,2},p{2,2});
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\pairwiseSUSU_Noise_CP.dat'];
% outfile = 'Z:\Users\Sam\Results\pairwiseSUSU_Noise_CP(plexonSort).dat';
%outfile = 'Z:\Users\Piyush\CP_coherence_blocked_for_congruency.dat';
for c=1:channelcount % for each neurons' cp
    buff = sprintf(sprint_txt, FILE, CP_all(c,:),p(c,:));
    outfile = 'Z:\Users\Yun\Alvin_CP.dat';
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid); 
end  

for c=1:num_pairs % for each pair's noise correlation
    buff = sprintf(sprint_txt, FILE, noise_r_coh(c,:),noise_p_coh(c,:));
    outfile = 'Z:\Users\Yun\Alvin_noisecorrelation.dat';
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid); 
end 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t    channel\t coherence\t Bias\t Thresh\t CP\t CPp\t');
%     %     fprintf(fid, 'FILE\t    LoNoise\t LoNoiseP\t LoBiasNu1\t LoThresh1\t LoBiasNu2\t LoThresh2\t LoCP1\t LoCPp1\t LoCP2\t LoCPp2\t HiNoise\t HiNoiseP\t HiBiasNu1\t HiBiasNu2\t HiThresh1\t HiThresh2\t HiCP1\t HiCPp1\t HiCP2\t HiCPp2\t');
%     %     fprintf(fid, 'FILE\t    LoNoise\t LoNoiseP\t LoCP1\t LoCPp1\t LoCP2\t LoCPp2\t HiNoise\t HiNoiseP\t HiCP1\t HiCPp1\t HiCP2\t HiCPp2\t');
%     fprintf(fid, '\r\n');
% end
% for h = 1:channelcount
%     for k = 1:length(unique_coherence)
%         buff = sprintf(sprint_txt, FILE, channel_analyze(h),unique_coherence(k),Bias_neu{h,k},Thresh_neu{h,k},CP_all{h,k},p{h,k},noise_r_coh(:,:));
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%     end
% end
% fclose(fid);
% -------------------------------------------------------------------------
%clear all;
return;

