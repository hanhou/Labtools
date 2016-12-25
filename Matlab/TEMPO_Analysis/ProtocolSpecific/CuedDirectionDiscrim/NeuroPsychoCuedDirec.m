%-----------------------------------------------------------------------------------------------------------------------
%-- NeuroPsychoCuedDirec.m -- Plots neurometric and psychometric curve for various cue types 
%--	VR, 9/19/05
%-----------------------------------------------------------------------------------------------------------------------
function NeuroPsychoCuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%parameters for bootstrapping to get confidence intervals
nboot = 80;
alpha = .05;

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direction = munique(direction');
Pref_direction = data.one_time_params(PREFERRED_DIRECTION);

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence');
signed_coherence = coherence.*(-1+2.*(direction==Pref_direction));
unique_signed_coherence = [-unique_coherence' unique_coherence'];

%get the cue validity: -1=Invalid; 0=Neutral; 1=Valid; 2=CueOnly
cue_val = data.cue_params(CUE_VALIDITY,:,PATCH2);
unique_cue_val = munique(cue_val');
cue_val_names = {'NoCue','Invalid','Neutral','Valid','CueOnly'};

%get the cue directions
cue_direc = data.cue_params(CUE_DIREC, :, PATCH1);
unique_cue_direc = munique(cue_direc');

%compute cue types - 0=neutral, 1=directional, 2=cue_only
cue_type = abs(cue_val); %note that invalid(-1) and valid(+1) are directional
unique_cue_type = munique(cue_type');

%get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

%get outcome for each trial: 0=incorrect, 1=correct
trials_outcomes = logical (data.misc_params(OUTCOME,:) == CORRECT);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (coherence == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );



%things to do: roc analysis for neurometrics, logistic fits of
%neurometrics/psychometrics, plot, bootstrap
% NeuroMarkers = {'bo','rs','g<','md','c>'}; %Neuron markers will be filled with the colors in NeuroLines
% NeuroLines = {'b','r','g','m','c'};
% NeuroShamLines = {'bo-','rs-','g<-','md-','c>-'}; %to allow both marker and line on legend
% PsychoMarkers = NeuroMarkers;
% PsychoLines = {'b--','r--','g--','m--','c--'};
% PsychoShamLines = {'bo--','rs--','g<--','md--','c>--'};
TempMarkers = {'bo','rs','g<'};
TempLines = {'b','r--','g:'};
TempShamLines = {'bo-','rs--','g<:'};
TempColors = {'b','r','g'};
NeuroMarkers = TempMarkers;
NeuroLines = TempLines;
NeuroShamLines = TempShamLines;
NeuroColors = TempColors;
PsychoMarkers = TempMarkers;
PsychoLines = TempLines;
PsychoShamLines = TempShamLines;
names = {'NoCue','Invalid','Neutral','Valid','CueOnly'};

hlist=figure; 
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Psychometric Function');
subplot(3, 1, 2); hold on;

%% ********* NEUROMETRIC ANALYSIS ********************
%loop through each coherence level per cue val, and do ROC analysis for each
ROC_values = []; N_obs = []; neuron_bootlog_CI = [];
neuron_legend_str = '';
tic

for j=1:sum(unique_cue_val~=2) %exclude CueOnly condition from plot
    for i=1:length(unique_coherence)
        CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE = 0;
        if (CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE)
            %Do a regression of spike rates against trial number for each coherence.
            trial_temp = trials((coherence == unique_coherence(i)) & (cue_val == unique_cue_val(j)) & select_trials);
            trial_temp = [trial_temp; ones(1,length(trial_temp))];
            spike_temp = spike_rates((coherence == unique_coherence(i)) & (cue_val == unique_cue_val(j)) & select_trials);
            [b, bint, r, rint, stats] = regress(spike_temp', trial_temp');
            spike_rates((coherence == unique_coherence(i)) & select_trials) = r';
        end
        pref_trials = ( (direction == Pref_direction) & (coherence == unique_coherence(i)) & (cue_val == unique_cue_val(j)) );
        pref_dist{i} = spike_rates(pref_trials & select_trials);
        null_trials = ( (direction ~= Pref_direction) & (coherence == unique_coherence(i)) & (cue_val == unique_cue_val(j)) );
        null_dist{i} = spike_rates(null_trials & select_trials);
        ROC_values{j}(i) = rocN(pref_dist{i}, null_dist{i}, 100);
        N_obs{j}(i) = length(pref_dist{i}) + length(null_dist{i});
        
        %data for logistic fit - For negative coherences, i'm using 1-ROC(coher); this unfortunately enforces symmetry
        fit_neuron_data{j}(i,1) = -unique_coherence(i);
        fit_neuron_data{j}(i+length(unique_coherence),1) = unique_coherence(i);
        fit_neuron_data{j}(i,2) = 1-ROC_values{j}(i);
        fit_neuron_data{j}(i+length(unique_coherence),2) = ROC_values{j}(i);
        fit_neuron_data{j}(i,3) = N_obs{j}(i);
        fit_neuron_data{j}(i+length(unique_coherence),3) = N_obs{j}(i);
    end
    
    %Plot stuff
    plot(fit_neuron_data{j}(:,1), fit_neuron_data{j}(:,2), NeuroMarkers{j}, 'MarkerFaceColor', NeuroColors{j});
    %plot([unique_bin_corr(1) unique_bin_corr(length(unique_bin_corr))],[0.5 0.5], 'k-.');   %make a line across the plot at y=0.5
    
    [neuron_alpha(j) neuron_beta(j)] = logistic_fit(fit_neuron_data{j});
    neuron_thresh(j) = get_logistic_threshold([neuron_alpha(j) neuron_beta(j)]);
    fit_x = -max(unique_coherence):0.1:max(unique_coherence);
    neuron_fit_y{j} = logistic_curve(fit_x, [neuron_alpha(j) neuron_beta(j)]);
    plot(fit_x, neuron_fit_y{j}, NeuroLines{j});
    n_Handl(j)=plot([-1 1],[-1 -1],NeuroShamLines{j},'MarkerFaceColor',NeuroColors{j}); %note this won't appear on plot
    YLim([0 1]);
    neuron_legend_str = strcat(neuron_legend_str,',''Neuron:',cue_val_names{unique_cue_val(j)+3},'''');
    
    %Bootstrap to generate new roc values for each coherence to generate logistic thresholds
    boot_roc = [];
    for i = 1:nboot
        for k = 1:length(unique_coherence)
            pd_select_boot{j,i,k} = find( select_trials & (cue_val == unique_cue_val(j)) & (coherence == unique_coherence(k)) & (squeeze_angle(direction) == squeeze_angle(Pref_direction)) );
            nd_select_boot{j,i,k} = find( select_trials & (cue_val == unique_cue_val(j)) & (coherence == unique_coherence(k)) & (squeeze_angle(direction) == squeeze_angle(Pref_direction-180)) );
            for m = 1:length(pd_select_boot{j,i,k})
                pd_boot_shuffle = pd_select_boot{j,i,k}(randperm(length(pd_select_boot{j,i,k})));
                pd_boot{j,i,k}(m) = pd_boot_shuffle(1);
                nd_boot_shuffle = nd_select_boot{j,i,k}(randperm(length(nd_select_boot{j,i,k})));
                nd_boot{j,i,k}(m) = nd_boot_shuffle(1);               
            end
            boot_roc(i,k) = rocN(spike_rates(pd_boot{j,i,k}), spike_rates(nd_boot{j,i,k}), 100);
            n_obs(i,k) = length(pd_boot{j,i,k}); %not multiplied by 2
        end
        boot_fit_roc = [1-boot_roc(i,:) boot_roc(i,:)];
        [neuron_bootlog_params{j,i}(1) neuron_bootlog_params{j,i}(2)] = logistic_fit([unique_signed_coherence' boot_fit_roc' [n_obs(i,:) n_obs(i,:)]']);
        neuron_bootlog_thresh(j,i) = get_logistic_threshold(neuron_bootlog_params{j,i}); 

    end
    %now compute confidence intervals
    sorted_thresh = sort(neuron_bootlog_thresh(j,:));
    neuron_bootlog_CI(j,:) = [sorted_thresh(floor( nboot*alpha/2 )) ...
            sorted_thresh(ceil( nboot*(1-alpha/2) ))];
end

xlabel('Coherence x Direction');
ylabel(sprintf('Fraction Choices in\nPreferred Direction'));
neuron_legend_str = strcat('legend(n_Handl',neuron_legend_str, ', ''Location'', ''SouthEast'');');
eval(neuron_legend_str); legend(gca,'boxoff');

toc

%for each combination of validities, use glm to find whether there is a
%significant interaction of cue_validity and coherence
%note that matlab's logistic uses a slightly different parameterization than the one in logistic_func
validity_combo = [ 1 0; 1 -1; 0 -1 ];
for i = 1:size(validity_combo,1)
    yy{i}=[];
    count = 1;
    for j = 1:length(unique_coherence)
        for k = 1:length(validity_combo(i,:))
            yy{i}(count,1) = unique_signed_coherence(j);
            yy{i}(count,2) = validity_combo(i,k);
            yy{i}(count,3) = 1-ROC_values{find(unique_cue_val==validity_combo(i,k))}(j);
			yy{i}(count,4) = 1;
            count = count+1;
            yy{i}(count,1) = -unique_signed_coherence(j);
            yy{i}(count,2) = validity_combo(i,k);
            yy{i}(count,3) = ROC_values{find(unique_cue_val==validity_combo(i,k))}(j);
			yy{i}(count,4) = 1;
%             yy{i}(count,3) = sum((pref_choices == 1) & (signed_coherence == unique_signed_coherence(j)) & (cue_val == validity_combo(i,k)) & select_trials);  % # preferred decisions
%             yy{i}(count,4) = sum((signed_coherence == unique_signed_coherence(j)) & (cue_val == validity_combo(i,k)) & select_trials);		% # trials
            count = count + 1;
        end
    end
    [n_b(i,:), n_dev(i), n_stats{i}] = glmfit([yy{i}(:,1) yy{i}(:,2) yy{i}(:,1).*yy{i}(:,2)],[yy{i}(:,3) yy{i}(:,4)],'binomial');
    P_n_bias(i) = n_stats{i}.p(3);  % P value for bias - always 1 since enforced symmetry eliminates bias
    P_n_slope(i) = n_stats{i}.p(4);	% P value for slope - well, an interaction between validity and signed_coherence
end
i = i+1;
yy{i}=[];
for j = 1:length(unique_coherence)
    for k = 1:sum(unique_cue_val~=2)
        yy{i}(count,1) = unique_coherence(j);
        yy{i}(count,2) = unique_cue_val(k);
        yy{i}(count,3) = ROC_values{k}(j);
        yy{i}(count,4) = 1; %dummy to allow feeding continuous ROC values into binomial
        count = count + 1;
    end
end
[n_b(i,:), n_dev(i), n_stats{i}] = glmfit([yy{i}(:,1) yy{i}(:,2) yy{i}(:,1).*yy{i}(:,2)],[yy{i}(:,3) yy{i}(:,4)],'binomial');
P_n_bias(i) = n_stats{i}.p(3);  % P value for bias - always 1 since enforced symmetry eliminates bias
P_n_slope(i) = n_stats{i}.p(4);	% P value for slope - well, an interaction between validity and signed_coherence

%%% ********* PSYCHOMETRIC ANALYSIS ********************

pct_correct = []; N_obs = []; fit_data = [];
monkey_bootlog_CI = []; monkey_legend_str = '';
subplot(3,1,3); hold on;
%this computes the percent of responses in the preferred direction
%combining across redundant conditions within each cue validity.
for i=1:sum(unique_cue_val~=2)
    for j=1:length(unique_direction)
        for k=1:length(unique_coherence)
            ind = k + (j-1)*length(unique_coherence);
            ok_values = logical( (direction == unique_direction(j)) & (coherence == unique_coherence(k)) ...
                & (cue_val == unique_cue_val(i)) & select_trials );
            pct_pd(i,ind) = sum(ok_values & (data.misc_params(OUTCOME, :) == CORRECT))/sum(ok_values);
            if (unique_direction(j) ~= Pref_direction)
                pct_pd(i,ind) = 1-pct_pd(i,ind);
            end
        end
    end
end

%plot the raw data
for i=1:sum(unique_cue_val~=2) %loop through cue val
    [sorted_coherence{i}, I{i}] = sort(unique_signed_coherence);
    plot(sorted_coherence{i}, pct_pd(i,I{i}), PsychoMarkers{i});
end
%keyboard
%now fit these data to logistic function and plot fits
for i=1:sum(unique_cue_val~=2)
    n_obs = sum(select_trials & (cue_val == unique_cue_val(i)))./length(unique_coherence).*ones(size(sorted_coherence{i}));
    [monkey_alpha(i) monkey_beta(i)] = logistic_fit([sorted_coherence{i}' pct_pd(i,I{i})' n_obs']);
    monkey_thresh(i) = get_logistic_threshold([monkey_alpha(i) monkey_beta(i)]);
    str = sprintf('%s cue: alpha(slope) = %5.3f, beta(bias) = %5.3f', cue_val_names{unique_cue_val(i)+3}, monkey_alpha(i), monkey_beta(i));
    plot([min(xlim):1:max(xlim)],logistic_curve([min(xlim):1:max(xlim)],[monkey_alpha(i) monkey_beta(i)]), PsychoLines{i});
    m_Handl(i) = plot([-1 1], [-1 -1], PsychoShamLines{i});
    monkey_legend_str = strcat(monkey_legend_str,',''Monkey:',cue_val_names{unique_cue_val(i)+3},'''');
end

xlabel('Coherence x Direction');
ylabel(sprintf('Fraction Choices in\nPreferred Direction'));
monkey_legend_str = strcat('legend(m_Handl',monkey_legend_str, ', ''Location'', ''SouthEast'');');
eval(monkey_legend_str); legend(gca,'boxoff');
ylim([0 1]);

% Bootstrap to get 95%CI around threshold behavior 
boot_outcomes = []; 
for i=1:sum(unique_cue_val~=2) %exclude CueOnly condition from plot
    for j=1:nboot
        for k = 1:length(unique_signed_coherence)
            if (k <= length(unique_signed_coherence)/2) %get direction
                direc = Pref_direction - 180;
            else
                direc = Pref_direction;
            end
            select_boot{i,j,k} = logical( select_trials & (cue_val == unique_cue_val(i)) & ...
                (signed_coherence == unique_signed_coherence(k)) );
            behav_select{i,j,k} = trials_outcomes(select_boot{i,j,k});
            for m = 1:length(behav_select)    %loop to generate bootstrap
                boot_shuffle = behav_select{i,j,k}(randperm(length(behav_select{i,j,k})));
                boot_outcomes{i,j,k}(m) = boot_shuffle(1);
            end
            boot_pct(j,k) = sum(boot_outcomes{i,j,k})./length(boot_outcomes{i,j,k});
            if (direc ~= Pref_direction) %for null use 1-pct
                boot_pct(j,k) = 1-boot_pct(j,k);
            end
            n_obs(j,k) = length(boot_outcomes{i,j,k});
        end
        [monkey_bootlog_params{i,j}(1) monkey_bootlog_params{i,j}(2)] = logistic_fit([unique_signed_coherence' boot_pct(j,:)' n_obs(j,:)']);
        monkey_bootlog_thresh(i,j) = get_logistic_threshold(monkey_bootlog_params{i,j});    
    end
    %now compute confidence intervals
    sorted_thresh = sort(monkey_bootlog_thresh(i,:));
    monkey_bootlog_CI(i,:) = [sorted_thresh(floor( nboot*alpha/2 )) ...
            sorted_thresh(ceil( nboot*(1-alpha/2) ))];
end


%for each combination of validities, use glm to find whether there is a
%significant interaction of cue_validity and coherence
%note that matlab's logistic uses a slightly different parameterization than the one in logistic_func
validity_combo = [ 1 0; 1 -1; 0 -1 ];
for i = 1:size(validity_combo,1)
    yy{i}=[];
    count = 1;
    pref_choices = trials_outcomes;
    pref_choices(direction~=Pref_direction) = 1-pref_choices(direction~=Pref_direction);
    for j = 1:length(unique_signed_coherence)
        for k = 1:length(validity_combo(i,:))
            yy{i}(count,1) = unique_signed_coherence(j);
            yy{i}(count,2) = validity_combo(i,k);
            yy{i}(count,3) = sum((pref_choices == 1) & (signed_coherence == unique_signed_coherence(j)) & (cue_val == validity_combo(i,k)) & select_trials);  % # preferred decisions
            yy{i}(count,4) = sum((signed_coherence == unique_signed_coherence(j)) & (cue_val == validity_combo(i,k)) & select_trials);		% # trials
            count = count + 1;
        end
    end
    [p_b(i,:), p_dev(i), p_stats{i}] = glmfit([yy{i}(:,1) yy{i}(:,2) yy{i}(:,1).*yy{i}(:,2)],[yy{i}(:,3) yy{i}(:,4)],'binomial');
    P_p_bias(i) = p_stats{i}.p(3);  % P value for bias
    P_p_slope(i) = p_stats{i}.p(4);	% P value for slope - well, an interaction between validity and signed_coherence
end
i=i+1; %now compute parameters for all three cue validity
yy{i}=[]; count=1;
pref_choices = trials_outcomes;
pref_choices(direction~=Pref_direction) = 1-pref_choices(direction~=Pref_direction);
for j = 1:length(unique_signed_coherence)
    for k = 1:sum(unique_cue_val~=2)
        yy{i}(count,1) = unique_signed_coherence(j);
        yy{i}(count,2) = unique_cue_val(k);
        yy{i}(count,3) = sum((pref_choices == 1) & (signed_coherence == unique_signed_coherence(j)) & (cue_val == unique_cue_val(k)) & select_trials);  % # preferred decisions
        yy{i}(count,4) = sum((signed_coherence == unique_signed_coherence(j)) & (cue_val == unique_cue_val(k)) & select_trials); % # trials
        count = count + 1;
    end
end
[p_b(i,:), p_dev(i), p_stats{i}] = glmfit([yy{i}(:,1) yy{i}(:,2) yy{i}(:,1).*yy{i}(:,2)],[yy{i}(:,3) yy{i}(:,4)],'binomial');
P_p_bias(i) = p_stats{i}.p(3);  % P value for bias
P_p_slope(i) = p_stats{i}.p(4);	% P value for slope - well, an interaction between validity and signed_coherence


bias=[]; slope=[]; crv=[]; se_thresh=[]; se_slope=[];
for i = 1:sum(unique_cue_val~=2)
    % compute the best-fitting curves - Note that the standard matlab logistic uses a slightly different parameterization than used above
    sc_ramp = linspace(min(unique_signed_coherence), max(unique_signed_coherence), 1000);
    
    crv(i,:) = 1./(1+exp(-1*(p_b(4,1)+p_b(4,2).*sc_ramp+p_b(4,3)*unique_cue_val(i)+p_b(4,4).*unique_cue_val(i).*sc_ramp)));
    bias(i) = -1*(p_b(4,1)+p_b(4,3)*unique_cue_val(i))/(p_b(4,2)+p_b(4,4)*unique_cue_val(i)); %signed_coherence where function crosses 0.5 pd choices
    slope(i) = 1/(1+exp(-1*(p_b(4,1)+p_b(4,2).*bias(i)+p_b(4,3)*unique_cue_val(i)+p_b(4,4).*unique_cue_val(i).*bias(i))))^-2 * ...
        exp(-1*(p_b(4,1)+p_b(4,2).*bias(i)+p_b(4,3)*unique_cue_val(i)+p_b(4,4).*unique_cue_val(i).*bias(i))) * ...
        (p_b(4,2)+p_b(4,4)*unique_cue_val(i)); 
        %slope(i) = evaluate derivative of logistic at bias(i) - this
        %disagrees with the 1/(4*alpha) estimate of slope from logistic_func above... but the rank order is correct.
    
%     thresh(1) = -1*b(1)/b(2);		% 50 pct PD threshold
%     slope(1) = b(2);
%    se_thresh = (1/b(2))*sqrt(stats.se(1)^2 + (b(1)^2/b(2)^2)*stats.se(2)^2  );
end

%compute fraction correct on cue_only trials
cue_only_trials = (select_trials & (cue_val==2));
cue_only_correct = (cue_only_trials & (data.misc_params(OUTCOME, :) == CORRECT) );
cue_only_pct_corr = sum(cue_only_correct)/sum(cue_only_trials);

for i = 1:length(unique_cue_val)
    pct_correct(i) = sum(trials_outcomes(cue_val(select_trials)==unique_cue_val(i))) ./ sum(cue_val(select_trials)==unique_cue_val(i));
end

%classify cell as good-bad-ugly based on SLOPE p-values from glm fits
cell_class_names = {'Good','Bad','Ugly','NoBehav'};
GOOD=1; BAD=2; UGLY=3; NB=4;
if ~sum(P_p_slope<.05) %if no significant differences in behavioral slopes
    cell_class_sl = NB;
elseif ~sum(P_n_slope<.05) %if no significant differences in neurometric slopes
    cell_class_sl = BAD;
else
    diff_pairs = find( (P_p_slope<.05) & (P_n_slope<.05) ); %indices of comparisons with significant differences
    if length(diff_pairs)==0  %when different comparisons show differences among neural and behavioral thresholds.
        [dummy, ind_n] = sort(neuron_thresh);
        [dummy, ind_p] = sort(monkey_thresh);
        if isequal(ind_p,ind_n) %in this case, only call it good if the order of all 3 validities are the same for monkey/neuron
            cell_class_sl = GOOD;
        else
            cell_class_sl = BAD;
        end
    else
        cell_class_sl = GOOD; %default
        for i = 1:length(diff_pairs)
            val1 = find(unique_cue_val==validity_combo(i,1));
            val2 = find(unique_cue_val==validity_combo(i,2));
            n_threshes = [neuron_thresh(val1) neuron_thresh(val2)];
            p_threshes = [monkey_thresh(val1) monkey_thresh(val2)];
            if xor(n_threshes(1)>n_threshes(2),p_threshes(1)>p_threshes(2)) %if the bigger value does NOT belong to the same validity,
                cell_class_sl = UGLY;                                          %then call the cell ugly
            end
        end
    end
end
cell_class_names{cell_class_sl};

%classify cell as good-bad-ugly based on BIAS p-values from glm fits
cell_class_names = {'Good','Bad','Ugly','NoBehav'};
GOOD=1; BAD=2; UGLY=3; NB=4;
if ~sum(P_p_bias<.05) %if no significant differences in behavioral slopes
    cell_class_bi = NB;
elseif ~sum(P_n_bias<.05) %if no significant differences in neurometric slopes
    cell_class_bi = BAD;
else
    diff_pairs = find( (P_p_bias<.05) & (P_n_bias<.05) ); %indices of comparisons with significant differences
    if length(diff_pairs)==0  %when different comparisons show differences among neural and behavioral thresholds.
        [dummy, ind_n] = sort(neuron_beta);
        [dummy, ind_p] = sort(monkey_beta);
        if isequal(ind_p,ind_n) %in this case, only call it good if the order of all 3 validities are the same for monkey/neuron
            cell_class_bi = GOOD;
        else
            cell_class_bi = BAD;
        end
    else
        cell_class_bi = GOOD; %default
        for i = 1:length(diff_pairs)
            val1 = find(unique_cue_val==cuedir_combo(i,1));
            val2 = find(unique_cue_val==cuedir_combo(i,2));
            n_biases = [neuron_beta(val1) neuron_beta(val2)];
            p_biases = [monkey_beta(val1) monkey_beta(val2)];
            if xor(n_biases(1)>n_biases(2),p_biases(1)>p_biases(2)) %if the bigger value does NOT belong to the same validity,
                cell_class_bi = UGLY;                                       %then call the cell ugly
            end
        end
    end
end
cell_class_names{cell_class_bi};

%% ********************** PRINT INFO *****************************
%now, print out some useful information in the upper subplot
subplot(3, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);


%now, print out some specific useful info.
xpos = -10; ypos = 25;
font_size = 8;
bump_size = 6;
% for j = 1:sum(unique_cue_val~=2)
%     line = sprintf('CueStatus = %s, slope = %6.2f, bias = %6.2f%%', ...
%         names{unique_cue_val(j)+3}, monkey_alpha(j), monkey_beta(j) );
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% end
% line = sprintf(FILE);
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
for j = 1:sum(unique_cue_val~=2)
    line = sprintf('Neuron: CueStatus = %s, thresh = %6.2f%%, %d%% CI = [%6.2f%% %6.2f%%]', ...
        cue_val_names{unique_cue_val(j)+3}, neuron_thresh(j), 100*(1-alpha), neuron_bootlog_CI(j,1), neuron_bootlog_CI(j,2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
% line = sprintf('Monkey Thresholds:');
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
for j = 1:sum(unique_cue_val~=2)
    line = sprintf('Monkey: CueStatus = %s, thresh = %6.2f%%, %d%% CI = [%6.2f%% %6.2f%%]', ...
        cue_val_names{unique_cue_val(j)+3}, monkey_thresh(j), 100*(1-alpha), monkey_bootlog_CI(j,1), monkey_bootlog_CI(j,2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
line = sprintf('Pct Correct:');
for j = 1:length(unique_cue_val)
    line = strcat(line, sprintf(' %s = %4.2f%%;',names{unique_cue_val(j)+3},pct_correct(j)*100));
end
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Directions tested: %6.3f, %6.3f deg', unique_direction(1), unique_direction(2) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% line = sprintf('Stimulus Duration: %5d', stim_duration );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% line = pct_str;
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

output = 0; %basic stuff... thresholds, alpha,beta,%correct,etc
output2 = 1; 
if (output)
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, VR 11/21/05
    outfile = [BASE_PATH 'ProtocolSpecific\CuedDirectionDiscrim\NeuroPsycho_Curve_summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t InvPct\t NeuPct\t ValPct\t CuePct\t P_InvTh\t P_NeuTh\t P_ValTh\t P_InvCILow\t P_InvCIHi\t P_NeuCILow\t P_NeuCIHi\t P_ValCILow\t P_ValCIHi\t P_InvSl\t P_InvBi\t P_NeuSl\t P_NeuBi\t P_ValSl\t P_ValBi\t N_InvTh\t N_NeuTh\t N_ValTh\t N_InvCILow\t N_InvCIHi\t N_NeuCILow\t N_NeuCIHi\t N_ValCILow\t N_ValCIHi\t N_InvSl\t N_InvBi\t N_NeuSl\t N_NeuBi\t N_ValSl\t N_ValBi\t MaxCoher\t Ntrials\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    inval = find(unique_cue_val==-1);   neu = find(unique_cue_val==0);   val = find(unique_cue_val==1);  cue = find(unique_cue_val==2);
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.1f\t %4d\t',...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        pct_correct(inval)*100,pct_correct(neu)*100,pct_correct(val)*100,pct_correct(cue)*100,...
        monkey_thresh(inval),monkey_thresh(neu),monkey_thresh(val), monkey_bootlog_CI(inval,1),monkey_bootlog_CI(inval,2),...
        monkey_bootlog_CI(neu,1),monkey_bootlog_CI(neu,2),monkey_bootlog_CI(val,1),monkey_bootlog_CI(val,2),...
        monkey_alpha(inval),monkey_beta(inval),monkey_alpha(neu),monkey_beta(neu),monkey_alpha(val),monkey_beta(val),...
        neuron_thresh(inval),neuron_thresh(neu),neuron_thresh(val), neuron_bootlog_CI(inval,1),neuron_bootlog_CI(inval,2),...
        neuron_bootlog_CI(neu,1),neuron_bootlog_CI(neu,2),neuron_bootlog_CI(val,1),neuron_bootlog_CI(val,2),...
        neuron_alpha(inval),neuron_beta(inval),neuron_alpha(neu),neuron_beta(neu),neuron_alpha(val),neuron_beta(val),...
        max(unique_coherence),(1+EndTrial-BegTrial) );
    %buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.4f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5d\t', ...
    %    FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
    %    monkey_alpha,monkey_beta,unique_direction(1), unique_direction(2), (1+ EndTrial - BegTrial), unique_coherence(length(unique_coherence)), stim_duration );
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %------------------------------------------------------------------------
end
if (output2) %categorization
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, VR 11/21/05
    outfile = [BASE_PATH 'ProtocolSpecific\CuedDirectionDiscrim\NeuroPsycho_Classify_summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t Class_Sl\t Class_Bi\t P_sl_pvals\t SlPpVN\t SlPpVI\t SlPpNI\t P_bi_pvals\t BiPpVN\t BiPpVI\t BiPpNI\t N_sl_pvals\t SlNpVN\t SlNpVI\t SlNpNI\t N_bi_pvals\t BiNpVN\t BiNpVI\t BiNpNI\t P_ValTh\t P_NeuTh\t P_InvTh\t P_ValSl\t P_ValBi\t P_NeuSl\t P_NeuBi\t P_InvSl\t P_InvBi\t N_ValTh\t N_NeuTh\t N_InvTh\t N_ValSl\t N_ValBi\t N_NeuSl\t N_NeuBi\t N_InvSl\t N_InvBi\t MaxCoher\t Ntrials\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    val = find(unique_cue_val==1);  neu = find(unique_cue_val==0);
    inv = find(unique_cue_val==-1); cue = find(unique_cue_val==2);
    buff = sprintf('%s\t %s\t %s\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.5f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.1f\t %4d\t',...
        FILE, cell_class_names{cell_class_sl}, cell_class_names{cell_class_bi},...
        circshift(P_p_slope,[1 1]), circshift(P_p_bias,[1 1]), circshift(P_n_slope,[1 1]), circshift(P_n_bias,[1 1]), ...
        monkey_thresh(val),monkey_thresh(neu),monkey_thresh(inv), ...
        monkey_alpha(val),monkey_beta(val), monkey_alpha(neu),monkey_beta(neu), monkey_alpha(inv),monkey_beta(inv), ...
        neuron_thresh(val),neuron_thresh(neu),neuron_thresh(inv), ...
        neuron_alpha(val),neuron_beta(val), neuron_alpha(neu),neuron_beta(neu), neuron_alpha(inv),neuron_beta(inv), ...
        max(unique_coherence),(1+EndTrial-BegTrial) );
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %------------------------------------------------------------------------
end


SAVE_FIGS = 0;
if (SAVE_FIGS)
    saveas(hlist, sprintf('%s_NP_curves.fig',FILE),'fig');
end

return