function bootstrap_curve_fits(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

global Data; %necessary to pass data to subroutines at bottom

Data = data;

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%alpha value for confidence intervals.
alpha = 0.95;

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direction = munique(direction');

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence');

Pref_direction = data.one_time_params(PREFERRED_DIRECTION);

%get the cue validity
cue_val = data.cue_params(CUE_VALIDITY,:,PATCH2);
cue_direc = data.cue_params(CUE_DIREC, :, PATCH1);
unique_cue_val = munique(cue_val');

%compute cue types - 0=neutral, 1=directional, 2=cue_only
cue_type = abs(cue_val); %note that invalid(-1) and valid(+1) are directional
unique_cue_type = munique(cue_type');

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (coherence == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%[direction' coherence' spike_rates' null_trials' select_trials']


figh = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Psychometric Function');
subplot(3, 1, 3);

symbols = {'bo', 'rx', 'g>'};
lines = {'b-', 'r--', 'g:'};
names = {'NoCue','Invalid','Neutral','Valid','CueOnly'};

%% *********** LOGISTIC ANALYSIS ****************************
    
pct_correct = []; N_obs = []; fit_data = [];
logistic_alpha = []; logistic_beta = []; logistic_thresh = [];
legend_str = 'legend(Handl, ';

%this computes the percent of responses in the preferred direction
%combining across redundant conditions within each cue validity.
for i=1:sum(unique_cue_val~=2)
    for j=1:length(unique_direction)
        for k=1:length(unique_coherence)
            ind = k + (j-1)*sum(unique_coherence~=2);
            ok_values = logical( (direction == unique_direction(j)) & (coherence == unique_coherence(k)) ...
                & (cue_val == unique_cue_val(i)) );
            pct_pd(i,ind) = sum(ok_values & (data.misc_params(OUTCOME, :) == CORRECT))/sum(ok_values);
            if (unique_direction(j) ~= Pref_direction)
                pct_pd(i,ind) = 1-pct_pd(i,ind);
            end
        end
    end
end

%plot the raw data
hold on;
for i=1:length(unique_cue_type) %loop through cue type
    signed_coherence = [-unique_coherence' unique_coherence'];
    [sorted_coherence{i}, I{i}] = sort(signed_coherence);
    plot(sorted_coherence{i}, pct_pd(i,I{i}), symbols{i});
end
%keyboard
%now fit these data to logistic function and plot fits
for i=1:sum(unique_cue_val~=2)
    n_obs = sum(cue_val == unique_cue_val(i))./length(unique_coherence).*ones(size(sorted_coherence{i}));
    [logistic_alpha(i) logistic_beta(i)] = logistic_fit([sorted_coherence{i}' pct_pd(i,I{i})' n_obs']);
    logistic_thresh(i) = get_logistic_threshold([logistic_alpha(i) logistic_beta(i)]);
    str = sprintf('%s cue: alpha(slope) = %5.3f, beta(bias) = %5.3f', names{unique_cue_val(i)+3}, logistic_alpha(i), logistic_beta(i));
    hold on
    Handl(i) = plot([min(xlim):1:max(xlim)],logistic_curve([min(xlim):1:max(xlim)],[logistic_alpha(i) logistic_beta(i)]), lines{i});
    legend_str = strcat(legend_str, sprintf(' ''%s'',',names{unique_cue_val(i   )+3}));
    disp(str)
end

xlabel('Coherence x Direction');
ylabel('Fraction Choices in Preferred Direction');
legend_str = strcat(legend_str, ' ''Location'', ''SouthEast'');');
eval(legend_str);
YLim([0 1]);
%comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
%set(gca, 'XScale', 'log');


%% *********** WEIBULL ANALYSIS ****************************

subplot(3, 1, 2);
pct_correct = []; N_obs = []; fit_data = [];
weibull_alpha = []; weibull_beta = []; weibull_gamma = []; weibull_thresh = [];
legend_str = 'legend(Handl, ';

for j = 1:sum(unique_cue_val~=2) %exclude cue_only trials
    for i=1:length(unique_coherence)
        ok_trials = ((coherence == unique_coherence(i)) & (cue_val == unique_cue_val(j))& select_trials);
        correct_trials = (ok_trials & (data.misc_params(OUTCOME, :) == CORRECT) );
        pct_correct(i) = sum(correct_trials)/sum(ok_trials);
        N_obs(i) = sum(ok_trials);
        % data for Weibull fit
        fit_data(i, 1) = unique_coherence(i);
        fit_data(i, 2) = pct_correct(i);
        fit_data(i, 3) = N_obs(i);
    end
    hold on;
    unique_coherence_hack = unique_coherence;
    unique_coherence_hack(1) = 1.01; %allows 0% coherence to appear on semilog plot at location 1.01%
    plot(unique_coherence_hack, pct_correct, symbols{j});
    hold off;
    
    fit_x = unique_coherence(1):0.1: unique_coherence(length(unique_coherence));
    [weibull_alpha(j) weibull_beta(j) weibull_gamma(j)]= weibull_bs_fit(fit_data);
    weibull_thresh(j) = weibull_bs_threshold([weibull_alpha(j) weibull_beta(j) weibull_gamma(j)]);
    monkey_fit_y = weibull_curve(fit_x, [weibull_alpha(j) weibull_beta(j) weibull_gamma(j)]);
 
    hold on;
    Handl(j) = plot(fit_x, monkey_fit_y, lines{j});
    hold off;
    legend_str = strcat(legend_str, sprintf(' ''%s'',',names{unique_cue_val(j)+3}));
end

xlabel('Coherence (% dots)');
ylabel('Fraction Correct');
legend_str = strcat(legend_str, ' ''Location'', ''SouthEast'');');
eval(legend_str);
%YLim([0.4 1]);
%comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
set(gca, 'XScale', 'log');
XLim([1 100]);

%compute fraction correct on cue_only trials
cue_only_trials = (select_trials & (cue_val==2));
cue_only_correct = (cue_only_trials & (data.misc_params(OUTCOME, :) == CORRECT) );
cue_only_pct_corr = sum(cue_only_correct)/sum(cue_only_trials);

% ********** PERFORM BOOSTRAP TESTS *******************************

nboot = 1000;
tic;
trials_outcomes = logical (data.misc_params(OUTCOME,:) == CORRECT);
for i = 1:sum(unique_cue_val~=2)
    for j = 1:nboot
        for k = 1:length(signed_coherence)
            boot_outcomes = []; 
            if (k <= length(signed_coherence)/2) %get direction
                direc = Pref_direction - 180;
            else
                direc = Pref_direction;
            end
            select_boot{i,j,k} = logical( (cue_val == unique_cue_val(i)) & (coherence == abs(signed_coherence(k))) & (direction == direc) );
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
        [bootlog_params{i,j}(1) bootlog_params{i,j}(2)] = logistic_fit([signed_coherence' boot_pct(j,:)' n_obs(j,:)']);
        bootlog_thresh(i,j) = get_logistic_threshold(bootlog_params{i,j});
    end
end    
toc    
keyboard
    
for i = 1:sum(unique_cue_val ~= 2)
    sorted_thresh = sort(bootlog_thresh(i,:));
    bootlog_CI{i} = [sorted_thresh(floor( length(sorted_thresh)*alpha )) ...
                     sorted_thresh(ceil( length(sorted_thresh)*(1-alpha) ))];
end
    
    
    
% %%%%%% Some code for comparing bootstrapped curves to the raw data curve; uncomment and execute.
% 
% fit_x = [-1.*max(unique_coherence):0.1:max(unique_coherence)];
% figure; 
% for i = 1:sum(unique_cue_val~=2)
%     subplot(sum(unique_cue_val~=2),1,i); hold on;
%     for j = 1:nboot
%         plot(fit_x, logistic_curve(fit_x,bootlog_params{i,j}),'k');
%     end
%     plot(fit_x, logistic_curve(fit_x,[logistic_alpha(i) logistic_beta(i)]),'r');
%     title(sprintf('%s',names{unique_cue_val(i)+3}));
% end

% ********** PRINT DATA ON UPPER SUBPLOT **************************

subplot(3, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

start_time = find(data.event_data(1, :, 1) == VSTIM_ON_CD);
stop_time = find(data.event_data(1, :, 1) == VSTIM_OFF_CD);
stim_duration = stop_time - start_time

%now, print out some specific useful info.
xpos = -5; ypos = 10;
font_size = 8;
bump_size = 7;
for j = 1:sum(unique_cue_val~=2)
    line = sprintf('CueStatus = %s, Fit: Weibull, [alpha beta gamma thresh] = [%6.3f %6.3f %6.3f %6.3f]', ...
        names{unique_cue_val(j)+3}, weibull_alpha(j), weibull_beta(j), weibull_gamma(j), weibull_thresh(j));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('CueStatus = %s, Fit: Logistic, alpha(slope) = %6.3f %%, beta(bias) = %6.3f', ...
        names{unique_cue_val(j)+3}, logistic_alpha(j), logistic_beta(j) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
line = sprintf('Monkey: Cue Only Trials, Pct Correct = %6.2f %%', cue_only_pct_corr*100 );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Directions tested: %6.3f, %6.3f deg', unique_direction(1), unique_direction(2) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% line = sprintf('Stimulus Duration: %5d', stim_duration );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%line = sprintf('(3,%0.5g) (6,%0.5g) (12,%0.5g) (24,%0.5g) (48,%0.5g) (96,%0.5g)', ...
%    pct_correct(1), pct_correct(2), pct_correct(3), pct_correct(4), pct_correct(5), pct_correct(6) );
%text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % bootlogistic takes in the bootstrap samples, and returns the slope
% % parameter from the fit.
% function [slope bias] = bootlogistic(bstrp_trials);
% 
% global data;
% TEMPO_Defs;		
% Path_Defs;
% ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
% 
% direction = data.dots_params(DOTS_DIREC, :, PATCH1);
% unique_direction = munique(direction(bstrp_trials)');
% Pref_direction = data.one_time_params(PREFERRED_DIRECTION);
% coherence = data.dots_params(DOTS_COHER, :, PATCH1);
% unique_coherence = munique(coherence');
% signed_coherence = [-munique(coherence'); munique(coherence')]';
% [sorted_coherence I] = sort(signed_coherence);
% 
% for j=1:length(unique_direction)
%     for k=1:length(unique_coherence)
%         ind = k + (j-1)*sum(unique_coherence~=2);
%         ok_values(ind,:) = logical( (direction == unique_direction(j)) & (coherence == unique_coherence(k)) & bstrp_trials');
%         pct_pd(ind) = sum(ok_values(ind) & (data.misc_params(OUTCOME, :) == CORRECT))/sum(ok_values(ind));
%         if (unique_direction(j) ~= Pref_direction)
%             pct_pd(ind) = 1-pct_pd(ind);
%         end
%     end
% end
% 
% [logistic_alpha logistic_beta] = logistic_fit([sorted_coherence' pct_pd(I)' sum(ok_values,2)]);
% slope = logistic_alpha;
% bias = logistic_beta;
% 
% return
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % bootweib takes in the bootstrap samples, and refits the sample to a
% % weibull_bs_fit.  The weibull threshold will be returned.
% function [thresh alpha beta gamma] = bootweib(bstrp_trials);
% 
% global data;
% TEMPO_Defs;		
% Path_Defs;
% ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
% 
% coherence = data.dots_params(DOTS_COHER, :, PATCH1);
% unique_coherence = munique(coherence(bstrp_trials)');
% 
% for i = 1:length(unique_coherence)
%     fit_data(i,1) = unique_coherence(i);
%     c = sum(logical ( bstrp_trials' & (coherence == unique_coherence(i)) & (data.misc_params(OUTCOME,:) == CORRECT) ));
%     t = sum(logical ( bstrp_trials' & (coherence == unique_coherence(i)) ));
%     fit_data(i,2) = c/t;
%     fit_data(i,3) = t;
% end
% 
% [alpha beta gamma] = weibull_bs_fit(fit_data);
% thresh = weibull_bs_threshold([alpha beta gamma]);
% 
% return;
% 
