%-----------------------------------------------------------------------------------------------------------------------
%-- PsychCuedDirec.m -- Plots psychometric curve for various cue types 
%--	VR, 6/2/05
%-----------------------------------------------------------------------------------------------------------------------
function PsychCuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);


TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

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
% cue_val = zeros(size(cue_direc));
% cue_val( (cue_direc == Pref_direction) | (cue_direc == Pref_direction - 360) ) = 1;
% cue_val(cue_direc == Pref_direction - 180) = -1;
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


figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Psychometric Function');
subplot(2, 1, 2);

symbols = {'bo', 'rx', 'g>'};
lines = {'b-', 'r--', 'g:'};
names = {'NoCue','Invalid','Neutral','Valid','CueOnly'};
%% *********** PSYCHOMETRIC ANALYSIS ****************************
    
pct_correct = []; N_obs = []; fit_data = [];
monkey_alpha = []; monkey_beta = [];
legend_str = 'legend(Handl, ';

% %this computes the percent of responses in the preferred direction
% %combining across redundant conditions within each cue validity.
% for i=1:sum(unique_cue_val~=2)
%     for j=1:length(unique_direction)
%         for k=1:length(unique_coherence)
%             ind = k + (j-1)*sum(unique_coherence~=2);
%             ok_values = logical( (direction == unique_direction(j)) & (coherence == unique_coherence(k)) ...
%                 & (cue_val == unique_cue_val(i)) );
%             pct_pd(i,ind) = sum(ok_values & (data.misc_params(OUTCOME, :) == CORRECT))/sum(ok_values);
%             if (unique_direction(j) ~= Pref_direction)
%                 pct_pd(i,ind) = 1-pct_pd(i,ind);
%             end
%         end
%     end
% end
% 
% %plot the raw data
% hold on;
% for i=1:length(unique_cue_type) %loop through cue type
%     signed_coherence = [-unique_coherence' unique_coherence'];
%     [sorted_coherence{i}, I{i}] = sort(signed_coherence);
%     plot(sorted_coherence{i}, pct_pd(i,I{i}), symbols{i});
% end
% %keyboard
% %now fit these data to logistic function and plot fits
% for i=1:sum(unique_cue_val~=2)
%     n_obs = sum(cue_val == unique_cue_val(i))./length(unique_coherence).*ones(size(sorted_coherence{i}));
%     [monkey_alpha(i) monkey_beta(i)] = logistic_fit([sorted_coherence{i}' pct_pd(i,I{i})' n_obs']);
%     str = sprintf('%s cue: alpha(slope) = %5.3f, beta(bias) = %5.3f', names{unique_cue_val(i)+3}, monkey_alpha(i), monkey_beta(i));
%     hold on
%     Handl(i) = plot([min(xlim):1:max(xlim)],logistic_curve([min(xlim):1:max(xlim)],[monkey_alpha(i) monkey_beta(i)]), lines{i});
%     legend_str = strcat(legend_str, sprintf(' ''%s'',',names{unique_cue_val(i   )+3}));
%     disp(str)
% end
% 
% xlabel('Coherence x Direction');
% ylabel('Fraction Choices in Preferred Direction');
% legend_str = strcat(legend_str, ' ''Location'', ''SouthEast'');');
% eval(legend_str);
% YLim([0 1]);
% %comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
% %set(gca, 'XScale', 'log');
% XLim([-100 100]);





for j = 1:sum(unique_cue_val~=2) %exclude cue_only trials
    for i=1:length(unique_coherence)
        trials = ((coherence == unique_coherence(i)) & (cue_val == unique_cue_val(j))& select_trials);
        correct_trials = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
        pct_correct(j,i) = sum(correct_trials)/sum(trials);
        N_obs(i) = sum(trials);
        % data for Weibull fit
        fit_data(i, 1) = unique_coherence(i);
        fit_data(i, 2) = pct_correct(j,i);
        fit_data(i, 3) = N_obs(i);
    end
    hold on;
    unique_coherence_hack = unique_coherence;
%     unique_coherence_hack(1) = 1.01; %allows 0% coherence to appear on semilog plot at location 1.01%
    plot(unique_coherence_hack, pct_correct(j,:), symbols{j});
    hold off;
    
    fit_x = unique_coherence(1):0.1: unique_coherence(length(unique_coherence));
    %[monkey_alpha(j) monkey_beta(j)]= weibull_fit(fit_data);
    %monkey_fit_y = weibull_curve(fit_x, [monkey_alpha(j) monkey_beta(j)]);
    [monkey_alpha(j) monkey_beta(j)]= logistic_fit(fit_data);
    monkey_fit_y(j,:) = logistic_curve(fit_x, [monkey_alpha(j) monkey_beta(j)]);
 
    hold on;
    Handl(j) = plot(fit_x, monkey_fit_y(j,:), lines{j});
    hold off;
    legend_str = strcat(legend_str, sprintf(' ''%s'',',names{unique_cue_val(j)+3}));
end

xlabel('Coherence (% dots)');
ylabel('Fraction Correct');
legend_str = strcat(legend_str, ' ''Location'', ''SouthEast'');');
eval(legend_str);
%YLim([0.4 1]);
%comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
% set(gca, 'XScale', 'log');
% XLim([1 100]);

%some temporary stuff to save out raw data to a text file so cleaner graphs
%can be made in origin
save_psych1d = 1;
if save_psych1d
    savename = sprintf('Z:\\Data\\Tempo\\Baskin\\Analysis\\1D-Psychophysics\\ShortDurDelayCells\\Psy1D-%s.txt',FILE);
    temp = unique_coherence';
    save(savename, 'temp', '-ascii'); %this saves a row containing the values of unique_coherence
    temp = pct_correct;
    save(savename, 'temp', '-ascii', '-append'); %this saves 3 lines - each containing the pct correct at each coherence for a single validity (invalid, neutral, valid)
    temp = fit_x;
    save(savename, 'temp', '-ascii', '-append'); %this saves the x-values of the best fit logistic curves
    temp = monkey_fit_y;
    save(savename, 'temp', '-ascii', '-append'); %this saves the three lines (one for each validity) of y-values of the best fit logistic curves
end

%compute and save RTs
for i = 1:length(trials)
    rt(i) = find(data.event_data(1,:,i)==SACCADE_BEGIN_CD) - find(data.event_data(1,:,i)==TARGS_ON_CD);
end
for i = 1:length(unique_coherence)
    for k = 1:length(unique_cue_val)
        if unique_cue_val(k) == 2, select = select_trials & (cue_val == 2); %kluge to combine across coher for cueonly 
        else, select = select_trials & (coherence == unique_coherence(i)) & (cue_val == unique_cue_val(k));
        end
        mean_rt(i,k) = mean(rt(select));
        std_rt(i,k) = std(rt(select));
    end
    select = select_trials & (coherence == unique_coherence(i));
    mean_rt(i,k+1) = mean(rt(select));
    std_rt(i,k+1) = std(rt(select));
end
save_rt = 1;
if save_rt
    savename = sprintf('Z:\\Data\\Tempo\\Baskin\\Analysis\\CuedDirecDiscrim_RT\\ShortDurDelayCells\\RT-%s.txt',FILE(1:8));
    temp = unique_coherence';
    save(savename, 'temp', '-ascii'); %this saves a row containing the values of unique_coherence
    temp = mean_rt';
    save(savename, 'temp', '-ascii', '-append'); %this saves 5 lines - mean RT given coherence for invalid, neutral, valid, cueonly, and group avg
    temp = std_rt';
    save(savename, 'temp', '-ascii', '-append'); %this saves 5 lines - std RT given coherence for invalid, neutral, valid, cueonly, and group avg
end

%compute fraction correct on cue_only trials
cue_only_trials = (select_trials & (cue_val==2));
cue_only_correct = (cue_only_trials & (data.misc_params(OUTCOME, :) == CORRECT) );
cue_only_pct_corr = sum(cue_only_correct)/sum(cue_only_trials);

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

start_time = find(data.event_data(1, :, 1) == VSTIM_ON_CD);
stop_time = find(data.event_data(1, :, 1) == VSTIM_OFF_CD);
stim_duration = stop_time - start_time

%now, print out some specific useful info.
xpos = 0; ypos = 10;
font_size = 10;
bump_size = 8;
for j = 1:sum(unique_cue_val~=2)
    line = sprintf('Monkey: CueStatus = %s, threshold = %6.3f %%, slope = %6.3f', names{unique_cue_val(j)+3}, monkey_alpha(j), monkey_beta(j) );
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

output = 0;
if (output)
    
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    outfile = [BASE_PATH 'ProtocolSpecific\CuedDirecDiscrim\Psycho_Curve_summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t MThr\t MSlp\t DspLo\t DspHi\t Ntrials\t HCorr\t Durat\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.4f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5d\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        monkey_alpha,monkey_beta,unique_direction(1), unique_direction(2), (1+ EndTrial - BegTrial), unique_coherence(length(unique_coherence)), stim_duration );
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %------------------------------------------------------------------------
end

return;