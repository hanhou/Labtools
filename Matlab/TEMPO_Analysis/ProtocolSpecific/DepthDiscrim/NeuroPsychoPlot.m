%-----------------------------------------------------------------------------------------------------------------------
%-- NeuroPsychoPlot.m -- Plots neurometric and psychometric functions on the same plot; used ROC analysis for neurometrics
%--	GCD, 5/26/00
%-----------------------------------------------------------------------------------------------------------------------
function [neuron_alpha, monkey_alpha] = NeuroPsychoPlot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of horiz. disparities in the dots_params matrix
h_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp = munique(h_disp');

%get the binocular correlations
binoc_corr = data.dots_params(DOTS_BIN_CORR, :, PATCH1);
unique_bin_corr = munique(binoc_corr');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%some stuff done to send data to Izumi: GCD, 2/21/01
%disp('saving...');
%spikes_mat = squeeze(data.spike_data(1,:,:))';
%size(spikes_mat)
%outfid = fopen('temp2.dat', 'w');
%for i=1:size(spikes_mat, 1)
%    fprintf(outfid, '%d ', find(spikes_mat(i,:)>0) );
%    fprintf(outfid, '\r\n');
%end
% fclose(outfid);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (binoc_corr == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(binoc_corr);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%[h_disp' binoc_corr' spike_rates' null_trials' select_trials']

Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

%get the random seed for each trial of the Patch1 dots
%check to see if there is a fixed seed and store this for later if there is.
if (size(data.dots_params,1) >= DOTS_BIN_CORR_SEED)  %for backwards compatibility with old files that lack this
    seeds = data.dots_params(DOTS_BIN_CORR_SEED, :, PATCH1);
    select_fixed_seeds = logical(seeds == data.one_time_params(FIXED_SEED));
else 
    select_fixed_seeds = [];
end
if (sum(select_fixed_seeds) >= 1)
    fixed_seed = data.one_time_params(FIXED_SEED);
else
    fixed_seed = NaN;
end

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Neurometric/Psychometric Comparison');
subplot(2, 1, 2);


% %--------------------------------------------------------
% % spit out mean and variance data for each stimulus condition to a
% % cumulative file, added by GCD 10/18/04
% outfile = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\PopVarMean_summary.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t BinCor\t HDisp\t Mean\t Var\t');
%     fprintf(fid, '\r\n');
%     printflag = 0;
% end
% 
% for i=1:length(unique_bin_corr)
%     for j=1:length(unique_hdisp)
%         temp_trials = ( (h_disp == unique_hdisp(j)) & (binoc_corr == unique_bin_corr(i)) );    
%         temp_mean = mean(spike_rates(temp_trials & select_trials));
%         temp_var = var(spike_rates(temp_trials & select_trials));
%         
%         buff = sprintf('%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', ...
%             FILE, unique_bin_corr(i), unique_hdisp(j), temp_mean, temp_var);
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%        end
% end
% 
% fclose(fid);
%--------------------------------------------------------

% %--------------------------------------------------------
% % spit out data for dprime calcs to a
% % cumulative file, added by GCD 10/19/04
% outfile = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\DPrime_summary.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t BinCor\t HDisp1\t Mean1\t Var1\t HDisp2\t Mean2\t Var2\t Dprime\t');
%     fprintf(fid, '\r\n');
%     printflag = 0;
% end
% 
% for i=1:length(unique_bin_corr)
%     buff = sprintf('%s\t %6.3f\t ', FILE, unique_bin_corr(i));
%     fprintf(fid, '%s', buff);
%     temp_mean=[]; temp_var=[];
%     for j=1:length(unique_hdisp)
%         temp_trials = ( (h_disp == unique_hdisp(j)) & (binoc_corr == unique_bin_corr(i)) );    
%         temp_mean(j) = mean(spike_rates(temp_trials & select_trials));
%         temp_var(j) = var(spike_rates(temp_trials & select_trials));
%         
%         buff = sprintf('%6.3f\t %6.3f\t %6.3f\t ', unique_hdisp(j), temp_mean(j), temp_var(j));
%         fprintf(fid, '%s', buff);
%        end
%     
%     dprime = abs(temp_mean(1)-temp_mean(2))/sqrt(mean([temp_var(1) temp_var(2)]));
%     buff = sprintf('%6.4f\t ', dprime);
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
% end
% 
% fclose(fid);
%--------------------------------------------------------
%return;   %%%!!!! TEMP!!!

%% ********* NEUROMETRIC ANALYSIS ********************
%loop through each binocular correlation levels, and do ROC analysis for each
ROC_values = []; N_obs = [];
for i=1:length(unique_bin_corr)
    
    CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE = 0;
    if (CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE)    
        %Do a regression of spike rates against trial number for each correlation.
        trial_temp = trials((binoc_corr == unique_bin_corr(i)) & select_trials);
        trial_temp = [trial_temp; ones(1,length(trial_temp))];
        spike_temp = spike_rates((binoc_corr == unique_bin_corr(i)) & select_trials);
        [b, bint, r, rint, stats] = regress(spike_temp', trial_temp');    
        spike_rates((binoc_corr == unique_bin_corr(i)) & select_trials) = r';
    end

    pref_trials = ( (h_disp == Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );    
    pref_dist{i} = spike_rates(pref_trials & select_trials);
    null_trials = ( (h_disp ~= Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );    
    null_dist{i} = spike_rates(null_trials & select_trials);
    ROC_values(i) = rocN(pref_dist{i}, null_dist{i}, 100);
    N_obs(i) = length(pref_dist) + length(null_dist);
    % data for Weibull fit
    fit_data(i, 1) = unique_bin_corr(i);
    fit_data(i, 2) = ROC_values(i);
    fit_data(i,3) = N_obs(i);
end

hold on;
Handl(1) = plot(unique_bin_corr, ROC_values, 'ko', 'MarkerFaceColor', 'k');
plot([unique_bin_corr(1) unique_bin_corr(length(unique_bin_corr))], [0.5 0.5], 'k-.');
YLim([0.4 1]);

[neuron_alpha neuron_beta] = weibull_fit(fit_data);
fit_x = unique_bin_corr(1):0.1: unique_bin_corr(length(unique_bin_corr));
neuron_fit_y = weibull_curve(fit_x, [neuron_alpha neuron_beta]);
hold on;
plot(fit_x, neuron_fit_y, 'k-');
hold off;

%% *********** PSYCHOMETRIC ANALYSIS ****************************
pct_correct = []; N_obs = [];
for i=1:length(unique_bin_corr)
    trials = (binoc_corr == unique_bin_corr(i)) & select_trials;
    correct_trials{i} = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
    correct_trials{i} = correct_trials{i}(trials);
    pct_correct(i) = sum(correct_trials{i})/sum(trials);
    N_obs(i) = sum(trials);
    % data for Weibull fit
    fit_data(i, 1) = unique_bin_corr(i);
    fit_data(i, 2) = pct_correct(i);
    fit_data(i,3) = N_obs(i);
end

hold on;
Handl(2) = plot(unique_bin_corr, pct_correct, 'ko');
hold off;

[monkey_alpha monkey_beta] = weibull_fit(fit_data);
monkey_fit_y = weibull_curve(fit_x, [monkey_alpha monkey_beta]);
hold on;
plot(fit_x, monkey_fit_y, 'k--');
hold off;
xlabel('Binocular Correlation (% dots)');
ylabel('Fraction Correct');
legend(Handl, 'Neuron', 'Monkey', 2);

YLim([0.4 1]);
%comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
set(gca, 'XScale', 'log');
XLim([1 100]);

%now calculate whether the slopes and thresholds are significantly different.
%[P_threshold, P_slope] = Weibull_fit_signif_test(ROC_values, pct_correct, unique_bin_corr, N_obs);
%[P_threshold, P_slope, P_conf_thres, P_conf_slope, N_conf_thres, N_conf_slope] = NP_signif_test(pct_correct, ROC_values, correct_trials, pref_dist, null_dist, unique_bin_corr, N_obs);

%titl = sprintf('P threshold = %6.4f, P slope = %6.4f', P_threshold, P_slope);
%title(titl);

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%now, print out some specific useful info.
xpos = 0; ypos = 10;
font_size = 11;
bump_size = 8;
line = sprintf('Neuron: threshold = %6.3f %%, slope = %6.3f', neuron_alpha, neuron_beta );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Monkey: threshold = %6.3f %%, slope = %6.3f', monkey_alpha, monkey_beta );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Disparities tested: %6.3f, %6.3f deg', unique_hdisp(1), unique_hdisp(2) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

%calculate mean and variance for NOVAR and VAR condition. TU 12/18/01
if (isnan(fixed_seed) == 0)
    novar_trials = ( (h_disp == Pref_HDisp) & (binoc_corr == unique_bin_corr(1)) );
    var_trials = ( (h_disp ~= Pref_HDisp) & (binoc_corr == unique_bin_corr(1)) );
    trial_length = ((find(data.event_data(:,:,:)==StopCode)) - (find(data.event_data(:,:,:)==StartCode)))/1000;
    novar_mean = mean(spike_rates(novar_trials & select_trials).*trial_length(novar_trials & select_trials)');
    novar_var = var(spike_rates(novar_trials & select_trials).*trial_length(novar_trials & select_trials)');
    var_mean = mean(spike_rates(var_trials & select_trials).*trial_length(var_trials & select_trials)');
    var_var = var(spike_rates(var_trials & select_trials).*trial_length(var_trials & select_trials)');   
end  


%Now calculate neuronal thresholds for different integration time. TU 12/18/01
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Integration time analysis');
subplot(2, 1, 2);

ROC_values_integ = [];
integ_step = 100;  % increment of integration time, ms
integ_time = integ_step: integ_step: 1500;
for j = 1:length(integ_time) %calculate spike rates for different integration time
    spike_allChan = ComputeSpikeRates(data, length(h_disp), StartCode, StartCode, 30, integ_time(j)+30);
    spike_rates = spike_allChan(1,:);

    for i=1:length(unique_bin_corr)%loop through each binocular correlation levels, and do ROC analysis for each
        pref_trials = ( (h_disp == Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );    
        pref_dist = spike_rates(pref_trials & select_trials);
        null_trials = ( (h_disp ~= Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );    
        null_dist = spike_rates(null_trials & select_trials);
        ROC_values_integ(i) = rocN(pref_dist, null_dist, 100);
        % data for Weibull fit
        fit_data(i, 1) = unique_bin_corr(i);
        fit_data(i, 2) = ROC_values_integ(i);
        fit_data(i,3) = N_obs(i);
    end

    [neuron_alpha_integ(j) neuron_beta_integ(j)] = weibull_fit(fit_data);
end
%normalize thresholds to the last (1500ms) data
norm_thres = neuron_alpha_integ/neuron_alpha_integ(length(neuron_alpha_integ));

hold on;
plot(integ_time, neuron_alpha_integ, 'k-');
hold off;
xlabel('Integration Time (ms)');
ylabel('Neuronal Threshold (%dots)');
%comment out the next 2 lines if you want the plot to be on a LINEAR Y-axis
set(gca, 'YScale', 'log');
YLim([5 100]);
XLim([0 1700]);
%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%Now calculate neuronal thresholds using sliding window. TU 12/03/02
ROC_values_sliding = [];
window_step = 100;  % increment of start of sliding window, ms
window_size = 200;  % size of sliding window, ms
start_time = 0: window_step: 1500 - window_size;
for j = 1:length(start_time) %calculate spike rates for different starting times
    spike_allChan = ComputeSpikeRates(data, length(h_disp), StartCode, StartCode, start_time(j) + 30, start_time(j) + window_size + 30);
    spike_rates = spike_allChan(1,:);
    
    for i=1:length(unique_bin_corr)%loop through each binocular correlation levels, and do ROC analysis for each
        pref_trials = ( (h_disp == Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );    
        pref_dist = spike_rates(pref_trials & select_trials);
        null_trials = ( (h_disp ~= Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );    
        null_dist = spike_rates(null_trials & select_trials);
        ROC_values_sliding(i) = rocN(pref_dist, null_dist, 100);
        % data for Weibull fit
        fit_data(i, 1) = unique_bin_corr(i);
        fit_data(i, 2) = ROC_values_sliding(i);
        fit_data(i,3) = N_obs(i);
    end

    [neuron_alpha_sliding(j) neuron_beta_sliding(j)] = weibull_fit(fit_data);
end
%normalize thresholds to the last (1300-1500ms) data
norm_thres_sliding = neuron_alpha_sliding/neuron_alpha_sliding(length(neuron_alpha_sliding));

%Simulate neuronal thresholds as if the experiments were done with identical dot pattern.
ROC_values_simul = []; 
variance_ratio = 1.5084;
for i=1:length(unique_bin_corr)
    pref_trials = ( (h_disp == Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) & select_trials);
    pref_dist = (spike_rates(pref_trials)-mean(spike_rates(pref_trials)))/sqrt(variance_ratio)+mean(spike_rates(pref_trials));
    null_trials = ( (h_disp ~= Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) & select_trials);    
    null_dist = (spike_rates(null_trials)-mean(spike_rates(null_trials)))/sqrt(variance_ratio)+mean(spike_rates(null_trials));
    ROC_values_simul(i) = rocN(pref_dist, null_dist, 100);
    % data for Weibull fit
    fit_data(i, 1) = unique_bin_corr(i);
    fit_data(i, 2) = ROC_values_simul(i);
    fit_data(i,3) = N_obs(i);
end
[neuron_alpha_simul neuron_beta_simul] = weibull_fit(fit_data);

%----------------------------------------------------------------------------------------
%now, print out data and fits to a file for external plotting purposes (e.g., in Origin)
char = size(PATH,2) - 1;
while PATH(char) ~='\'	%Analysis directory is one branch below Raw Data Dir
    char = char - 1;
end   
PATHOUT = [PATH(1:char) 'Analysis\NeuroPsychoCurves\'];
char = size(FILE,2) - 1;
while FILE(char) ~='.'
    char = char - 1;
end

FILEOUT = [FILE(1:char) 'np_curves'];
fileid = [PATHOUT FILEOUT];
fwriteid = eval(['fopen(fileid, ''w'')']);
fprintf(fwriteid,'BCorrI\tN_fit\tP_fit\tBCorrR\tN_raw\tP_raw\n');
for i=1:length(fit_x)
    fprintf(fwriteid, '%6.3f\t%6.3f\t%6.3f\t', fit_x(i), neuron_fit_y(i), monkey_fit_y(i));
    if (i <= length(unique_bin_corr))
        fprintf(fwriteid, '%6.3f\t%6.3f\t%6.3f\n', unique_bin_corr(i), ROC_values(i), pct_correct(i) );
    else
        fprintf(fwriteid, '\n');
    end
end
fclose(fwriteid);
%----------------------------------------------------------------------------------

%------------------------------------------------------------------------
%write out all relevant parameters to a cumulative text file, GCD 8/08/01
outfile = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\NeuroPsycho_Curve_summary.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    %fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Nthr\t NSlp\t Mthr\t MSlp\t Nthr95\t Nthr95\t Nslp95\t Nslp95\t Pthr95\t Pthr95\t Pslp95\t Pslp95\t Pthr\t Pslp\t DspLo\t DspHi\t Ntrials\t HCorr\t ROChCorr\t NovarMean\t NovarVar\t VarMeam\t VarVar\t');
    fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Nthr\t NSlp\t Mthr\t MSlp\t DspLo\t DspHi\t Ntrials\t HCorr\t ROChCorr\t NovarMean\t NovarVar\t VarMeam\t VarVar\t');
    fprintf(fid, '\r\n');
    printflag = 0;
end
if isnan(fixed_seed)
%    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %10.8f\t %10.8f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5.3f\t %6s\t %6s\t %6s\t %6s\t', ...
%        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%        neuron_alpha,neuron_beta,monkey_alpha,monkey_beta, N_conf_thres(1), N_conf_thres(2), N_conf_slope(1), N_conf_slope(2), P_conf_thres(1), P_conf_thres(2), P_conf_slope(1), P_conf_slope(2), P_threshold, P_slope, unique_hdisp(1), unique_hdisp(2), (1+ EndTrial - BegTrial), unique_bin_corr(length(unique_bin_corr)), ROC_values(length(unique_bin_corr)), '--', '--', '--', '--');
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %7.5f\t %7.5f\t %10.8f\t %10.8f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5.3f\t %6s\t %6s\t %6s\t %6s\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        neuron_alpha,neuron_beta,monkey_alpha,monkey_beta, unique_hdisp(1), unique_hdisp(2), (1+ EndTrial - BegTrial), unique_bin_corr(length(unique_bin_corr)), ROC_values(length(unique_bin_corr)), '--', '--', '--', '--');
else
%    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %10.8f\t %10.8f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', ...
%        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%        neuron_alpha,neuron_beta,monkey_alpha,monkey_beta, N_conf_thres(1), N_conf_thres(2), N_conf_slope(1), N_conf_slope(2), P_conf_thres(1), P_conf_thres(2), P_conf_slope(1), P_conf_slope(2), P_threshold, P_slope, unique_hdisp(1), unique_hdisp(2), (1+ EndTrial - BegTrial), unique_bin_corr(length(unique_bin_corr)), ROC_values(length(unique_bin_corr)), novar_mean, novar_var, var_mean, var_var);
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %7.5f\t %7.5f\t %10.8f\t %10.8f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        neuron_alpha,neuron_beta,monkey_alpha,monkey_beta, unique_hdisp(1), unique_hdisp(2), (1+ EndTrial - BegTrial), unique_bin_corr(length(unique_bin_corr)), ROC_values(length(unique_bin_corr)), novar_mean, novar_var, var_mean, var_var);
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%write out data for integration time analysis TU 12/18/01
outfile2 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\IntegrationTime_summary.dat'];
printflag = 0;
if (exist(outfile2, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile2, 'a');
if (printflag)
    fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
    fprintf(fid, '\r\n');
    printflag = 0;
end

buff = sprintf('%s\t ', FILE);
for i=1:length(norm_thres)
    buff = sprintf('%s %6.4f\t', buff, norm_thres(i));
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

outfile3 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\IntegrationTime_raw.dat'];
printflag = 0;
if (exist(outfile3, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile3, 'a');
if (printflag)
    fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
    fprintf(fid, '\r\n');
    printflag = 0;
end

buff = sprintf('%s\t ', FILE);
for i=1:length(neuron_alpha_integ)
    buff = sprintf('%s %6.4f\t', buff, neuron_alpha_integ(i));
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%write out data for sliding window analysis TU 12/03/02
outfile2 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\SlidingNThres_summary.dat'];
printflag = 0;
if (exist(outfile2, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile2, 'a');
if (printflag)
    fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t');
    fprintf(fid, '\r\n');
    printflag = 0;
end

buff = sprintf('%s\t ', FILE);
for i=1:length(norm_thres_sliding)
    buff = sprintf('%s %6.4f\t', buff, norm_thres_sliding(i));
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

outfile3 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\SlidingNThres_raw.dat'];
printflag = 0;
if (exist(outfile3, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile3, 'a');
if (printflag)
    fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t');
    fprintf(fid, '\r\n');
    printflag = 0;
end

buff = sprintf('%s\t ', FILE);
for i=1:length(neuron_alpha_sliding)
    buff = sprintf('%s %6.4f\t', buff, neuron_alpha_sliding(i));
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%write out parameters for simulation on identical dot patterns, TU 12/28/01
outfile4 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\Simulation_summary.dat'];
printflag = 0;
if (exist(outfile4, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile4, 'a');
if (printflag)
    fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Nthr\t NSlp\t Mthr\t MSlp\t NthrSim\t NslpSim\t DspLo\t DspHi\t Ntrials\t');
    fprintf(fid, '\r\n');
    printflag = 0;
end
buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %7.5f\t %6.3f\t %6.3f\t %4d\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        neuron_alpha,neuron_beta,monkey_alpha,monkey_beta, neuron_alpha_simul, neuron_beta_simul, unique_hdisp(1), unique_hdisp(2), (1+ EndTrial - BegTrial));
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%------------------------------------------------------------------------

return;