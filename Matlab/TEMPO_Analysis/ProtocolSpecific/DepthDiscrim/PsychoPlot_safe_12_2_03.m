%-----------------------------------------------------------------------------------------------------------------------
%-- PsychoPlot.m -- Plots psychometric function
%--	GCD, 7/24/02
%-----------------------------------------------------------------------------------------------------------------------
function [monkey_alpha] = PsychoPlot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of horiz. disparities in the dots_params matrix
h_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp = munique(h_disp');

%get the binocular correlations
binoc_corr = data.dots_params(DOTS_BIN_CORR, :, PATCH1);
unique_bin_corr = munique(binoc_corr');

%get the patch X-center location
x_ctr = data.dots_params(DOTS_AP_XCTR, :, PATCH1);
unique_x_ctr = munique(x_ctr');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (binoc_corr == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(binoc_corr);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%[h_disp' binoc_corr' spike_rates' null_trials' select_trials']

Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Psychometric Function');
subplot(2, 1, 2);

symbols = {'ko', 'kx'};
lines = {'k-', 'k--'};

%% *********** PSYCHOMETRIC ANALYSIS ****************************
pct_correct = []; N_obs = [];
for i=1:length(unique_bin_corr)
    trials = (binoc_corr == unique_bin_corr(i)) & select_trials;
    correct_trials = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
    pct_correct(i) = sum(correct_trials)/sum(trials);
    N_obs(i) = sum(trials);
    % data for Weibull fit
    fit_data(i, 1) = unique_bin_corr(i);
    fit_data(i, 2) = pct_correct(i);
    fit_data(i,3) = N_obs(i);
end

hold on;
Handl(1) = plot(unique_bin_corr, pct_correct, symbols{1});
hold off;

fit_x = unique_bin_corr(1):0.1: unique_bin_corr(length(unique_bin_corr));
[monkey_alpha monkey_beta] = weibull_fit(fit_data);
monkey_fit_y = weibull_curve(fit_x, [monkey_alpha monkey_beta]);
hold on;
plot(fit_x, monkey_fit_y, lines{1});
hold off;


xlabel('Binocular Correlation (% dots)');
ylabel('Fraction Correct');

YLim([0.4 1]);
%comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
set(gca, 'XScale', 'log');
XLim([1 100]);

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

start_time = find(data.event_data(1, :, 1) == VSTIM_ON_CD);
stop_time = find(data.event_data(1, :, 1) == VSTIM_OFF_CD);
stim_duration = stop_time - start_time

%now, print out some specific useful info.
xpos = 0; ypos = 10;
font_size = 11;
bump_size = 8;
line = sprintf('Monkey: threshold = %6.3f %%, slope = %6.3f', monkey_alpha, monkey_beta );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Disparities tested: %6.3f, %6.3f deg', unique_hdisp(1), unique_hdisp(2) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Stimulus Duration: %5d', stim_duration );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

output = 1;
if (output)
    
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    outfile = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\Psycho_Curve_summary.dat'];
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
        monkey_alpha,monkey_beta,unique_hdisp(1), unique_hdisp(2), (1+ EndTrial - BegTrial), unique_bin_corr(length(unique_bin_corr)), stim_duration );
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %------------------------------------------------------------------------
end

return;