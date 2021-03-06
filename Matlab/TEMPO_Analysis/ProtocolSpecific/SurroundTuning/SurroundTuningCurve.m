%-----------------------------------------------------------------------------------------------------------------------
%-- Surround Tuning curve will plot a polar plot for each set of surrounds -JDN 3/28/00
%-----------------------------------------------------------------------------------------------------------------------
function SurroundTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

%get the column of disparity values for the surround patches
disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (disp == data.one_time_params(NULL_VALUE)) );
patch_off_condition = logical((disp == data.one_time_params(PATCH_OFF)));

unique_disp = munique(disp(~null_trials & ~patch_off_condition)');

%get the column of values of offset angle of the surround patches
direc = data.dots_params(DOTS_DIREC,BegTrial:EndTrial,PATCH1);
unique_direc = munique(direc(~null_trials & ~patch_off_condition)');

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (disp == LEYE_CONTROL) | (disp == REYE_CONTROL) | (disp == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stringarray = [];

%now, print out some useful information in the upper subplot
%subplot(2, 1, 1);
%PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

figure
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Surround Tuning Curve');
subplot(212);
if (length(unique_disp) ~= 1)
    for i = 1:length(unique_direc)

        direc_select = logical(direc == unique_direc(i));

        plot_x = disp(direc_select & ~null_trials & ~control_trials & select_trials);
        plot_y = spike_rates(direc_select & ~null_trials & ~control_trials & select_trials);

        %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);

        hold off;

    end
else
    for i = 1:length(unique_disp)
        disp_select = logical(disp == unique_disp(i));

        plot_x = direc(disp_select & ~null_trials & ~control_trials & select_trials);
        plot_y = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);

        %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);

        hold off;
    end
end

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(px) max(px)];
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');

%now, get the firing rate for RF only condition trials and add spontaneous rate to plot
rf_x = [min(px) max(px)];
rf_rate = mean(data.spike_rates(SpikeChan, patch_off_condition & select_trials));
rf_y = [rf_rate rf_rate];
hold on;
plot(rf_x, rf_y, 'b--');

hold off;


grid on

%height = axis;
%yheight = height(4);
%string = sprintf('Ap Size = %4.1f', unique_ap_size(select_ap_size));
%text(height(1)+2, 0.95*yheight, string, 'FontSize', 8);
%string = sprintf('Mag Disp = %1.3f', unique_disp(k));
%text(height(1)+2, 0.85*yheight, string, 'FontSize', 8);
%hold off

XLabel('Surround Horizontal Disparity(deg)');
YLabel('Response (spikes/sec)');

p_value = spk_anova(plot_y, plot_x, unique_disp);

%% ********************** PRINT INFO *****************************
%now, print out some useful information in the upper subplot
subplot(211);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%now, print out some specific useful info.
xpos = -10; ypos = 25;
font_size = 8;
bump_size = 6;

line = sprintf('ANOVA p value: %6.5f', p_value );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

return;
