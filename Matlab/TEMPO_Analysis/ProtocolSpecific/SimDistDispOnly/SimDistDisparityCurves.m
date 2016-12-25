%-----------------------------------------------------------------------------------------------------------------------
% SimDistDisparityCurves.m -- Module to display center disparity tuning for different surround
%	simulated depths, or vice-versa.  Starting 7/2/02 MLM modified from RelativeDisparityCurves.m
%-----------------------------------------------------------------------------------------------------------------------

function SimDistDisparityCurves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;

symbols = {'ko' 'r*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'r--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity of center in the dots_params matrix
hor_disp_ctr = data.dots_params(DOTS_HDISP,:,PATCH1);

%get the column of values of simulated depth of surround in the dots_params matrix
dist_sim_surr = data.dots_params(DEPTH_DIST_SIM,:,PATCH2);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp_ctr == data.one_time_params(NULL_VALUE)) );

% Surr_Control = no dots in surround, ctr disp varies
surr_control = logical(dist_sim_surr == data.one_time_params(PATCH_OFF));

% Ctr_Control = no dots in ctr, surround varies
ctr_control = logical(hor_disp_ctr == data.one_time_params(PATCH_OFF));


unique_dist_sim_surr = munique(dist_sim_surr(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp_ctr == LEYE_CONTROL) | (hor_disp_ctr == REYE_CONTROL) | (hor_disp_ctr == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp_ctr);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );


figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 150 500 473], 'Name', 'Simulated Depth Tuning Curve');
subplot(2, 1, 2);

all_resp_data = [];

for i=1:length(unique_dist_sim_surr)	%for each different surround simulated depth value, plot a separate disparity tuning curve
    surr_dist_sim_select = logical( (dist_sim_surr == unique_dist_sim_surr(i)) );
    
    plot_x = hor_disp_ctr(surr_dist_sim_select & ~null_trials & ~ctr_control & ~control_trials & select_trials);
    plot_y = spike_rates(surr_dist_sim_select & ~null_trials & ~ctr_control & ~control_trials & select_trials); 
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    hold on;
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
    
    all_resp_data{i}.cdisp = px;
    all_resp_data{i}.resp = py;
    all_resp_data{i}.resp_err = perr;
    
    if unique_dist_sim_surr(i) ~= data.one_time_params(PATCH_OFF)
        surr_control_trials = logical(ctr_control & surr_dist_sim_select);
        if ( sum(surr_control_trials) > 0 )
            surr_resp = spike_rates(surr_control_trials);
            hold on;
            errorbar(max(px)*1.07, mean(surr_resp), std(surr_resp)/sqrt(sum(surr_control)), std(surr_resp)/sqrt(sum(surr_control)), symbols{i});
            text(max(px)*1.12, mean(surr_resp), num2str(unique_disp_surr(i)));
        end
    end
    
    p_value(i) = spk_anova(plot_y, plot_x, px);
    avg_resp(i) = mean(plot_y); 
    hold on;
    H(i) = plot(px, py, symbols{i});
    legend_string{i} = sprintf('%5.2f', unique_dist_sim_surr(i));
end

%write out data in a form for 2-way ANOVA
unique_disp_ctr = munique(hor_disp_ctr(~null_trials & ~ctr_control)');
for i=1:length(unique_dist_sim_surr)	%for each different surround simulated distance value, plot a separate disparity tuning curve
    for j=1:length(unique_disp_ctr)
        select1 = ( logical(dist_sim_surr == unique_dist_sim_surr(i)) & ~surr_control);
        select2 = ( logical(hor_disp_ctr == unique_disp_ctr(j)) & ~ctr_control);
        %[hor_disp_ctr(select1&select2)' hor_disp_surr(select1&select2)' spike_rates(select1&select2)']
        temp = [];
        temp = spike_rates(select1&select2);
        for k=1:length(temp)
            buff = sprintf('%d %d %6.3f', j, i, temp(k));
            %disp(buff);
        end
    end
end   

%perform a 2-way ANOVA in MATLAB using the anovan function
hdisp_group = hor_disp_ctr(~null_trials & ~ctr_control & ~control_trials & select_trials);
simdist_group = dist_sim_surr(~null_trials & ~ctr_control & ~control_trials & select_trials);
resp_var = spike_rates(~null_trials & ~ctr_control & ~control_trials & select_trials); 
[p, atab] = anovan(sqrt(resp_var), {hdisp_group simdist_group}, 2, 3, strvcat('HDISP', 'SIMDIST'));


%write out data in a form for plotting with Origin, etc.
for i=1:length(unique_dist_sim_surr)	
    line = '';
    for j=1:length(unique_disp_ctr)
        buff = sprintf('%7.3f %7.2f %7.3f ', all_resp_data{i}.cdisp(j), all_resp_data{i}.resp(j), all_resp_data{i}.resp_err(j) );
        line = [line buff];
    end
    %disp(line);
end

yl = YLim;
YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
XLabel('Horizontal Disparity of Center(deg)');
YLabel('Response (spikes/sec)');
legend(H, legend_string, -1);

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(px) max(px)];
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');
hold off;

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%now print out useful values for H.Disp specific 
% pmax, pmin, py, 
PrintRelDispData(data, p_value, avg_resp, pmax, pmin, px, null_rate, unique_dist_sim_surr, PATH, FILE);

%------------------------------------------------------------------------
%write out all relevant parameters to a cumulative text file, GCD 10/18/01
outfile = [BASE_PATH 'ProtocolSpecific\SimDistDispOnly\SimDistANOVASummary.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t PvalHDISP\t PvalSDIST\t PvalINTER\t');
    fprintf(fid, '\r\n');
    printflag = 0;
end
buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %10.8f\t %10.8f\t %10.8f\t', ...
    FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
    p);
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%------------------------------------------------------------------------


return;

