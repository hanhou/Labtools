%-----------------------------------------------------------------------------------------------------------------------
%-- NeuroPsychoPlot.m -- Plots neurometric and psychometric functions on the same plot; used ROC analysis for neurometrics
%--	GCD, 5/26/00
%-----------------------------------------------------------------------------------------------------------------------
function NeuroPsychoPlot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%defns like CORRECT
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

%get the column of values of horiz. disparities in the dots_params matrix
h_disp_p1 = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp_p1 = munique(h_disp_p1');

%get the column of values of horiz. disparities in the dots_params matrix
h_disp_p4 = data.dots_params(DOTS_HDISP,:,PATCH4);
unique_hdisp_p4 = munique(h_disp_p4');

%compute the unsigned horizontal disparity
unsigned_hdisp = abs(h_disp_p1-h_disp_p4);
unsigned_hdisp = (round(10000 * unsigned_hdisp)) / 10000;
unique_unsigned_hdisp = munique(unsigned_hdisp');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (h_disp_p1 == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(h_disp_p1);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%[h_disp' binoc_corr' spike_rates' null_trials' select_trials']

Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 150 500 473], 'Name', 'Neurometric/Psychometric Comparison');
subplot(2, 1, 2);

%% ********* NEUROMETRIC ANALYSIS ********************
%loop through each unsigned disparity levels, and do ROC analysis for each
ROC_values = []; N_obs = []; table = [];
for i=1:length(unique_unsigned_hdisp)
    if (unique_unsigned_hdisp(i) ~= 0)
        pref_trials = ( (sign(h_disp_p1-h_disp_p4) == sign(Pref_HDisp-h_disp_p4)) & (unsigned_hdisp == unique_unsigned_hdisp(i)) );
        pref_dist{i} = spike_rates(pref_trials & select_trials);
        null_trials = ( (sign(h_disp_p1-h_disp_p4) ~= sign(Pref_HDisp-h_disp_p4)) & (unsigned_hdisp == unique_unsigned_hdisp(i)) );
        null_dist{i} = spike_rates(null_trials & select_trials);
        ROC_values(i) = rocN(pref_dist{i}, null_dist{i}, 100);
        N_obs(i) = length(pref_dist{i}) + length(null_dist{i});
        % data for Weibull fit
        fit_data(i, 1) = unique_unsigned_hdisp(i);
        fit_data(i, 2) = ROC_values(i);
        fit_data(i,3) = N_obs(i);
        %table = [table pref_dist' null_dist'];
    end
end

%now, print out raw data (sorted by disparity) for external plotting purposes (e.g., in Origin)
%i = size(PATH,2) - 1;
%while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
%    i = i - 1;
%end   
%PATHOUT = [PATH(1:i) 'Analysis\NeuroPsychoCurves\'];
%i = size(FILE,2) - 1;
%while FILE(i) ~='.'
%    i = i - 1;
%end
%FILEOUT = [FILE(1:i) 'resp_by_disp'];
%fileid = [PATHOUT FILEOUT];
%fwriteid = eval(['fopen(fileid, ''w'')']);
%for i=1:length(unique_unsigned_hdisp)
%    fprintf(fwriteid,'Pref%d\tNull%d\t', i, i);
%end
%fprintf(fwriteid,'\n');
%for j=1:size(table,1)  %loop through trials
%    junk = table(j,:)
%    fprintf(fwriteid,'%6.2f ', junk);
%    fprintf(fwriteid,'\n');
%end
%fclose(fwriteid);
%save 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\Stereoacuity\Resp_by_disp.dat' table -ASCII

hold on;
Handl(1) = plot(unique_unsigned_hdisp, ROC_values, 'ko', 'MarkerFaceColor', 'k');
plot([unique_unsigned_hdisp(1) unique_unsigned_hdisp(length(unique_unsigned_hdisp))], [0.5 0.5], 'k-.');
YLim([0.4 1]);

[neuron_alpha neuron_beta] = weibull_fit(fit_data);
fit_x = unique_unsigned_hdisp(1):0.0001: unique_unsigned_hdisp(length(unique_unsigned_hdisp));
neuron_fit_y = weibull_curve(fit_x, [neuron_alpha neuron_beta]);
hold on;
plot(fit_x, neuron_fit_y, 'k-');
hold off;

%% *********** PSYCHOMETRIC ANALYSIS ****************************
pct_correct = []; N_obs = [];
for i=1:length(unique_unsigned_hdisp)
    if (unique_unsigned_hdisp(i) ~= 0)
        trials = (unsigned_hdisp == unique_unsigned_hdisp(i)) & select_trials;
        correct_trials{i} = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
        correct_trials{i} = correct_trials{i}(trials);
        pct_correct(i) = sum(correct_trials{i})/sum(trials);
        N_obs(i) = sum(trials);
        % data for Weibull fit
        fit_data(i, 1) = unique_unsigned_hdisp(i);
        fit_data(i, 2) = pct_correct(i);
        fit_data(i,3) = N_obs(i);
    end
end

hold on;
Handl(2) = plot(unique_unsigned_hdisp, pct_correct, 'ko');
hold off;
[monkey_alpha monkey_beta] = weibull_fit(fit_data);
monkey_fit_y = weibull_curve(fit_x, [monkey_alpha monkey_beta]);
hold on;
plot(fit_x, monkey_fit_y, 'k--');
hold off;
temp = [fit_x' monkey_fit_y'];
save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\Stereoacuity\weibul-fit.txt', '-ascii', 'temp');
xlabel('Binocular Disparity (deg)');
ylabel('Fraction Correct');
legend(Handl, 'Neuron', 'Monkey', 2);
YLim([0.4 1]);
%comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
set(gca, 'XScale', 'log');
XLim([.001 1]);

%now calculate whether the slopes and thresholds are significantly different.
%[P_threshold, P_slope] = Weibull_fit_signif_test(ROC_values, pct_correct, unique_unsigned_hdisp, N_obs);
%[P_threshold, P_slope, P_conf_thres, P_conf_slope, N_conf_thres, N_conf_slope] = NP_signif_test(pct_correct, ROC_values, correct_trials, pref_dist, null_dist, unique_unsigned_hdisp, N_obs);
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
line = sprintf('Pedestal disparity: %6.3f deg', unique_hdisp_p4(1) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = '  '; 
for i = 1:length(unique_unsigned_hdisp)
    line = [line sprintf('(%5.4f,%4.2f) ',unique_unsigned_hdisp(i), pct_correct(i))];
end
text(xpos,ypos,line,'FontSize',8);		ypos = ypos - bump_size;


%Now calculate neuronal thresholds for different integration time. TU 08/27/03
ROC_values_integ = [];
integ_step = 100;  % increment of integration time, ms
integ_time = integ_step: integ_step: 1500;
for j = 1:length(integ_time) %calculate spike rates for different integration time
    spike_allChan = ComputeSpikeRates(data, length(h_disp_p1), StartCode, StartCode, 30, integ_time(j)+30);
    spike_rates = spike_allChan(1,:);

    for i=1:length(unique_unsigned_hdisp)%loop through each unsigned center disparity, and do ROC analysis for each
        if (unique_unsigned_hdisp(i) ~= 0)
            pref_trials = ( (sign(h_disp_p1-h_disp_p4) == sign(Pref_HDisp-h_disp_p4)) & (unsigned_hdisp == unique_unsigned_hdisp(i)) );
            pref_dist{i} = spike_rates(pref_trials & select_trials);
            null_trials = ( (sign(h_disp_p1-h_disp_p4) ~= sign(Pref_HDisp-h_disp_p4)) & (unsigned_hdisp == unique_unsigned_hdisp(i)) );
            null_dist{i} = spike_rates(null_trials & select_trials);
            ROC_values(i) = rocN(pref_dist{i}, null_dist{i}, 100);
            N_obs(i) = length(pref_dist{i}) + length(null_dist{i});
            % data for Weibull fit
            fit_data(i, 1) = unique_unsigned_hdisp(i);
            fit_data(i, 2) = ROC_values(i);
            fit_data(i,3) = N_obs(i);
        end
    end

    [neuron_alpha_integ(j) neuron_beta_integ(j)] = weibull_fit(fit_data);
end
%normalize thresholds to the last (1500ms) data
norm_thres = neuron_alpha_integ/neuron_alpha_integ(length(neuron_alpha_integ));


%--------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------
%write out all relevant parameters to a cumulative text file, GCD 8/08/01
outfile = [BASE_PATH 'ProtocolSpecific\Stereoacuity\NeuroPsycho_Curve_summary.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    %fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Nthr\t NSlp\t Mthr\t MSlp\t Pthr\t Pslp\t SurDisp\t Ntrials\t HDisp\t ROChCorr\t');
    fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Nthr\t NSlp\t Mthr\t MSlp\t SurDisp\t Ntrials\t HDisp\t ROChCorr\t');
    fprintf(fid, '\r\n');
    printflag = 0;
end
%buff = sprintf('%s  %4.1f %5.2f %5.3f %4.2f %4.2f %5.2f %7.6f %7.4f %7.6f %7.4f %7.4f %7.4f %5.3f %4d %5.3f %5.3f', FILE, data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED),...
%    data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER), ...
%    neuron_alpha, neuron_beta, monkey_alpha, monkey_beta, P_threshold, P_slope, unique_hdisp_p4(1), (1+ EndTrial - BegTrial), unique_unsigned_hdisp(length(unique_unsigned_hdisp)), ROC_values(length(unique_unsigned_hdisp)) );
buff = sprintf('%s  %4.1f %5.2f %5.3f %4.2f %4.2f %5.2f %7.6f %7.4f %7.6f %7.4f %5.3f %4d %5.3f %5.3f', FILE, data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED),...
    data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER), ...
    neuron_alpha, neuron_beta, monkey_alpha, monkey_beta, unique_hdisp_p4(1), (1+ EndTrial - BegTrial), unique_unsigned_hdisp(length(unique_unsigned_hdisp)), ROC_values(length(unique_unsigned_hdisp)) );
%disp(buff);
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

%-------------------------------------------------------------------------------------------------------------
%now, print out data and fits to a file for external plotting purposes (e.g., in Origin)
i = size(PATH,2) - 1;
while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
    i = i - 1;
end   
PATHOUT = [PATH(1:i) 'Analysis\NeuroPsychoCurves\'];
i = size(FILE,2) - 1;
while FILE(i) ~='.'
    i = i - 1;
end
FILEOUT = [FILE(1:i) 'np_st_curves'];

fileid = [PATHOUT FILEOUT];
fwriteid = eval(['fopen(fileid, ''w'')']);

fprintf(fwriteid,'DispI\tN_fit\tP_fit\tDispR\tN_raw\tP_raw\n');
for i=1:length(fit_x)
    fprintf(fwriteid, '%8.5f\t%6.3f\t%6.3f\t', fit_x(i), neuron_fit_y(i), monkey_fit_y(i));
    if (i <= length(unique_unsigned_hdisp))
        fprintf(fwriteid, '%8.5f\t%6.3f\t%6.3f\n', unique_unsigned_hdisp(i), ROC_values(i), pct_correct(i) );
    else
        fprintf(fwriteid, '\n');
    end
end

fclose(fwriteid);

%------------------------------------------------------------------------
%write out data for integration time analysis TU 08/27/03
outfile2 = [BASE_PATH 'ProtocolSpecific\Stereoacuity\IntegrationTime_summary.dat'];
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

outfile3 = [BASE_PATH 'ProtocolSpecific\Stereoacuity\IntegrationTime_raw.dat'];
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

return;