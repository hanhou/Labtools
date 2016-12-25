%-----------------------------------------------------------------------------------------------------------------------
%-- GaborFit_RelDisp.m -- Fits disparity tuning curves with full Gabors, for each individual surround disparity.
%-- NOTE: responses are square-rooted in the error functions, to help homogenize variance (a la Cumming).  
%-- A constrained version of Levenberg-Marquardt optimization is used, as implemented using 'fmincon'.
%--	TU 07/21/01
%-----------------------------------------------------------------------------------------------------------------------
function GaborFit_TransRelDisp(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

%determine which parameters are going to be shared. shared = 1, not shared = 0.
Base_Rate = 1;
Amplitude = 1;
Gauss_Ctr = 0;
Gauss_SD = 1;
Sine_Freq = 1;
Sine_Phase = 1;

%determine whether shift ratios will be calculated for all combination. All combination = 1, relative to zero-disparity surround = 0.
all_combination = 2;

% not implemented yet to select output
output = 0;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp_ctr = data.dots_params(DOTS_HDISP,:,PATCH1);
hor_disp_sur = data.dots_params(DOTS_HDISP,:,PATCH4);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp_ctr == data.one_time_params(NULL_VALUE)) & (hor_disp_sur == data.one_time_params(NULL_VALUE)) );

% Ctr_Control_trials = no dots in ctr, surround varies
ctr_control_trials = logical(hor_disp_ctr == data.one_time_params(PATCH_OFF));
   
% Sur_Control_trials = no dots in surround, ctr disp varies
sur_control_trials = logical(hor_disp_sur == data.one_time_params(PATCH_OFF));
   
control_trials = (ctr_control_trials | sur_control_trials);

unique_hor_disp_sur = munique(hor_disp_sur(~null_trials & ~sur_control_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp_ctr);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));

%get the average eye positions to calculate vergence
if (data.eye_calib_done == 1)
    Leyex_positions = data.eye_positions_calibrated(1, :);
    Leyey_positions = data.eye_positions_calibrated(2, :);
    Reyex_positions = data.eye_positions_calibrated(3, :);
    Reyey_positions = data.eye_positions_calibrated(4, :);
     
    vergence_h = Leyex_positions - Reyex_positions;
    vergence_v = Leyey_positions - Reyey_positions;
else     
    Leyex_positions = data.eye_positions(1, :);
    Leyey_positions = data.eye_positions(2, :);
    Reyex_positions = data.eye_positions(3, :);
    Reyey_positions = data.eye_positions(4, :);
    
    vergence_h = Leyex_positions - Reyex_positions;
    vergence_v = Leyey_positions - Reyey_positions;
end
    
%check if there is more than two surround disparities
if (length(unique_hor_disp_sur)<2) 
    disp('This analysis only for runs with more than two surround disparities interleaved');
    return;
end

% First fit separate gabors for each surround disparity.
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'HDisp Tuning Curve w/ Best Fit');

for i=1:length(unique_hor_disp_sur)	%for each different surround disparity, plot a separate disparity tuning curve.
    
    hor_disp_sur_select = logical( (hor_disp_sur == unique_hor_disp_sur(i)) );
    
    plot_x = hor_disp_ctr(hor_disp_sur_select & ~null_trials & ~control_trials & select_trials);
    plot_y = spike_rates(hor_disp_sur_select & ~null_trials & ~control_trials & select_trials);  
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    p_value(i) = spk_anova(plot_y, plot_x, px);
    
    means{i} = [px py];
    raw{i} = [plot_x' plot_y'];
    sem{i} = [perr];
    verg{i} = vergence_h(hor_disp_sur_select & ~null_trials & ~control_trials & select_trials); 

    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    
    %fit the GABOR finction
    [pars{i},freq(i)] = gaborfit(means{i},raw{i},fixed_param_flags,fixed_param_values);

    subplot(5,2,2*i);
    
    %plot the raw data
    PlotRawData(plot_x', plot_y', symbols{i}, null_rate);
  
    x_interp = (px(1)): .01 : (px(length(px)));
    y_gabor = gaborfunc(x_interp, pars{i});
    y_gabor(y_gabor < 0) = 0;
    y_gauss =  pars{i}(1) + pars{i}(2)*exp(-0.5*((x_interp - pars{i}(3))/ pars{i}(4)).^2);
    y_sine =  pars{i}(1) + 0.5*pars{i}(2)*cos(2*pi*pars{i}(5)*(x_interp - pars{i}(3))+pars{i}(6) );
    
    %Compute peak and trough disparity
    [peak_rate(i) peak_index] = max(y_gabor);
    peak_disp(i) = x_interp(peak_index);
    [trough_rate(i) trough_index] = min(y_gabor);
    trough_disp(i) = x_interp(trough_index);
    
    % Now plot fitted curves; square the fitted values to counteract sqrt() above
    hold on;
    plot(x_interp, y_gabor, 'k-');
    plot(x_interp, y_gauss, 'k--');
    plot(x_interp, y_sine, 'k:');
    Gabor_SSE(i) = gaborerr(pars{i});
    buff = sprintf('Full Gabor fit (sse = %8.5f)', gaborerr(pars{i}) );
    title(buff);
    hold off;   
    
    y_fit = gaborfunc(px, pars{i});
    y_fit(y_fit < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit),1) y_fit];
    [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
    
    y_fit_raw = gaborfunc(plot_x', pars{i});
    y_fit_raw(y_fit_raw < 0) = 0;
    y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
    [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);
end

DDI_index = zeros(1,length(unique_hor_disp_sur));  % store DDI values for each surround disparity

outstr2 = [0];
for i=1:length(unique_hor_disp_sur)	%for each different surround disparity, plot some information.
    
    %print fit parameters to the screen
    subplot(5,2,2*i-1);
    axis([0 100 0 100]);    axis('off');
    xpos = -10; ypos = 100;
    font_size = 8;  bump_size = 10;

    if(i == 1)
        temp = strcat(PATH, FILE);
        temp(temp == '\') = '/';
        % this prevents a stupid error from appearing on the screen
        line = sprintf('File: %s', temp);
        text(xpos-15,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    end
    line = sprintf('Pedestal: %6.2f degree', unique_hor_disp_sur(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Base Rate: %8.4f', pars{i}(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Amplitude: %8.4f', pars{i}(2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss Ctr: %8.4f', pars{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss  SD: %8.4f', pars{i}(4));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sine Freq: %8.4f', pars{i}(5));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sin Phase: %8.4f', pars{i}(6));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('FT Freq: %8.4f', freq(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    
    %calculate shift ratio from individual gabor fits
    %shift ratio = (peak or trough difference - vergence difference)/(surround disparity difference - vergence difference)
    if(all_combination == 1) %calculate shift ratio for all combination
        if(p_value(i) < 0.01) %calculate shift ratio only if tuning is significant
            line = sprintf('Shift Ratio:');
            for j = i+1:length(unique_hor_disp_sur) %calculate shift ratio for each other surround disparity
                if(p_value(j) < 0.01)
                    %calculate shift ratio from peaks or troughs of individual gabor fits
                    shift_ratio_temp(1) = ((trough_disp(i)-trough_disp(j))-(mean(verg{i})-mean(verg{j})))/((unique_hor_disp_sur(i)-unique_hor_disp_sur(j))-((mean(verg{i})-mean(verg{j}))));
                    shift_ratio_temp(2) = ((peak_disp(i)-peak_disp(j))-(mean(verg{i})-mean(verg{j})))/((unique_hor_disp_sur(i)-unique_hor_disp_sur(j))-((mean(verg{i})-mean(verg{j}))));
                    if(logical((abs(pars{i}(6))>0.5*pi) & (abs(pars{i}(6))<1.5*pi)) == logical((abs(pars{i}(6))>0.5*pi) & (abs(pars{i}(6))<1.5*pi))) %use smaller shift ratio if phase for both curves do not agree
                        [abs_shift_ratio shift_ratio_index] = min(abs(shift_ratio_temp));
                        shift_ratio = shift_ratio_temp(shift_ratio_index);
                    elseif((abs(pars{i}(6))>0.5*pi) & (abs(pars{i}(6))<1.5*pi)) %use shift ratio calculated from trough if phase for both curves are near 180 deg
                        shift_ratio = shift_ratio_temp(1);   
                    else %use shift ratio calculated from peak if phase for both curves are near 0 deg
                        shift_ratio = shift_ratio_temp(2);
                    end    
                    
                    %significance test for tuning shape. fit two tuning curves with a single fit, allowing base_rate and amplitude to vary; 
                    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
                    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
                    shared_param_flags = [0 0 1 1 1 1]; %share everything except base_rate and amplitude
                    means_signif{1} = means{i};
                    means_signif{2} = means{j};
                    raw_signif{1} = raw{i};
                    raw_signif{2} = raw{j};
                    [pars_shared_signif,freq_shared_signif] = multi_gaborfit(means_signif,raw_signif,fixed_param_flags,fixed_param_values,shared_param_flags);
                    
                    %do Sequential F-test to test if individual gabor fits are better than paired fits with same shape gabor
                    hor_disp_sur_select = ( logical( (hor_disp_sur == unique_hor_disp_sur(i)) ) | logical( (hor_disp_sur == unique_hor_disp_sur(j)) ) );
                    n_sur_disp = 2;
                    independent_err = Gabor_SSE(i) + Gabor_SSE(j);
                    shared_err = multi_gaborerr(pars_shared_signif);
                    n_indep_params = length(shared_param_flags)*n_sur_disp;
                    n_shared_params = n_indep_params-(n_sur_disp-1)*sum(shared_param_flags);
                    Npts = length(spike_rates(~null_trials & ~control_trials & select_trials & hor_disp_sur_select));
                    Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
                    Pseq = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );
                    
                    %Now, calculate shift ratio from paired gabor fits allowing the gaussCtr to vary, ala Thomas
                    %Try varying the amplitude and baserate as well 12/14/01 TU
                    shared_param_flags = [0 0 0 1 1 1]; %share everything except gaussCtr, baserate and amplitude
                    [pars_shared_signif,freq_shared_signif] = multi_gaborfit(means_signif,raw_signif,fixed_param_flags,fixed_param_values,shared_param_flags);
                    shared_err_gaussCtr = multi_gaborerr(pars_shared_signif);
                    sum_shared_params_gaussCtr = sum(shared_param_flags);
                    shift_ratio_paired = ((pars_shared_signif(3)-pars_shared_signif(9))-(mean(verg{i})-mean(verg{j})))/((unique_hor_disp_sur(i)-unique_hor_disp_sur(j))-((mean(verg{i})-mean(verg{j}))));
                    
                    %do sequential F-test to test if shift ratio is significant
                    shared_param_flags = [0 0 1 1 1 1]; %share everything except baserate and amplitude
                    [pars_shared_signif,freq_shared_signif] = multi_gaborfit(means_signif,raw_signif,fixed_param_flags,fixed_param_values,shared_param_flags);
                    shared_err_all = multi_gaborerr(pars_shared_signif);
                    sum_shared_params_all = sum(shared_param_flags);                   
                    independent_err = shared_err_gaussCtr;
                    shared_err = shared_err_all;
                    n_indep_params = length(shared_param_flags)*n_sur_disp-(n_sur_disp-1)*sum_shared_params_gaussCtr;
                    n_shared_params = length(shared_param_flags)*n_sur_disp-(n_sur_disp-1)*sum_shared_params_all;
                    Npts = length(spike_rates(~null_trials & ~control_trials & select_trials & hor_disp_sur_select));
                    Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
                    Pseq_paired = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );
                    
                    %do additional sequential F-test to test if individual gabor fits are better than paired gabor fits allowing the gaussCtr to vary
                    independent_err = Gabor_SSE(i) + Gabor_SSE(j);
                    shared_err = shared_err_gaussCtr;
                    n_indep_params = length(shared_param_flags)*n_sur_disp;
                    n_shared_params = n_indep_params-(n_sur_disp-1)*sum_shared_params_gaussCtr;
                    Npts = length(spike_rates(~null_trials & ~control_trials & select_trials & hor_disp_sur_select));
                    Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
                    Pseq_gaussCtr = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );
                    
                    line = sprintf('%s %8.4f', line, shift_ratio);
                    outstr2 = sprintf('%s %s %8.4f %8.4f %8.4f %8.4f %8.4f\r\n', outstr2, FILE, shift_ratio, Pseq, shift_ratio_paired, Pseq_paired, Pseq_gaussCtr);
                else
                    line = sprintf('%s n/a', line);
                end
            end
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        else
            line = sprintf('Shift Ratio: n/a');
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        end
        
    elseif all_combination == 0 %calculate shift ratio relative to zero-disparity surround
        disp_index = find(unique_hor_disp_sur == 0.0);
        if(p_value(disp_index) < 0.01) %calculate shift ratio only if tuning is significant
            if(i == disp_index) %cannot calculate shift ratio for zero-degree surround
                line = sprintf('Shift Ratio: n/a');
                text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            else
                if(p_value(i) < 0.01)
                    if(abs(pars{disp_index}(6))>0.5*pi & abs(pars{disp_index}(6))<1.5*pi) %use trough disparity to calculate shift ratio if phase is near 180 deg
                        shift_ratio = ((trough_disp(i)-trough_disp(disp_index))-(mean(verg{i})-mean(verg{disp_index})))/((unique_hor_disp_sur(i)-unique_hor_disp_sur(disp_index))-((mean(verg{i})-mean(verg{disp_index}))));
                    else
                        shift_ratio = ((peak_disp(i)-peak_disp(disp_index))-(mean(verg{i})-mean(verg{disp_index})))/((unique_hor_disp_sur(i)-unique_hor_disp_sur(disp_index))-((mean(verg{i})-mean(verg{disp_index}))));
                    end    
                    line = sprintf('Shift Ratio: %8.4f', shift_ratio);
                    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size; 
                    outstr2 = sprintf('%s %8.4f ', outstr2, shift_ratio);
                else
                    line = sprintf('Shift Ratio: n/a');
                    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
                end
            end
        else
            line = sprintf('Shift Ratio: n/a');
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        end
    elseif all_combination == 2
        [dummy, choose_DDI] = max(DDI_index);
        if i ~= choose_DDI
            line = sprintf('Shift Ratio:');
            j = choose_DDI;
                    if(p_value(i) < 0.01 & p_value(j) <.01) %calculate shift ratio only if tuning is significant
                        
                    %significance test for tuning shape. fit two tuning curves with a single fit, allowing base_rate and amplitude to vary; 
                    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
                    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
                    %shared_param_flags = [0 0 1 1 1 1]; %share everything except base_rate and amplitude
                    means_signif{1} = means{i};
                    means_signif{2} = means{j};
                    raw_signif{1} = raw{i};
                    raw_signif{2} = raw{j};
                    %[pars_shared_signif,freq_shared_signif] = multi_gaborfit(means_signif,raw_signif,fixed_param_flags,fixed_param_values,shared_param_flags);
                    
                    %do Sequential F-test to test if individual gabor fits are better than paired fits with same shape gabor
                    hor_disp_sur_select = ( logical( (hor_disp_sur == unique_hor_disp_sur(i)) ) | logical( (hor_disp_sur == unique_hor_disp_sur(j)) ) );
                    n_sur_disp = 2;
                    %independent_err = Gabor_SSE(i) + Gabor_SSE(j);
                    %shared_err = multi_gaborerr(pars_shared_signif);
                    %n_indep_params = length(shared_param_flags)*n_sur_disp;
                    %n_shared_params = n_indep_params-(n_sur_disp-1)*sum(shared_param_flags);
                    %Npts = length(spike_rates(~null_trials & ~control_trials & select_trials & hor_disp_sur_select));
                    %Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
                    %Pseq = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );
                    
                    %Now, calculate shift ratio from paired gabor fits allowing the gaussCtr to vary, ala Thomas
                    %Try varying the amplitude and baserate as well 12/14/01 TU
                    shared_param_flags = [0 0 0 1 1 1]; %share everything except gaussCtr, baserate and amplitude
                    [pars_shared_signif,freq_shared_signif] = multi_gaborfit(means_signif,raw_signif,fixed_param_flags,fixed_param_values,shared_param_flags);
                    shared_err_gaussCtr = multi_gaborerr(pars_shared_signif);
                    sum_shared_params_gaussCtr = sum(shared_param_flags);
                    shift_ratio_paired = ((pars_shared_signif(3)-pars_shared_signif(9))-(mean(verg{i})-mean(verg{j})))/((unique_hor_disp_sur(i)-unique_hor_disp_sur(j))-((mean(verg{i})-mean(verg{j}))));
                    
                    %do sequential F-test to test if shift ratio is significant
                    shared_param_flags = [0 0 1 1 1 1]; %share everything except baserate and amplitude
                    [pars_shared_signif,freq_shared_signif] = multi_gaborfit(means_signif,raw_signif,fixed_param_flags,fixed_param_values,shared_param_flags);
                    shared_err_all = multi_gaborerr(pars_shared_signif);
                    sum_shared_params_all = sum(shared_param_flags);                   
                    independent_err = shared_err_gaussCtr;
                    shared_err = shared_err_all;
                    n_indep_params = length(shared_param_flags)*n_sur_disp-(n_sur_disp-1)*sum_shared_params_gaussCtr;
                    n_shared_params = length(shared_param_flags)*n_sur_disp-(n_sur_disp-1)*sum_shared_params_all;
                    Npts = length(spike_rates(~null_trials & ~control_trials & select_trials & hor_disp_sur_select));
                    Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
                    Pseq_paired = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );
                    
                    %do additional sequential F-test to test if individual gabor fits are better than paired gabor fits allowing the gaussCtr to vary
                    independent_err = Gabor_SSE(i) + Gabor_SSE(j);
                    shared_err = shared_err_gaussCtr;
                    n_indep_params = length(shared_param_flags)*n_sur_disp;
                    n_shared_params = n_indep_params-(n_sur_disp-1)*sum_shared_params_gaussCtr;
                    Npts = length(spike_rates(~null_trials & ~control_trials & select_trials & hor_disp_sur_select));
                    Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
                    Pseq_gaussCtr = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );
                    
                    line = sprintf('%s %8.4f', line, shift_ratio_paired);
                    outstr2 = sprintf('%s %s %8.4f %8.4f %8.4f\r\n', outstr2, FILE,  shift_ratio_paired, Pseq_paired, Pseq_gaussCtr);
                    else
                        line = sprintf('%s n/a', line);
                    end
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        else
            line = sprintf('Shift Ratio: n/a');
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        end
    end
end

%Next, fit multiple gabors with shared parameters.
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'HDisp Tuning Curve w/ Best Fit');
    
fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
shared_param_flags = [Base_Rate Amplitude Gauss_Ctr Gauss_SD Sine_Freq Sine_Phase]; %shared parameters are defined above

%fit gabor functions to multiple curves with shared parameters
[pars_shared,freq_shared] = multi_gaborfit(means,raw,fixed_param_flags,fixed_param_values,shared_param_flags);

%re-structure the parameters
k = 1;
params{1} = pars_shared(1:6);
for i=2:length(unique_hor_disp_sur)
    params{i} = pars_shared(1:6);
    for j=1:6
        if(shared_param_flags(j) == 0)
            params{i}(j) = pars_shared(6 + k);
            k = k + 1;
        end
    end
end

for i=1:length(unique_hor_disp_sur)	%for each different surround disparity, plot a separate disparity tuning curve.
    
    subplot(5,2,2*i);
    
    %plot the raw data
    PlotRawData(raw{i}(:,1), raw{i}(:,2), symbols{i}, null_rate);

    x_interp = (means{i}(1,1)): .01 : (means{i}(length(means{i}),1));
    y_gabor = gaborfunc(x_interp, params{i});
    y_gabor(y_gabor < 0) = 0;
    y_gauss =  params{i}(1) + params{i}(2)*exp(-0.5*((x_interp - params{i}(3))/ params{i}(4)).^2);
    y_sine =  params{i}(1) + 0.5*params{i}(2)*cos(2*pi*params{i}(5)*(x_interp - params{i}(3))+params{i}(6) );
    
    % Now plot fitted curves; square the fitted values to counteract sqrt() above
    hold on;
    plot(x_interp, y_gabor, 'k-');
    plot(x_interp, y_gauss, 'k--');
    plot(x_interp, y_sine, 'k:');
    hold off;   
    
    y_fit = gaborfunc(means{i}(:,1), params{i});
    y_fit(y_fit < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit),1) y_fit];
    [b, bint, r, rint, stats1{i}] = regress(means{i}(:,2), y_fit);
    
    y_fit_raw = gaborfunc(raw{i}(:,1), params{i});
    y_fit_raw(y_fit_raw < 0) = 0;
    y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
    [b, bint, r, rint, stats2{i}] = regress(raw{i}(:,2), y_fit_raw);
    
    %print fit parameters to the screen
    subplot(5,2,2*i-1);
    axis([0 100 0 100]);    axis('off');
    xpos = -10; ypos = 100;
    font_size = 8;  bump_size = 10;

    if(i == 1)
        temp = strcat(PATH, FILE);
        temp(temp == '\') = '/';
        % this prevents a stupid error from appearing on the screen
        line = sprintf('File: %s', temp);
        text(xpos-15,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    end
    line = sprintf('Pedestal: %6.2f degree', unique_hor_disp_sur(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Base Rate: %8.4f', params{i}(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Amplitude: %8.4f', params{i}(2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss Ctr: %8.4f', params{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss  SD: %8.4f', params{i}(4));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sine Freq: %8.4f', params{i}(5));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sin Phase: %8.4f', params{i}(6));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('FT Freq: %8.4f', freq_shared(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
 
    %calculate shift ratio
    disp_index = find(unique_hor_disp_sur == 0.0);
    if(p_value(disp_index) < 0.05)
        if(i == disp_index)
            line = sprintf('Shift Ratio: n/a');
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        else
            if(p_value(i) < 0.05)
                shift_ratio = ((params{i}(3)-params{disp_index}(3))-(mean(verg{i})-mean(verg{disp_index})))/((unique_hor_disp_sur(i)-unique_hor_disp_sur(disp_index))-((mean(verg{i})-mean(verg{disp_index}))));
                line = sprintf('Shift Ratio: %8.4f', shift_ratio);
                text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size; 
%                outstr3 = sprintf('%s %8.4f ', outstr3, shift_ratio);
            else
                line = sprintf('Shift Ratio: n/a');
                text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            end
        end
    else
        line = sprintf('Shift Ratio: n/a');
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    end    
end

%do Sequential F-test
independent_err = sum(Gabor_SSE);
shared_err = multi_gaborerr(pars_shared);
n_indep_params = length(shared_param_flags)*length(unique_hor_disp_sur);
n_shared_params = n_indep_params-(length(unique_hor_disp_sur)-1)*sum(shared_param_flags);
Npts = length(spike_rates(~null_trials & ~control_trials & select_trials));
Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
Pseq = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );

ypos = ypos - bump_size;
line = sprintf('Sequential F-test: F = %7.3f, P = %7.5f', Fseq, Pseq);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    
%save some values to file
outstr1 = sprintf('%s %8.4f %8.6f', FILE, Fseq, Pseq);

%----------------------------------------------------------------------------------------------------------------------------------------------------------
%write out data in form suitable for plotting tuning curve with Origin.
null_x = [min(means{1}(:,1)) max(means{1}(:,1))];
null_y = [null_rate null_rate];

i = size(PATH,2) - 1;
while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
    i = i - 1;
end   
PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];

i = size(FILE,2) - 1;
while FILE(i) ~='.'
    i = i - 1;
end
FILEOUT = [FILE(1:i) 'reldisp_curv'];
fileid = [PATHOUT FILEOUT];
proffid = eval(['fopen(fileid, ''w'')']);
 
fprintf(proffid,'FitX ');
for i=1:length(unique_hor_disp_sur)
    fprintf(proffid,'GaborFit%d ', i);
end    
fprintf(proffid,'HDisp ');
for i=1:length(unique_hor_disp_sur)
    fprintf(proffid,'AvgResp%d StdErr%d ', i, i);
end    
fprintf(proffid,'HDisp2 Spon\n');

for i=1:length(x_interp)
    fprintf(proffid,'%6.2f ', x_interp(i));
    for j=1:length(unique_hor_disp_sur)
        y_gabor = gaborfunc(x_interp, pars{j});
        y_gabor(y_gabor < 0) = 0;
        fprintf(proffid,'%6.2f ', y_gabor(i));
    end    
    if (i<=length(sem{1}))    
        fprintf(proffid,'%6.2f ', means{1}(i,1));
        for j=1:length(unique_hor_disp_sur)
            fprintf(proffid,'%6.2f %6.2f ', means{j}(i,2), sem{j}(i));
        end
        if (i == 1)
            fprintf(proffid,'%6.2f %6.2f\n',null_x(i),null_y(i));
        elseif (i==2)		
            fprintf(proffid,'%6.2f %6.2f\n',null_x(i),null_y(i));
        else
            fprintf(proffid,'\n');
        end
    else
        fprintf(proffid,'\n');
    end
end
    
fclose(proffid);
    

%also write out data for summary.
outfile1 = [BASE_PATH 'ProtocolSpecific\TransRelativeDisparity\Seq_Ftest_GCtr.dat'];
outfile2 = [BASE_PATH 'ProtocolSpecific\TransRelativeDisparity\Shift_ratio_summary.dat'];

%printflag = 0;
%if (exist(outfile1, 'file') == 0)    %file does not yet exist
%    printflag = 1;
%end
%fid = fopen(outfile1, 'a');
%if (printflag)
%    fprintf(fid, 'File Fseq Pseq');
%    fprintf(fid, '\r\n');
%end
%fprintf(fid, '%s', outstr1);
%fprintf(fid, '\r\n');
%fclose(fid);


printflag = 0;
if (exist(outfile2, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile2, 'a');
if (printflag)
    fprintf(fid, 'File ShiftRatioPaired PPaired PGaussCtr');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', outstr2);
%fprintf(fid, '\r\n');
fclose(fid);

return;


%----------------------------------------
function    PlotRawData(x, y, symb1, null_resp)

plot(x, y, symb1);

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(x) max(x)];
null_y = [null_resp null_resp];
hold on;
plot(null_x, null_y, 'k--');
hold off;

return;