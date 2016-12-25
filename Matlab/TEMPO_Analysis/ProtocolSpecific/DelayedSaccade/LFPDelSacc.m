%-----------------------------------------------------------------------------------------------------------------------
%-- LFPDelSacc.m -- Computes the response for various LFP bands and plots
%--                 1d or 2d-tuning curves for the various bands.
%--	VR, 1/23/07
%-----------------------------------------------------------------------------------------------------------------------

function LFPDelSacc(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE)

TEMPO_defs;
ProtocolDefs;
Path_Defs;

%get delay period target brightness
targ_dimmer = data.one_time_params(TARG_DIMMER);

%get angle values (deg)
angle = data.dots_params(DOTS_DIREC,:,PATCH1);
angle = squeeze_angle(angle);
unique_angle = munique(angle');

%get eccentricity (radius) values - stored in dummy variable DOTS_AP_XCTR
rad = data.dots_params(DOTS_AP_XCTR,:,PATCH1);
unique_rad = munique(rad');

%make list of trials
trials = 1:length(angle);
select_trials = logical (data.misc_params(OUTCOME,BegTrial:EndTrial) == CORRECT);

timing_offset = 71; %in ms, time between TARGS_ON_CD and detectable target on screen; only use this for aligning to target onset
%get firing rates for delay period (VSTIM_OFF:FP_OFF) and saccade (FP:IN_T1_WIN)
delay_rates = data.spike_rates(SpikeChan,:);
for i = trials
    trialdata = data.event_data(1,:,i);
    if ( sum(trialdata == TARGS_ON_CD) & sum(trialdata == VSTIM_OFF_CD) & ...
            sum(trialdata == FP_OFF_CD) & sum(trialdata == IN_T1_WIN_CD) )
        targ_on(i) = find(data.event_data(1,:,i) == TARGS_ON_CD) + timing_offset;
        delay_start(i) = find(data.event_data(1,:,i) == VSTIM_OFF_CD) + timing_offset;
        fix_off(i) = find(data.event_data(1,:,i) == FP_OFF_CD) + timing_offset;
        in_T1(i) = find(data.event_data(1,:,i) == IN_T1_WIN_CD,1,'last');
        delay_rates(i) = sum(data.spike_data(SpikeChan, delay_start(i):fix_off(i), i)) / length(delay_start(i):fix_off(i)) * 1000;
        sacc_rates(i) = sum(data.spike_data(SpikeChan, fix_off(i):in_T1(i), i)) / length(fix_off(i):in_T1(i)) * 1000;
    else
        select_trials(i) == 0;
    end
end

%note some times
delay_period = mean(fix_off-targ_on);

%now, get the LFP power for all the trials
if (isempty(data.lfp_data)) %in case the lfp data wasn't saved, then don't run anything else here
    disp('No LFP Data saved... quitting.');
    return;
end
for i = trials(select_trials)
    %note that lfp is sampled at half the frequency as spikes, so divide bins by 2
    delay_start(i) = ceil((timing_offset+find(data.event_data(1,:,i) == VSTIM_OFF_CD))/2);
    delay_end(i) = floor((timing_offset+find(data.event_data(1,:,i) == FP_OFF_CD))/2);
    n_samp(i) = delay_end(i)-delay_start(i)+1;

    %now filter lfp signal using Jing Liu's noise canceller.  This filters the 60Hz signal and its first three harmonics
    %(which are visible in some datasets).  The fifth argument (.01) specifies the width of the filter and was determined empirically.
    prefilter_delay_lfp{i} = data.lfp_data(1,delay_start(i):delay_end(i),i);
    filtered_delay_lfp{i} = noise_canceller(prefilter_delay_lfp{i}', 60, [1 1 1], 500, 0.01, 0);
end
%do the following to get the power of lfp between 0 and 200Hz, sampled at 500Hz
nsamp = 400; %only use 400 samples because i only care to filter up to 200Hz
for i = trials
    lfp_delay_powerspect{i} = abs(fft(filtered_delay_lfp{i},nsamp)).^2 ./ nsamp;
end


if (length(unique_angle) == 1)

elseif (length(unique_rad) == 1)
    %first compute DDI, ChiP, and p_anova from spike data
    %first compute the mean and stddev power trials at each angle
    for k = 1:length(unique_angle)
        py = delay_rates( (angle == unique_angle(k)) & select_trials );
        mean_delay_rate(k) = mean(py);
        std_delay_rate(k) = std(py);
    end
    plot_x = angle(select_trials);
    spk_plot_y = delay_rates(select_trials);
    %now shift the data so that the angle with the highest power
    %is in the middle of the distribution
    [m, max_ind] = max(mean_delay_rate);
    adjust_angle = unique_angle - 360.*(unique_angle-unique_angle(max_ind) > 180) + 360.*(unique_angle-unique_angle(max_ind) < -180);
    adjust_plot_x = plot_x - 360.*(plot_x-unique_angle(max_ind) > 180) + 360.*(plot_x-unique_angle(max_ind) < -180);
    %fit the shifted data with a gaussian
    means = [adjust_angle 1e10.*mean_delay_rate'];
    raw = [adjust_plot_x' 1e10.*delay_rates(select_trials)'];
    [pars] = gaussfit(means, raw, 0);   %last arg: allow positive going fit only
    % Do chi-square goodness of fit test
    [spk_chi2, spk_chiP] = Chi2_Test(adjust_plot_x, 1e10.*delay_rates(select_trials), 'gaussfunc', pars, length(pars));
    % Compute direction discrimination index
    spk_delay_DDI = Compute_DDI(angle(select_trials), delay_rates(select_trials));
    % calculate some metrics and stats 
    spk_p_value = spk_anova(spk_plot_y, plot_x, unique_angle);
    

    %now for each possible band, compute the LFPs and use that to compute
    %the same metrics
    lo_bnd = 10:10:190; hi_bnd = 20:10:200; 
    spk_bnd_corr = -1.*ones(length(lo_bnd)); %prepare spkbndcorr matrix
    for i = 1:length(lo_bnd)
        for j = i:length(hi_bnd)
            clear plot_* adjust_* 
            %first define the band, and compute the mean power in each trial
            band = find( (500*(0:nsamp/2)./nsamp >= lo_bnd(i)) & (500*(0:nsamp/2)./nsamp <= hi_bnd(j)) );
            for k = 1:length(lfp_delay_powerspect) %for each trial
                if select_trials(k)
                    lfp_delay_pwr{i,j}(k) = mean(lfp_delay_powerspect{k}(band));
                end
            end
            %    plot_y = lfp_delay_pwr{i,j}(~null_trials & select_trials);

            %first compute the mean and stddev power trials at each angle
            for k = 1:length(unique_angle)
                py = lfp_delay_pwr{i,j}( (angle == unique_angle(k)) & select_trials );
                mean_delay_pwr{i,j}(k) = mean(py);
                std_delay_pwr{i,j}(k) = std(py);
            end
            plot_x = angle(select_trials);
            plot_y = lfp_delay_pwr{i,j}(select_trials);

            %now shift the data so that the angle with the highest power
            %is in the middle of the distribution
            [m, max_ind] = max(mean_delay_pwr{i,j});
            adjust_angle = unique_angle - 360.*(unique_angle-unique_angle(max_ind) > 180) + 360.*(unique_angle-unique_angle(max_ind) < -180);
            adjust_plot_x = plot_x - 360.*(plot_x-unique_angle(max_ind) > 180) + 360.*(plot_x-unique_angle(max_ind) < -180);

            %fit the shifted data with a gaussian
            means = [adjust_angle 1e10.*mean_delay_pwr{i,j}'];
            raw = [adjust_plot_x' 1e10.*lfp_delay_pwr{i,j}(select_trials)'];
            [pars] = gaussfit(means, raw, 0);   %last arg: allow positive going fit only

            % Do chi-square goodness of fit test
            [chi2(i,j), chiP(i,j)] = Chi2_Test(adjust_plot_x, 1e10.*lfp_delay_pwr{i,j}(select_trials), 'gaussfunc', pars, length(pars));
            % Compute direction discrimination index
            lfp_delay_DDI(i,j) = Compute_DDI(angle(select_trials), lfp_delay_pwr{i,j}(select_trials));
            % calculate some metrics and stats 
            pref_angle = pars(3);
            p_value(i,j) = spk_anova(plot_y, plot_x, unique_angle);
            % compute correlation between lfp power at that band and spk rates.
            spk_bnd_corr(i,j) = corr(lfp_delay_pwr{i,j}(select_trials)',spk_plot_y');
            
            %     base_rate = pars(1);
            %     amplitude = pars(2);
            %     max_rate = base_rate + amplitude;
            %     width(i,j) = sqrt(-(log(.5)))*pars(4)*2*sqrt(2);
            %     DSI(i,j) = 1 - (base_rate - null_rate)/(max_rate - null_rate);
            %     %Calculate modulation index using sqrt raw responses and subtracting spontaneous
            %     DMI(i,j) = Compute_ModIndex(plot_x, plot_y, null_resp);
        end
    end


%     %now plot
%     hold on;
%     plot(adjust_plot_x, plot_y, 'r*');
%     handl1 = errorbar(adjust_angle, mean_delay_pwr, std_delay_pwr, 'ro');
%     ax1 = gca;
%     set(ax1,'XColor','k','YColor','r','YAxisLocation','right');
%     xlabel('Angle (deg)'); 	ylabel(ax1, 'LFP Power (a.u.)');
%     axis('tight');
%     title(sprintf('%s: Full Delay Period tuning curve',FILE));

    figh = figure;
    set(figh,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', sprintf('%s: LFP at different bands',FILE));
    colormap(gray);
    
    subplot(211); %first plot out the delay DDI
    x = repmat(lfp_delay_DDI, [1 1 3]); %to use image in gray scale in need to replicate the matrix
    x(isnan(x)) = 0;
    image([20:10:200], [10:10:190], x);
    cb = colorbar('Location','EastOutside');
    yl = get(cb,'YLim');
    set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
    xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
    title('Direction Discrimination Index');
    
    subplot(212); %now plot the SU-MU correlation
    x = repmat(spk_bnd_corr./2+.5, [1 1 3]); %to use image in gray scale, replicate the matrix
    x(isnan(x)) = 0;
    image([20:10:200], [10:10:190], x);
    cb = colorbar('Location','EastOutside');
    yl = get(cb,'YLim');
    set(cb,'YTick',[min(yl):range(yl)/4:max(yl)],'YTickLabel',{'-1.0' '-0.5' '0.0' '0.5' '1.0'});
    xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
    title('SU-LFP correlation');
    
    savefig = 1;
    if (savefig)
        savename = sprintf('Z:\\Data\\Tempo\\Baskin\\Analysis\\LFP_tuning_DelSacc\\%s_LIP_lfp',FILE(1:8));
        saveas(figh, savename);
    end
        


    output = 1;
    if output %save DirDI, anova p_val, chi2 p_val, and filenames to cumulative matfiles

        outfile = [BASE_PATH 'ProtocolSpecific\DelayedSaccade\LFP_DelSacc_summary.mat'];
        if (exist(outfile, 'file') == 0)
            cum_DirDI{1} = lfp_delay_DDI;
            cum_p_anova{1} = p_value;
            cum_chiP{1} = chiP;
            cum_spk_bnd_corr{1} = spk_bnd_corr;
            cum_spk_metrics = [spk_delay_DDI spk_chiP spk_p_value];
            files{1} = FILE;
            save('-v6', outfile, 'cum_DirDI', 'cum_p_anova', 'cum_chiP', 'cum_spk_bnd_corr', 'cum_spk_metrics', 'files');
        else
            load(outfile);
            cum_DirDI{length(cum_DirDI)+1} = lfp_delay_DDI;
            cum_p_anova{length(cum_p_anova)+1} = p_value;
            cum_chiP{length(cum_chiP)+1} = chiP;
            cum_spk_bnd_corr{length(cum_spk_bnd_corr)+1} = spk_bnd_corr;
            cum_spk_metrics = [cum_spk_metrics; spk_delay_DDI spk_chiP spk_p_value];
            files{length(files)+1} = FILE;
            save('-v6', outfile, 'cum_DirDI', 'cum_p_anova', 'cum_chiP', 'cum_spk_bnd_corr', 'cum_spk_metrics', 'files');
        end

    end

else %2d gaussian fit

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% After running all the individual trials and saving the data to
% LFP_indices.mat, run the following code to look at population data.
load LFP_DelSacc_summary.mat
n = length(cum_DirDI);
DirDI = zeros(19);
chiP = zeros(19);
p_anova = zeros(19);
spk_bnd_corr = zeros(19);
for i = 1:n
    DirDI = DirDI + cum_DirDI{i}./n;
    chiP = chiP + (cum_chiP{i}>0.05)./n;
    p_anova = p_anova + (cum_p_anova{i}<.05)./n;
    spk_bnd_corr = spk_bnd_corr + cum_spk_bnd_corr{i}./n;
end

figure
colormap(gray)
x = repmat(DirDI,[1 1 3]);
x(isnan(x)) = 0;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title(sprintf('Dir Discrim Index (spk mean: %3.2f)',mean(cum_spk_metrics(:,1))));

figure
colormap(gray)
x = repmat(chiP,[1 1 3]);
x(isnan(x)) = 0;
x((x>1) & (x<1.000001)) = 1;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title(sprintf('Frxn Sites w/ChiP>0.05 (spk mean: %3.2f)',mean(cum_spk_metrics(:,2)>.05)));

figure
colormap(gray)
x = repmat(p_anova,[1 1 3]);
x(isnan(x)) = 0;
x((x>1) & (x<1.000001)) = 1;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title(sprintf('Frxn Sites w/sig tuning (ANOVA) (spk mean: %3.2f)',mean(cum_spk_metrics(:,3)<.05)));

figure
colormap(gray)
x = repmat(spk_bnd_corr./2+.5,[1 1 3]);
x(isnan(x)) = 0;
x((x<0) & (x>-0.000001)) = 0;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/4:max(yl)],'YTickLabel',{'-1.0' '-0.5' '0.0' '0.5' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('Power-SpkRate Correlation');
