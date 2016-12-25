%-----------------------------------------------------------------------------------------------------------------------
%-- LFPDirectionTuningCurve.m -- Computes the Direction Tuning Index for various bands of the power spectrum.
%--	VR, 5/1/06
%-----------------------------------------------------------------------------------------------------------------------
function LFPDirectionTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

tic 

TEMPO_defs;
Path_Defs;

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

%now, get the LFP power for all the trials
if (isempty(data.lfp_data)) %in case the lfp data wasn't saved, then don't run anything else here
    disp('No LFP Data saved... quitting.');
    return;
end
for i = 1:length(direction)
    %note that lfp is sampled at half the frequency as spikes, so divide bins by 2
    start_stim(i) = ceil(find(data.event_data(1,:,i) == VSTIM_ON_CD)/2);
    end_stim(i) = floor(find(data.event_data(1,:,i) == VSTIM_OFF_CD)/2);
    n_samp(i) = end_stim(i)-start_stim(i)+1;
    
    %now filter lfp signal using Jing Liu's noise canceller.  This filters the 60Hz signal and its first three harmonics 
    %(which are visible in some datasets).  The fifth argument (.01) specifies the width of the filter and was determined empirically.
    prefilter_stim_lfp{i} = data.lfp_data(1,start_stim(i):end_stim(i),i);
    filtered_stim_lfp{i} = noise_canceller(prefilter_stim_lfp{i}', 60, [1 1 1], 500, 0.01, 0);
end


%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (direction == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(direction);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%do the following to get the power of lfp between 0 and 200Hz, sampled at 500Hz
nsamp = 400; %only use 400 samples because i only care to filter up to 200Hz
for i = 1:length(direction)
    lfp_stim_powerspect{i} = abs(fft(filtered_stim_lfp{i},nsamp)).^2 ./ nsamp;
end


plot_x = direction(~null_trials & select_trials);

lo_bnd = 10:10:190;
% lo_bnd = 50; %FOR TESTING PURPOSES... DELETE THIS LINE
for i = 1:length(lo_bnd)
    hi_bnd = 20:10:200;
    %     hi_bnd = 150; %FOR TESTING PURPOSES.... DELETE THIS LINE
    for j = 1:length(hi_bnd)
        if (lo_bnd(i) < hi_bnd(j))
            band = find( (500*(0:nsamp/2)./nsamp >= lo_bnd(i)) & (500*(0:nsamp/2)./nsamp <= hi_bnd(j)) );
            for k = 1:length(lfp_stim_powerspect)
                lfp_pwr{i,j}(k) = mean(lfp_stim_powerspect{k}(band));
            end
            plot_y = lfp_pwr{i,j}(~null_trials & select_trials);

            %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
            [px, py, perr] = PlotTuningCurve(plot_x', plot_y', 'ko', 'k-', 1, 0);	%note: last arg=0 means just get output, no plot

            %keep a copy of the original data before shifting below
            px_orig = px;
            py_orig = py;

            unique_dirs = px; % store category groups for ANOVAs

            % Compute a direction discrimination index analogous to the DDI
            [DirDI(i,j), var_term(i,j)] = Compute_DDI(plot_x, plot_y);

            %now, fit the data with a Gaussian curve and plot this as well
            means = [px py];
            raw = [plot_x' plot_y'];
            [pars] = gaussfit(means, raw, 0);   %last arg: allow positive going fit only
            x_interp = (px(1)): 0.5 : (px(length(px)));
            y_interp = gaussfunc(x_interp, pars);

            %now, get the LFP power for NULL condition trials and add spontaneous rate to plot
            null_x = [min(px) max(px)];
            null_resp = lfp_pwr{i,j}(null_trials & select_trials);
            null_rate = mean(null_resp);
            null_y = [null_rate null_rate];

            %Compute R^2 of the fit for both means and raw values
            y_fit = gaussfunc(px, pars);
            y_fit(y_fit < 0) = 0;
            %add a column of ones to yfit to make regress happy
            y_fit = [ones(length(y_fit),1) y_fit];
            [b, bint, r{i,j}, rint, stats1{i,j}] = regress(py, y_fit);

            y_fit_raw = gaussfunc(plot_x', pars);
            y_fit_raw(y_fit_raw < 0) = 0;
            y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
            [b, bint, r{i,j}, rint, stats2{i,j}] = regress(plot_y', y_fit_raw);

            % Do chi-square goodness of fit test
            [chi2(i,j), chiP(i,j)] = Chi2_Test(plot_x, plot_y, 'gaussfunc', pars, length(pars));


            % calculate some metrics and stats then print them in plot
            pref_dir = pars(3);
            p_value(i,j) = spk_anova(plot_y, plot_x, unique_dirs);
            base_rate = pars(1);
            amplitude = pars(2);
            max_rate = base_rate + amplitude;
            width = sqrt(-(log(.5)))*pars(4)*2*sqrt(2);
            DSI(i,j) = 1 - (base_rate - null_rate)/(max_rate - null_rate);

            %Calculate modulation index using sqrt raw responses and subtracting spontaneous
            DMI(i,j) = Compute_ModIndex(plot_x, plot_y, null_resp);

            %         keyboard

            if (lo_bnd==50) & (hi_bnd==150)
                figure;
                set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', sprintf('%s: Direction Tuning Curve',FILE));
                subplot(2, 1, 2);

                %now, shift the px, py, and perr vectors such that the peak of tuning curve is in middle of axis range
                % now, we need to shift these vectors so that the peak response is always in middle of vector
                ctr_indx = round(length(px)/2 - rem(length(px),2)) + 1;
                [max_val max_indx] = max(py);
                shift = max_indx - ctr_indx;
                if (shift > 0)
                    px = [ px(shift+1 : length(px)) ; px(1 : shift)+360];
                    px_orig = [ px_orig(shift+1 : length(px_orig)) ; px_orig(1 : shift)];
                    py = [ py(shift+1 : length(py)) ; py(1 : shift)];
                    perr = [ perr(shift+1 : length(perr)) ; perr(1 : shift)];
                end
                if (shift < 0)
                    px = [ px(length(px)+shift+1 : length(px))-360 ; px(1 : length(px)+shift)];
                    px_orig = [ px_orig(length(px_orig)+shift+1 : length(px_orig)) ; px_orig(1 : length(px_orig)+shift)];
                    py = [ py(length(py)+shift+1 : length(py)) ;  py(1 : length(py)+shift)];
                    perr = [ perr(length(perr)+shift+1 : length(perr)) ;  perr(1 : length(perr)+shift)];
                end

                %now apply the shift to the raw data (plot_x and plot_y)
                for k=1:length(px_orig)
                    select = logical(plot_x == px_orig(k));
                    plot_x(select) = px(k);
                end

                % Since direction is circular, duplicate  the lowest value (px) as px+ 360, and change spike arrays accordingly
                px = [px; px(1)+360];
                py = [py; py(1)];
                perr = [perr; perr(1)];

                hold on;
                %plot the data, after shifting as necessary above
                plot(plot_x, plot_y, 'k.');
                errorbar(px, py, perr, perr, 'ko');
                hold on;
                plot(x_interp, y_interp, 'k-');
                hold on;
                plot(null_x, null_y, 'k--');
                hold off;
                yl = YLim;
                YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
                XLabel('Direction of Motion (deg)');
                YLabel('Response LFP');

                %now, print out some useful information in the upper subplot
                subplot(2, 1, 1);
                PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
                PrintDirectionData(p_value(i,j), base_rate, null_rate, amplitude, pref_dir, max_rate, width, DSI(i,j), stats1{i,j}, stats2{i,j}, DirDI(i,j), chi2(i,j), chiP(i,j));
            end
        end
    end
end

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', sprintf('%s: LFP at different bands',FILE));
subplot(2, 1, 2);
colormap(gray);
x = repmat(DirDI, [1 1 3]); %to use image in gray scale in need to replicate the matrix
x(isnan(x)) = 0; 
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('Direction Discrimination Index');

subplot(2, 1, 1);
colormap(gray);
x = repmat(DMI, [1 1 3]); %to use image in gray scale in need to replicate the matrix
x(isnan(x)) = 0; 
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('Direction Modulation Index');

toc
% keyboard


output = 0;
output2 = 1;
%output tuning curve metrics
if (output == 1)
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i) 'dir'];
    
    fileid = [PATHOUT FILEOUT];
    fwriteid = eval(['fopen(fileid, ''w'')']);
    %fprintf(fwriteid, '%%Base Rate (q1)	Amplitude (q2)	Pref dir (q3)	q4	Width(FWHM)	Max Resp	Spont Resp	DSI	Curve ANOVA	Mapped Pref Dir	Pref Speed	Pref H Disp	RF X-Ctr	RF Y-Ctr	Diam\n');
    fprintf(fwriteid, '%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f\n', pars(1), pars(2), pars(3), pars(4), width, max_rate, null_rate, DSI, p_value, data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED), data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER));

	fclose(fwriteid);

   %---------------------------------------------------------------------------------------
   %also write out data in form suitable for plotting tuning curve with Origin.
    FILEOUT2 = [FILE(1:i) 'direc_curv_fit'];
    fileid = [PATHOUT FILEOUT2];
    proffid = fopen(fileid, 'w');
    fprintf(proffid,'DirIn\tFit\tDirec\tAvgResp\tStdErr\tDir2\tSpon\n');
    for kk=1:length(x_interp)
        fprintf(proffid,'%6.2f\t%6.2f\t', x_interp(kk), y_interp(kk));
        if (kk <= length(px))
            fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t', px(kk), py(kk), perr(kk));
        else
            fprintf(proffid,'\t\t\t');
        end
        if (kk <= 2)
            fprintf(proffid,'%6.2f\t%6.2f\n',null_x(kk),null_y(kk));
        else
            fprintf(proffid,'\t\n');
        end
    end
    fclose(proffid);
    
    %---------------------------------------------------------------------------------------
    %ALso, write out summary data to a cumulative summary file
    [pdir360, pdir180, pdir90] = AngleWrap(pref_dir);
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.4f\t %6.3f\t %8.5f\t %10.8f\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        null_rate, DMI, base_rate, amplitude, pdir360, pdir180, pdir90, width, DSI, p_value, stats1(1), stats1(3), stats2(1), stats2(3), DirDI, var_term, chi2, chiP);
    outfile = [BASE_PATH 'ProtocolSpecific\DirectionTuning\DirectionTuningSummary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Spont\t DirMI\t BRate\t Ampl\t PD360\t PD180\t PD90\t FWHM\t DSI\t AnovaP\t\t Rmeans\t Pmeans\t\t Rraw\t Praw\t\t DirDI\t VarTrm\t Chi2\t\t ChiP\t\t');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %---------------------------------------------------------------------------------------
    
    %output a cumulative file of the Gaussian fit parameters
    outfile = [BASE_PATH 'ProtocolSpecific\DirectionTuning\DirectionParams.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fsummid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fsummid, 'FILE\t\t q(1)\t q(2)\t q(3)\t q(4)\t spont\t');
        fprintf(fsummid, '\r\n');
    end
    fprintf(fsummid, '%s\t %7.5f %7.5f %7.5f %7.5f %7.5f', FILE, pars(1), pars(2), pars(3), pars(4), null_rate);
    fprintf(fsummid, '\r\n');
    fclose(fsummid);
    
end

if output2 %save DirDI, DMI, anova p_val, chi2 p_val, and filenames to cumulative matfiles
    
    outfile = [BASE_PATH 'ProtocolSpecific\DirectionTuning\LFP_indices.mat'];
    if (exist(outfile, 'file') == 0)
        cum_DirDI{1} = DirDI;
        cum_DMI{1} = DMI;
        cum_p_anova{1} = p_value;
        cum_chiP{1} = chiP;
        files{1} = FILE;
        save(outfile, 'cum_DirDI', 'cum_DMI', 'cum_p_anova', 'cum_chiP', 'files');
    else
        load(outfile);
        cum_DirDI{length(cum_DirDI)+1} = DirDI;
        cum_DMI{length(cum_DMI)+1} = DMI;
        cum_p_anova{length(cum_p_anova)+1} = p_value;
        cum_chiP{length(cum_chiP)+1} = chiP;
        files{length(files)+1} = FILE;
        save(outfile, 'cum_DirDI', 'cum_DMI', 'cum_p_anova', 'cum_chiP', 'files');
    end

end
        
%print(2); % Uncomment for autoprinting.  JWN 081605
%close(2); % Uncomment for autoprinting.

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% After running all the individual trials and saving the data to
% LFP_indices.mat, run the following code to look at population data.
load LFP_indices.mat
DMI = zeros(19);
DirDI = zeros(19);
chiP = zeros(19);
p_anova = zeros(19);
n = length(cum_DirDI);
for i = 1:34
    DMI = DMI + cum_DMI{i}./34;
    DirDI = DirDI + cum_DirDI{i}./34;
    chiP = chiP + (cum_chiP{i}>0.05)./34;
    p_anova = p_anova + (cum_p_anova{i}<.05)./34;
end

figure
colormap(gray)
x = repmat(DMI,[1 1 3]);
x(isnan(x)) = 0;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('Direction Modulation Index');

figure
colormap(gray)
x = repmat(DirDI,[1 1 3]);
x(isnan(x)) = 0;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('Direction Discrimination Index');

figure
colormap(gray)
x = repmat(chiP,[1 1 3]);
x(isnan(x)) = 0;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('Fraction Sites with ChiP > 0.05');

figure
colormap(gray)
x = repmat(p_anova,[1 1 3]);
x(isnan(x)) = 0;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/5:max(yl)],'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('Fraction Sites with Significant Tuning (ANOVA)');



%%%%%%%%%%%%%%%%%%%%%%
% To compare the single unit tuning to the tuning in the various lfp bands,
% run the following code
su_metric = SU_p_anova;
lfp_matrix = cum_p_anova;
clear lfp_metric bandcorr;
for i = 1:size(lfp_matrix{1},1)
    for j = 1:size(lfp_matrix{1},2)
        for k = 1:length(files)
            lfp_metric{i,j}(k) = lfp_matrix{k}(i,j);
        end
        bandcorr(i,j) = corr(lfp_metric{i,j}(:),su_metric');
    end
end
figure
colormap(gray)
x = repmat(bandcorr./2 + 0.5,[1 1 3]);
x(isnan(x)) = 0;
image([20:10:200], [10:10:190], x);
cb = colorbar('Location','EastOutside');
yl = get(cb,'YLim');
set(cb,'YTick',[min(yl):range(yl)/4:max(yl)],'YTickLabel',{'-1' '-0.5' '0' '0.5' '1.0'});
xlabel('Upper Bound (Hz)'); ylabel('Lower Bound (Hz)');
title('LFP-SU Correlation for p_anova');

    