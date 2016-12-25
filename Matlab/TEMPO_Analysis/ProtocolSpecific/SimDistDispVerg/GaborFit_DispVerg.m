%-----------------------------------------------------------------------------------------------------------------------
%-- GaborFit_SimVerg.m -- Fits disparity tuning curves with full Gabors, for each fixation depth.
%-- NOTE: responses are square-rooted in the error functions, to help homogenize variance (a la Cumming).  
%-- A constrained version of Levenberg-Marquardt optimization is used, as implemented using 'fmincon'.
%--	VR 08/04/04
%-- Modified from GaborFit_SimDist by TU 
%-- Modified from GaborFit_RelDisp by MLM
%-----------------------------------------------------------------------------------------------------------------------
function GaborFit_SimVerg(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

%determine which parameters are going to be shared. shared = 1, not shared = 0.
%note that these variables are redefined later, when multiple multigabor
%fits are performed.
Base_Rate = 1;
Amplitude = 1;
Gauss_Ctr = 1;
Gauss_SD = 1;
Sine_Freq = 1;
Sine_Phase = 1;

%determine whether disparity values will be adjusted to compensate error in
%actual vergence angles. 0=don't adjust; 1=adjust
adjust_disparities = 0;
targ_vergence = [-3.51816 0 1.75908]; %ideal angles according to Tempo, for 28.5 57 114 respectively

%determine whether shift ratios will be calculated for all combination. All combination = 1, relative to zero-disparity surround = 0
% all_combination == 2 -> compute DDI of tuning for all surround disparities, use largest DDI as reference for shift ratio JPM 3-21-02.
all_combination = 1;

% not implemented yet to select output
output = 0;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hdisp = data.dots_params(DOTS_HDISP,:,PATCH1);

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hdisp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%get the depth of fixation dot
depth_fix_real = data.dots_params(DEPTH_FIX_REAL,:,PATCH2);
unique_depth_fix_real = munique(depth_fix_real');

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hdisp == data.one_time_params(NULL_VALUE)));

unique_hdisp = munique(hdisp(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hdisp == data.one_time_params(NULL_VALUE)));

% Calculate spontaneous rate
null_rate = mean(spike_rates(null_trials & select_trials));

%get the average eye positions to calculate vergence
if (data.eye_calib_done == 1)
    Leyex_positions = data.eye_positions_calibrated(1, :);
    Leyey_positions = data.eye_positions_calibrated(2, :);
    Reyex_positions = data.eye_positions_calibrated(3, :);
    Reyey_positions = data.eye_positions_calibrated(4, :);
     
    vergence_h = Reyex_positions(trials) - Leyex_positions(trials);
    vergence_v = Reyey_positions(trials) - Leyey_positions(trials);
else     
    Leyex_positions = data.eye_positions(1, :);
    Leyey_positions = data.eye_positions(2, :);
    Reyex_positions = data.eye_positions(3, :);
    Reyey_positions = data.eye_positions(4, :);
    
    vergence_h = Reyex_positions(trials) - Leyex_positions(trials);
    vergence_v = Reyey_positions(trials) - Leyey_positions(trials);
end
    
%check if there are least 2 fixation distances
if (length(unique_depth_fix_real)<2) 
    disp('This analysis only for runs with at least two fixation distances interleaved');
    return;
end

%calculate vergence errors for each distance
% First fit separate gabors for each fixation distance.

DDI_index = zeros(1,length(unique_depth_fix_real));  % store DDI values for each fixation distance

for i=1:length(unique_depth_fix_real)	%for each different fixation distance, plot a separate disparity tuning curve.
    
    depth_fix_real_select = logical( (depth_fix_real == unique_depth_fix_real(i)) );
    verg_x_error(i) = mean(vergence_h(depth_fix_real_select & ~null_trials & select_trials)) + targ_vergence(i);
    if (adjust_disparities)
        plot_x = hdisp(depth_fix_real_select & ~null_trials & select_trials) + verg_x_error(i);
    else
        plot_x = hdisp(depth_fix_real_select & ~null_trials & select_trials);
    end
    plot_y = spike_rates(depth_fix_real_select & ~null_trials & select_trials);  

    %Calculate a mean firing rate for this condition.
    avg_resp(i)=mean(plot_y);
    
    [ DDI_index(i), var_term(i) ] = Compute_DDI(plot_x, plot_y);
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    p_value(i) = spk_anova(plot_y, plot_x, px);
    
    means{i} = [px py];
    raw{i} = [plot_x' plot_y'];
    sem{i} = [perr];
    verg{i} = vergence_h(depth_fix_real_select & ~null_trials & select_trials); 

%     fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
%     fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
%     
%     %fit the GABOR function
%     [pars{i},freq(i)] = gaborfit(means{i},raw{i},fixed_param_flags,fixed_param_values);
% 
%     subplot(5,2,2*i);
%     
%     %plot the raw data
%     PlotRawData(plot_x', plot_y', symbols{i}, null_rate);
%   
%     x_interp = (px(1)): .01 : (px(length(px)));
%     y_gabor = gaborfunc(x_interp, pars{i});
%     y_gabor(y_gabor < 0) = 0;
%     y_gauss =  pars{i}(1) + pars{i}(2)*exp(-0.5*((x_interp - pars{i}(3))/ pars{i}(4)).^2);
%     y_sine =  pars{i}(1) + 0.5*pars{i}(2)*cos(2*pi*pars{i}(5)*(x_interp - pars{i}(3))+pars{i}(6) );
%     
%     %Compute peak and trough disparity
%     [peak_rate(i) peak_index] = max(y_gabor);
%     peak_disp(i) = x_interp(peak_index);
%     [trough_rate(i) trough_index] = min(y_gabor);
%     trough_disp(i) = x_interp(trough_index);
%     
%     % Now plot fitted curves; square the fitted values to counteract sqrt() above
%     hold on;
%     plot(x_interp, y_gabor, 'k-');
%     plot(x_interp, y_gauss, 'k--');
%     plot(x_interp, y_sine, 'k:');
%     Gabor_SSE(i) = gaborerr(pars{i});
%     buff = sprintf('Full Gabor fit (sse = %8.5f)', gaborerr(pars{i}) );
%     title(buff);
%     hold off;   
%     
%     y_fit = gaborfunc(px, pars{i});
%     y_fit(y_fit < 0) = 0;
%     %add a column of ones to yfit to make regress happy
%     y_fit = [ones(length(y_fit),1) y_fit];
%     [b, bint, r, rint, stats1{1}{i}] = regress(py, y_fit);
%     
%     y_fit_raw = gaborfunc(plot_x', pars{i});
%     y_fit_raw(y_fit_raw < 0) = 0;
%     y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
%     [b, bint, r, rint, stats2{1}{i}] = regress(plot_y', y_fit_raw);
end


% for i=1:length(unique_depth_fix_real)	%for each different fixation distance, plot some information.
%     
%     %print fit parameters to the screen
%     subplot(5,2,2*i-1);
%     axis([0 100 0 100]);    axis('off');
%     xpos = -10; ypos = 100;
%     font_size = 8;  bump_size = 10;
% 
%     if(i == 1)
%         temp = strcat(PATH, FILE);
%         temp(temp == '\') = '/';
%         % this prevents a stupid error from appearing on the screen
%         line = sprintf('File: %s', temp);
%         text(xpos-15,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     end
%     line = sprintf('Fixation Distance: %6.2f cm', unique_depth_fix_real(i));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Base Rate: %8.4f', pars{i}(1));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Amplitude: %8.4f', pars{i}(2));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Gauss Ctr: %8.4f', pars{i}(3));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Gauss  SD: %8.4f', pars{i}(4));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Sine Freq: %8.4f', pars{i}(5));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Sin Phase: %8.4f', pars{i}(6));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('FT Freq: %8.4f', freq(i));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{1}{i}(1), stats1{1}{i}(3));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
%     line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{1}{i}(1), stats2{1}{i}(3));
%     text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% 
%     %Omitted the section about shift calculations MLM
%     
% end

%To perform multi gabor fits with different sets of constrained parameters,
%loop over the following code and store parameters in different subscripts.
%1 indicates parameter is shared; 0 indicates parameter is free to vary.
% Fit 1 = share all parameters; 
% Fit 2 = allow Ctr (i.e., Gaussian Ctr) to vary; 
% Fit 3 = allow Amplitude to vary
% Fit 4 = allow Ctr and Amplitude to vary
num_fits = 5;
Base_Rate = [1 1 1 1 0]; Amplitude = [1 1 0 0 0]; Gauss_Ctr = [1 0 1 0 0]; Gauss_SD = [1 1 1 1 0]; Sine_Freq = [1 1 1 1 0]; Sine_Phase = [1 1 1 1 0];
Fit_labels = [{'ShAll'} {'Ctr'} {'Amp'} {'CtrAmp'} {'Uncons'}];

for h=1:num_fits    
    %Next, fit multiple gabors with shared parameters.
    figure;
    set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'HDisp Tuning Curve w/ Best Fit');
    
    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    shared_param_flags{h} = [Base_Rate(h) Amplitude(h) Gauss_Ctr(h) Gauss_SD(h) Sine_Freq(h) Sine_Phase(h)]; %shared parameters are defined above
    
    %fit gabor functions to multiple curves with shared parameters
    [pars_fit{h},freq_fit{h}] = multi_gaborfit(means,raw,fixed_param_flags,fixed_param_values,shared_param_flags{h});
    fit_err(h) = multi_gaborerr(pars_fit{h});
    
    %re-structure the parameters
    k = 1;
    params{h}{1} = pars_fit{h}(1:6);
    for i=2:length(unique_depth_fix_real)
        params{h}{i} = pars_fit{h}(1:6);
        for j=1:6
            if(shared_param_flags{h}(j) == 0)
                params{h}{i}(j) = pars_fit{h}(6 + k);
                k = k + 1;
            end
        end
    end
    
    for i=1:length(unique_depth_fix_real)	%for each different fixation distance, plot a separate disparity tuning curve.
        
        subplot(5,2,2*i);
        
        %plot the raw data
        PlotRawData(raw{i}(:,1), raw{i}(:,2), symbols{i}, null_rate);
        
        x_interp = (means{i}(1,1)): .01 : (means{i}(length(means{i}),1));
        y_gabor = gaborfunc(x_interp, params{h}{i});
        y_gabor(y_gabor < 0) = 0;
        y_gauss =  params{h}{i}(1) + params{h}{i}(2)*exp(-0.5*((x_interp - params{h}{i}(3))/ params{h}{i}(4)).^2);
        y_sine =  params{h}{i}(1) + 0.5*params{h}{i}(2)*cos(2*pi*params{h}{i}(5)*(x_interp - params{h}{i}(3))+params{h}{i}(6) );
        
        % Now plot fitted curves; square the fitted values to counteract sqrt() above
        hold on;
        plot(x_interp, y_gabor, 'k-');
        plot(x_interp, y_gauss, 'k--');
        plot(x_interp, y_sine, 'k:');
        buff = sprintf('%s Gabor fit (sse = %8.5f)', Fit_labels{h}, multi_gaborerr(pars_fit{h}) );
        title(buff);
        hold off;   
        
        y_fit = gaborfunc(means{i}(:,1), params{h}{i});
        y_fit(y_fit < 0) = 0;
        %add a column of ones to yfit to make regress happy
        y_fit = [ones(length(y_fit),1) y_fit];
        [b, bint, r, rint, stats1{h}{i}] = regress(means{i}(:,2), y_fit);
        
        y_fit_raw = gaborfunc(raw{i}(:,1), params{h}{i});
        y_fit_raw(y_fit_raw < 0) = 0;
        y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
        [b, bint, r, rint, stats2{h}{i}] = regress(raw{i}(:,2), y_fit_raw);
        
        %calculate chi-squared goodness of fit
        [chi2{h}(i) chiP{h}(i)] = chi2_test(raw{i}(:,1), raw{i}(:,2), 'gaborfunc', ...
            params{h}{i}, (sum(shared_param_flags{h})));
        
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
        line = sprintf('Fixation Dist: %6.2f cm', unique_depth_fix_real(i));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Base Rate: %8.4f', params{h}{i}(1));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Amplitude: %8.4f', params{h}{i}(2));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Gauss Ctr: %8.4f', params{h}{i}(3));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Gauss  SD: %8.4f', params{h}{i}(4));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Sine Freq: %8.4f', params{h}{i}(5));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Sin Phase: %8.4f', params{h}{i}(6));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('FT Freq: %8.4f', freq_fit{h}(i));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{h}{i}(1), stats1{h}{i}(3));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{h}{i}(1), stats2{h}{i}(3));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Goodness of Fit: Chi2: %6.3f, P=%8.6f', chi2{h}(i), chiP{h}(i));
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        
        %again omit shift ratio here
        
    end
    
end

% Sequential F-tests

% Recall Fit #s from above, with addition of totally unconstrained fits as the final fit.
% Fit 1 = share all parameters; 
% Fit 2 = allow Ctr (i.e., Gaussian Ctr) to vary; 
% Fit 3 = allow Amplitude to vary
% Fit 4 = allow Ctr and Amplitude to vary
% Fit 5 = all parameters free.

% Define contrasts as cells within a matrix each containing an ordered pair, where
% the ordered pair refers to the two Fits.  e.g., [ {[1,5]} ] represents a
% F-test of a fully constrained fit against a fully free fits

Ftests = [ {[1,5]} {[2,5]} {[3,5]} {[4,5]} {[2,4]} {[3,4]} {[1,2]} {[1,3]} ];
Ftest_labels = [ {'AllVsUncons'} {'CtrVsUncons'} {'AmpVsUncons'} {'CtrAmpVsUncons'} {'CtrVsCtrAmp'} ...
        {'AmpVsCtrAmp'} {'AllVsCtr'} {'AllVsAmp'} ]; 
num_Ftests = length(Ftests);

%n_free_fit_params(num_shared_fits+1) = ...                          % moves total # independent parameters  
%    length(shared_param_flags{1})*length(unique_depth_fix_real);    % to correspond to unconstrained fit
for i=1:num_fits
%    shared_err(i) = multi_gaborerr(pars_shared{i}); %now calculated
%    immediately after gaborfit
    n_free_fit_params(i) = ( length(shared_param_flags{1})*length(unique_depth_fix_real) ) - ...
        ( (length(unique_depth_fix_real)-1)*sum(shared_param_flags{i}) );
end
%fit_err(num_shared_fits+1) = sum(Gabor_SSE); % moves independent error to end of fit_err

Npts = length(spike_rates(~null_trials & select_trials));

%compute and print sequential F-test for each defined contrast.
for i=1:num_Ftests
    Fit1 = Ftests{i}(1); Fit2 = Ftests{i}(2);
    Fseq(i) = ( (fit_err(Fit1) - fit_err(Fit2) )/(n_free_fit_params(Fit2)-n_free_fit_params(Fit1)) ) / ...
        (fit_err(2)/(Npts - n_free_fit_params(Fit2)));
    Pseq(i) = 1 - fcdf(Fseq(i), (n_free_fit_params(Fit2)-n_free_fit_params(Fit1)), (Npts-n_free_fit_params(Fit2)));

    ypos = ypos - bump_size;    ypos = ypos - bump_size;
    line = sprintf('Sequential F-test - %s', Ftest_labels{i});
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('      F = %7.3f, P = %7.5f', Fseq(i), Pseq(i)); 
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('      N-indep = %2d, N-constr = %2d', n_free_fit_params(Fit2), n_free_fit_params(Fit1));
    text(xpos,ypos,line,'FontSize',font_size);		%ypos = ypos - bump_size;
end

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
FILEOUT = [FILE(1:i) 'fixdist_curv'];
fileid = [PATHOUT FILEOUT];
proffid = eval(['fopen(fileid, ''w'')']);
 
fprintf(proffid,'FitX ');

for i=1:(num_fits)
    for j=1:length(unique_depth_fix_real)
        fprintf(proffid,'%s%d ',Fit_labels{i},j);
    end
end

fprintf(proffid,'HDisp ');
for i=1:length(unique_depth_fix_real)
    fprintf(proffid,'AvgResp%d StdErr%d ', i, i);
end    
fprintf(proffid,'FixDist2 Spon\n');

% calculate fit values for each fit and fixation distance
for i=1:(num_fits) 
    for j=1:length(unique_depth_fix_real)
        y_gaborfit{i}{j} = gaborfunc(x_interp, params{i}{j});
        y_gaborfit{i}{j}(y_gaborfit{i}{j} < 0) = 0;
    end
end

for i=1:length(x_interp)
    fprintf(proffid,'%6.2f ', x_interp(i));
    for j=1:(num_fits)
        for k=1:length(unique_depth_fix_real)
            fprintf(proffid,'%6.2f ', y_gaborfit{j}{k}(i));
        end
    end    
    if (i<=length(sem{1}))    
        fprintf(proffid,'%6.2f ', means{1}(i,1));
        for j=1:length(unique_depth_fix_real)
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

%------------------------------------------------------------------------
%write out all relevant parameters to a cumulative text file, GCD 8/08/01
%adapted for this file, MLM 7/18/2002
%write out one line for each fixation distance for each neuron.
outfile = [BASE_PATH 'ProtocolSpecific\SimDistDispVerg\GaborFitSummary.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILENAME\t PrDir\t PrSpd\t PrHDsp\t RFXCtr\t RFYCtr\t RFDiam\t FixDst\t VergErrX\t MaxDsp\t MaxRsp\t MinDsp\t MinRsp\t AvgRsp\t Spont\t DDI\t VarTrm\t'); 
    for i=1:(num_fits)
        fprintf(fid, '%s:BasRt\t %s:Amp\t %s:Ctr\t %s:Siz\t %s:Freq\t %s:Phase\t %s:FFTFr\t %s:Rsq_Means\t %s:Rsq_Raw\t %s:Chi2\t %s:ChiP\t',Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i},Fit_labels{i});
    end
    for i=1:num_Ftests
        fprintf(fid, 'Fseq-%s\t Pseq-%s\t', Ftest_labels{i}, Ftest_labels{i});
    end
    %fprintf(fid, 'Peak\t Trough');
    fprintf(fid, '\r\n');
    printflag = 0;
end

for j = 1:length(unique_depth_fix_real)
    outbuff = [];
    outbuff{1} = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t ', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1) );
    outbuff{2} = sprintf('%6.1f\t %6.3f\t %6.3f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.3f\t ', ...
        unique_depth_fix_real(j), verg_x_error(j), pmax(j).x, pmax(j).y, pmin(j).x, pmin(j).y, avg_resp(j), null_rate, DDI_index(j), var_term(j));
    %outbuff{3} = sprintf('%7.2f\t %7.2f\t %6.3f\t %7.4f\t %6.3f\t %6.3f\t %6.3f\t ', ...
    %    pars_fit{j}(1), pars_fit{j}(2), pars_fit{j}(3), pars_fit{j}(4), pars_fit{j}(5), pars_fit{j}(6), freq_fit(j) );
    for i=1:num_fits
        outbuff{2+i} = sprintf('%7.2f\t %7.2f\t %6.3f\t %7.4f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %8.6f\t', ...
            params{i}{j}(1), params{i}{j}(2), params{i}{j}(3), params{i}{j}(4), params{i}{j}(5), params{i}{j}(6), freq_fit{i}(j), ...
            stats1{i}{j}(1), stats2{i}{j}(1), chi2{i}(j), chiP{i}(j) );
    end
    for i=1:num_Ftests
        outbuff{3+num_fits+i} = sprintf('%7.4f\t %7.4f\t', Fseq(i), Pseq(i));
    end
    %outbuff(5+num_shared_fits+num_Ftests) = sprintf('%7.4f\t %7.4f\t', peak_disp(j), trough_disp(j));
    fprintf(fid, '%s', [ outbuff{:} ]);
    fprintf(fid, '\r\n');
end
fclose(fid);
%------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------
%also write out data for summary.
outfile1 = [BASE_PATH 'ProtocolSpecific\SimDistDispVerg\Seq_Ftest_GCtr.dat'];

printflag = 0;
if (exist(outfile1, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile1, 'a');
%save some values to file
outstr1 = sprintf('File');
outstr2 = sprintf('%s', FILE);
for i = 1:num_Ftests
    outstr1 = sprintf('%s Fseq-%s Pseq-%s', outstr1, Ftest_labels{i}, Ftest_labels{i});
    outstr2 = sprintf('%s %8.4f %8.6f', outstr2, Fseq(i), Pseq(i));
end

if (printflag)
    fprintf(fid, '%s', outstr1);
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', outstr2);
fprintf(fid, '\r\n');
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

