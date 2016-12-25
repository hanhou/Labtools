%-----------------------------------------------------------------------------------------------------------------------
%-- GaborFit_SimDist.m -- Fits disparity tuning curves with full Gabors, for each individual surround distance.
%-- NOTE: responses are square-rooted in the error functions, to help homogenize variance (a la Cumming).  
%-- A constrained version of Levenberg-Marquardt optimization is used, as implemented using 'fmincon'.
%--	TU 07/21/01
%-- Modified from GaborFit_RelDisp by MLM
%-----------------------------------------------------------------------------------------------------------------------
function GaborFit_SimDist(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

%determine which parameters are going to be shared. shared = 1, not shared = 0.
Base_Rate = 1;
Amplitude = 1;
Gauss_Ctr = 1;
Gauss_SD = 1;
Sine_Freq = 1;
Sine_Phase = 1;

%determine whether shift ratios will be calculated for all combination. All combination = 1, relative to zero-disparity surround = 0
% all_combination == 2 -> compute DDI of tuning for all surround disparities, use largest DDI as reference for shift ratio JPM 3-21-02.
all_combination = 1;

% not implemented yet to select output
output = 0;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp_ctr = data.dots_params(DOTS_HDISP,:,PATCH1);
%get the simulated depth of the surround
dist_sim_sur = data.dots_params(DEPTH_DIST_SIM,:,PATCH2);

%get indices of any NULL conditions (for measuring spontaneous activity)
%changed & to | by MLM
null_trials = logical( (hor_disp_ctr == data.one_time_params(NULL_VALUE)) | (dist_sim_sur == data.one_time_params(NULL_VALUE)) );

% Ctr_Control_trials = no dots in ctr, surround varies
ctr_control_trials = logical(hor_disp_ctr == data.one_time_params(PATCH_OFF));
   
% Sur_Control_trials = no dots in surround, ctr disp varies
sur_control_trials = logical(dist_sim_sur == data.one_time_params(PATCH_OFF));
   
control_trials = (ctr_control_trials | sur_control_trials);

unique_dist_sim_sur = munique(dist_sim_sur(~null_trials & ~sur_control_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp_ctr);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rate
null_rate = mean(spike_rates(null_trials & select_trials));

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
    
%check if there are least surround simulated distances
if (length(unique_dist_sim_sur)<2) 
    disp('This analysis only for runs with at least two surround simulated distances interleaved');
    return;
end

% First fit separate gabors for each surround simulated distance.
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'HDisp Tuning Curve w/ Best Fit');

DDI_index = zeros(1,length(unique_dist_sim_sur));  % store DDI values for each surround simulated distance

for i=1:length(unique_dist_sim_sur)	%for each different surround simulated distance, plot a separate disparity tuning curve.
    
    dist_sim_sur_select = logical( (dist_sim_sur == unique_dist_sim_sur(i)) );
    
    plot_x = hor_disp_ctr(dist_sim_sur_select & ~null_trials & ~control_trials & select_trials);
    plot_y = spike_rates(dist_sim_sur_select & ~null_trials & ~control_trials & select_trials);  

    %Calculate a mean firing rate for this condition.
    avg_resp(i)=mean(plot_y);
    
    [ DDI_index(i), var_term(i) ] = Compute_DDI(plot_x, plot_y);
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    p_value(i) = spk_anova(plot_y, plot_x, px);
    
    means{i} = [px py];
    raw{i} = [plot_x' plot_y'];
    sem{i} = [perr];
    verg{i} = vergence_h(dist_sim_sur_select & ~null_trials & ~control_trials & select_trials); 

    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    
    %fit the GABOR function
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

outstr2 = [0];
for i=1:length(unique_dist_sim_sur)	%for each different surround disparity, plot some information.
    
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
    line = sprintf('Dist Sim: %6.2f cm', unique_dist_sim_sur(i));
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

    %Omitted the section about shift calculations MLM
    
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
for i=2:length(unique_dist_sim_sur)
    params{i} = pars_shared(1:6);
    for j=1:6
        if(shared_param_flags(j) == 0)
            params{i}(j) = pars_shared(6 + k);
            k = k + 1;
        end
    end
end

for i=1:length(unique_dist_sim_sur)	%for each different surround simulated distance, plot a separate disparity tuning curve.
    
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
    line = sprintf('Dist Sim: %6.2f cm', unique_dist_sim_sur(i));
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
 
    %again omit shift ratio here
    
end

%do Sequential F-test
independent_err = sum(Gabor_SSE);
shared_err = multi_gaborerr(pars_shared);
n_indep_params = length(shared_param_flags)*length(unique_dist_sim_sur);
n_shared_params = n_indep_params-(length(unique_dist_sim_sur)-1)*sum(shared_param_flags);
Npts = length(spike_rates(~null_trials & ~control_trials & select_trials));
Fseq = ( (shared_err - independent_err )/(n_indep_params-n_shared_params) ) / ( independent_err/(Npts - n_indep_params) );
Pseq = 1 - fcdf(Fseq, (n_indep_params-n_shared_params), (Npts-n_indep_params) );

ypos = ypos - bump_size;
line = sprintf('Sequential F-test: F = %7.3f, P = %7.5f', Fseq, Pseq);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('      N-indep = %2d, N-constr = %2d', n_indep_params, n_shared_params);
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
FILEOUT = [FILE(1:i) 'simdist_curv'];
fileid = [PATHOUT FILEOUT];
proffid = eval(['fopen(fileid, ''w'')']);
 
fprintf(proffid,'FitX ');
for i=1:length(unique_dist_sim_sur)
    fprintf(proffid,'GaborFit%d ', i);
end    
fprintf(proffid,'HDisp ');
for i=1:length(unique_dist_sim_sur)
    fprintf(proffid,'AvgResp%d StdErr%d ', i, i);
end    
fprintf(proffid,'SimDist2 Spon\n');

for i=1:length(x_interp)
    fprintf(proffid,'%6.2f ', x_interp(i));
    for j=1:length(unique_dist_sim_sur)
        y_gabor = gaborfunc(x_interp, pars{j});
        y_gabor(y_gabor < 0) = 0;
        fprintf(proffid,'%6.2f ', y_gabor(i));
    end    
    if (i<=length(sem{1}))    
        fprintf(proffid,'%6.2f ', means{1}(i,1));
        for j=1:length(unique_dist_sim_sur)
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
%write out one line for each surround simulated distance for each neuron.
outfile = [BASE_PATH 'ProtocolSpecific\SimDistDispOnly\GaborFitSummary.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILENAME\t PrDir\t PrSpd\t PrHDsp\t RFXCtr\t RFYCtr\t RFDiam\t SimDst\t MaxDsp\t MaxRsp\t MinDsp\t MinRsp\t AvgRsp\t Spont\t DDI  \t VarTrm\t BasRt\t Amp\t Ctr\t Siz\t Freq\t Phase\t FFTFr\t ShBasRt\t ShAmp\t ShCtr\t ShSiz\t ShFr\t ShPh\t ShFFTFr\t Fseq\t Pseq\t Peak\t Trough\t');
    fprintf(fid, '\r\n');
    printflag = 0;
end
for j = 1:length(unique_dist_sim_sur)
    outbuff = [];
    outbuff{1} = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t ', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1) );
    outbuff{2} = sprintf('%6.1f\t %6.3f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.3f\t ', ...
        unique_dist_sim_sur(j), pmax(j).x, pmax(j).y, pmin(j).x, pmin(j).y, avg_resp(j), null_rate, DDI_index(j), var_term(j));
    outbuff{3} = sprintf('%7.2f\t %7.2f\t %6.3f\t %7.4f\t %6.3f\t %6.3f\t %6.3f\t ', ...
        pars{j}(1), pars{j}(2), pars{j}(3), pars{j}(4), pars{j}(5), pars{j}(6), freq(j) );
    outbuff{4} = sprintf('%7.2f\t %7.2f\t %6.3f\t %7.4f\t %6.3f\t %6.3f\t %6.3f\t', ...
        params{j}(1), params{j}(2), params{j}(3), params{j}(4), params{j}(5), params{j}(6), freq_shared(j));
    outbuff{5} = sprintf('%7.4f\t %7.4f\t %7.4f\t %7.4f\t', Fseq, Pseq, peak_disp(j), trough_disp(j));
    fprintf(fid, '%s', [ outbuff{:} ]);
    fprintf(fid, '\r\n');
end
fclose(fid);
%------------------------------------------------------------------------

%----------------------------------------------------------------------------------------------------------------------------------------------------------
%also write out data for summary.
outfile1 = [BASE_PATH 'ProtocolSpecific\SimDistDispOnly\Seq_Ftest_GCtr.dat'];

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