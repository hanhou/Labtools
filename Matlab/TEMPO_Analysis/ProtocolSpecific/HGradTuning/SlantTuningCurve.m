function SlantTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

%get the x_ctr and y_ctr to calculate eccentricity
x_ctr = data.one_time_params(RF_XCTR);
y_ctr = data.one_time_params(RF_YCTR);

eccentricity = sqrt((x_ctr^2) + (y_ctr^2));

%get the column of values of horiz. disparity magnitude in the dots_params matrix
mag_disp = data.dots_params(DOTS_HGRAD_MAG,BegTrial:EndTrial,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (mag_disp == data.one_time_params(NULL_VALUE)) );

unique_mag_disp = munique(mag_disp(~null_trials)');	
num_mag_disp = length(unique_mag_disp);

speed_dots = data.dots_params(DOTS_SPEED, BegTrial:EndTrial,PATCH1);
unique_speed = munique(speed_dots(~null_trials)');

coh_dots = data.dots_params(DOTS_COHER, BegTrial:EndTrial,PATCH1);
unique_coh = munique(coh_dots(~null_trials)');

%get the column of values of horiz. disparity angle of orientation in the dots_params matrix
disp_ang = data.dots_params(DOTS_HGRAD_ANGLE,BegTrial:EndTrial,PATCH1);
unique_disp_ang = munique(disp_ang(~null_trials)');

if ((length(unique_mag_disp)==1) & (length(unique_disp_ang) == 1))
    HDispTuningCurveSize(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    break
end


%get the column of mean disparity values
mean_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (mean_disp == LEYE_CONTROL) | (mean_disp == REYE_CONTROL) | (mean_disp == UNCORR_CONTROL) );

%display monoc control switch
no_display_monoc = 0;

%display monoc or not?
if no_display_monoc == 1
    unique_mean_disp = munique(mean_disp(~null_trials & ~control_trials)');
else
    unique_mean_disp = munique(mean_disp(~null_trials)');
end

num_mean_disp = length(unique_mean_disp);

%get the column of different aperture sizes
ap_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH1);
%do all sizes
all_sizes = 0
unique_ap_size = munique(ap_size(~null_trials)');
if all_sizes ~= 1
    unique_ap_size = unique_ap_size(length(unique_ap_size));
    num_ap_size = length(unique_ap_size);
else
    num_ap_size = length(unique_ap_size);
end

%get the average horizontal eye positions to calculate vergence
Leyex_positions = data.eye_positions(1, :);
Reyex_positions = data.eye_positions(3, :);

vergence = Leyex_positions - Reyex_positions;

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(mag_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stringarray = [];

%now, print out some useful information in the upper subplot
gen_data_fig = figure;
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);



%prepare the main graphing window where the graphs will be drawn
%graph = figure;
%set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [500 50 500 773], 'Name', 'Mag Disp Tuning Curve');
%ver_graph = figure;
%set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 50 500 773], 'Name', 'Tilt vs. Horizontal Vergence');

data_string = '';   
angle_out = '';
TDI_mdisp_out = '';
IndTDI_vs_ALLTDI = '';
%resp_data = []; verg_data=[];
p_val = zeros(length(unique_ap_size), length(unique_mean_disp));
pref_tilt = zeros(length(unique_ap_size), length(unique_mean_disp));

spike_rates = data.spike_rates(SpikeChan, :);

%calculate TDI here
for i=1:num_mag_disp
    response_select = logical(mag_disp == unique_mag_disp(i));
    plot_y = spike_rates(response_select & ~null_trials & ~control_trials & select_trials);
    plot_x = disp_ang(response_select & ~null_trials & ~control_trials & select_trials);
    
    [TDI(i), var_term] = Compute_DDI(plot_x, plot_y);
end



for i=1:length(unique_disp_ang)
    response_select = logical(disp_ang == unique_disp_ang(i));
    AvgTiltResp(i) = mean(spike_rates(response_select & ~null_trials & ~control_trials & select_trials));
end

[PrefTiltVal PrefTiltInd] = max(AvgTiltResp);

PrefTilt = unique_disp_ang(PrefTiltInd);

response_select = logical(disp_ang == PrefTilt);
plot_x = mag_disp(response_select & ~null_trials & ~control_trials & select_trials);
plot_y = spike_rates(response_select & ~null_trials & ~control_trials & select_trials);

%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
figure(gen_data_fig);
subplot(2, 1, 2);
[px, py, perr, spk_max, spk_min] = PlotTuningCurve(plot_x', plot_y', symbols{1}, lines{1}, 1, 1);

[PrefSlantVal PrefSlantInd] = max(py);
PrefSlant = px(PrefSlantInd);

norm_py = py/PrefSlantVal;

[SDI, var_term] = Compute_DDI(plot_x, plot_y);

null_x = [min(px) max(px)];
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
null_y = [null_rate null_rate];
hold on
plot(null_x, null_y, 'k--');

cell_x = {plot_x'};

[p,T,STATS,TERMS] = anovan(plot_y', cell_x, 'full', 3, {'MagDisp'}, 'off');
%MS_error = T{4, 6};
%MS_treatment = T{2, 6};
%F_index = MS_treatment/ (MS_error + MS_treatment);

XLabel('Mag Disp (deg/deg)');
YLabel('Response (spikes/sec)');

height = axis;
yheight = height(4);
string = sprintf('Ap Size = %4.1f', unique_ap_size(1));
text(height(1), 0.9*yheight, string, 'FontSize', 8);
string = sprintf('Pref Tilt = %3.2f', PrefTilt);
text(height(1), 0.8*yheight, string, 'FontSize', 8);  
string = sprintf('Pref Slant = %2.3f', PrefSlant);
text(height(1), 0.7*yheight, string, 'FontSize', 8);  
string = sprintf('p val = %1.4f', p);
text(height(1), 0.6*yheight, string, 'FontSize', 8);  
string = sprintf('SDI = %2.3f', SDI);
text(height(1), 0.5*yheight, string, 'FontSize', 8);  
slant_string = '';

%this prints out the data curves for plotting in origin
printcurves = 0;
if printcurves == 1
    %print out each individual tuning curve for origin
    pathsize = size(PATH,2) - 1;
    while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
        pathsize = pathsize - 1;
    end   
    PATHOUT = 'Z:\Users\Jerry\GradAnalysis\data_curves\';
    filesize = size(FILE,2) - 1;
    while FILE(filesize) ~='.'
        filesize = filesize - 1;
    end
    FILEOUT = [FILE(1:filesize) 'slant_curve'];
    fileid = [PATHOUT FILEOUT];
    printflag = 0;
    if (exist(fileid, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    proffid = fopen(fileid, 'a');
    if (printflag)
        fprintf(proffid,'Slant\tMDisp\tAvgResp\tStdErr\tSponX\tSponY\n');
        printflag = 0;
    end
    
    for go=1:length(px)
        if (go<=2)
            fprintf(proffid,'%6.2f\t%2.4f\t%6.2f\t%6.3f\t%6.2f\t%6.2f\n', px(go), unique_mean_disp(1), py(go), perr(go), null_x(go),null_y(go));
        else
            fprintf(proffid,'%6.2f\t%2.4f\t%6.2f\t%6.3f\t\t\n', px(go), unique_mean_disp(1), py(go), perr(go));
        end
    end
    fclose(proffid);
end %end printcurves;

line = sprintf('\t%3.1d\t%3.2f\t%2.2f\t%2.2f\t%3.2f\t%1.3f\t%2.4f\t%1.4f\t%1.3f\n', unique_ap_size(1), data.neuron_params(PREFERRED_DIRECTION, 1), x_ctr, y_ctr, PrefTilt, PrefSlant, unique_mean_disp(1), p, SDI);
slant_string = strcat(slant_string, line);


%this prints out important values for analysis
printme = 0;
if (printme == 1)
    PATHOUT = 'Z:\Users\jerry\gradanalysis\';
    
    line = sprintf('%s\t%d\t', FILE, SpikeChan);
    slant_string = strcat(line, slant_string);
    
    %print slant metrics
    outfile = [PATHOUT 'Slant_metrics_05.21.03.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid,'File\tChannel\tApSize1\tPrefDir1\txctr1\tyctr1\tPrefTilt\tPrefSlant\tMdisp\tp\tSDI\n')
        printflag = 0;
    end
    
    fprintf(fid, '%s', [slant_string]);
    fprintf(fid, '\r\n');
    fclose(fid);
end

%this prints out the TDI vs peak tilt response at each slant
printTDI = 0;
if printTDI == 1
    PATHOUT = 'Z:\Users\jerry\gradanalysis\';
    
    line = sprintf('%s\t%d\t', FILE, SpikeChan);
    %print slant metrics
    outfile = [PATHOUT 'Slant_TDI_metrics_05.22.03.dat'];
    printflag = 0;
    if (exist(fileid, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    %if (printflag)
        %fprintf(fid,'File\tTDI\tPeakResp\tMagDisp\tMdisp\n')
        %printflag = 0;
    %end
    
    for i=1:num_mag_disp
        fprintf(fid,'%1.4f\t%3.2f\t%1.4f\t%1.4f\n', TDI(i), norm_py(i), unique_mag_disp(i), unique_mean_disp(1));
    end
    fclose(fid);
end %end printTDI;