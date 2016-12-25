%-----------------------------------------------------------------------------------------------------------------------
%-- Surround Mapping curve will plot a polar plot for each set of surrounds -JDN 3/28/00
%-----------------------------------------------------------------------------------------------------------------------
function SurroundMappingCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

%get the column of disparity values for the surround patches
disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (disp == data.one_time_params(NULL_VALUE)) );
patch_off_condition = logical((disp == data.one_time_params(PATCH_OFF)));

unique_disp = munique(disp(~null_trials)');

%get the column of values of offset angle of the surround patches
ang = data.dots_params(DOTS_AP_OFF_ANG,BegTrial:EndTrial,PATCH1);
unique_ang = munique(ang(~null_trials)');

%get the column of different aperture sizes
ap_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH1);
unique_ap_size = munique(ap_size(~null_trials)');
ap1_size = unique_ap_size(2:length(unique_ap_size));

ap4_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH4);
ap4_size = munique(ap4_size(~null_trials)');


%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (disp == LEYE_CONTROL) | (disp == REYE_CONTROL) | (disp == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

stringarray = [];


control_px = unique_ang(2:length(unique_ang));

print_circle = 0;

if print_circle == 1
    figure
    set(gcf, 'Name', 'Surround Response');
    for i = 1:length(unique_disp)
        
        disp_select = logical(disp == unique_disp(i));
        
        plot_x = ang(disp_select & ~null_trials & ~control_trials & select_trials);
        plot_y = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);
        
        %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
        %the following code creates a circle with radius equal to the center only condition
        if px == -7777
            px = 0:.25:359.75;
            for j = 1:length(px)-1
                py = [py;py(1)];
            end
            px = [px'; px(1)];
            py = [py; py(1)];
        else
            px = [px; px(1)];
            py = [py; py(1)];
        end
        
        
        %figure
        hold on
        px = px * pi/180;
        polar(px, py, lines{i});
        leg_disp{i} = sprintf('%1.3f', unique_disp(i));

    end
    
    axhandles = get(gca, 'Children');
    numlines = length(unique_disp);
    leghandles = [];
    %legindex = 1;
    %for legendloop = numlines:-1:1
    %   leghandles(legindex) = axhandles((legendloop*3));
    %   legindex = legindex + 1;
    %end
    legend(leg_disp);
    
    %now, get the firing rate for NULL condition trials and add null circle to plot
    null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
    null_circ = [0 60 120 180 240 300 0]; 
    null_array = [null_rate null_rate null_rate null_rate null_rate null_rate null_rate];
    null_circ = null_circ * (pi/180);
    hold on;
    polar(null_circ, null_array, 'k*-');
    circ = findobj(gca, 'Type', 'line', 'Color', 'black', 'Marker', '*');
    set(circ, 'LineWidth', 2);
    plot(0, 0, 'k*');
    axis equal
    
    grid on	
end

plot_sur_tuning = 0;
if plot_sur_tuning == 1
    surr_tuning_fig=figure;
    
    %crf_axis = axes('position', [.35 .32 .3 .3]);
    surr_axis(1) = axes('position', [.68 .32 .3 .3]);
    surr_axis(2) = axes('position', [.52 .64 .3 .3]);
    surr_axis(3) = axes('position', [.18 .64 .3 .3]);
    surr_axis(4) = axes('position', [.02 .32 .3 .3]);
    surr_axis(5) = axes('position', [.18 0 .3 .3]);
    surr_axis(6) = axes('position', [.52 0 .3 .3]);
end

%now plot a horizontal disparity tuning curve for each point
sur_fig = figure
set(sur_fig,'PaperPosition', [.2 .2 8 10.7], 'Position', [450 50 560 673], 'Name', 'H. Disp Tuning Curves for each Surround Patch');
subplot(2, 1, 2)
hold on;

list_disp = [];
total_angles = [];

unique_ang = unique_ang(2:length(unique_ang));

for i=1:length(unique_ang)
    disptune{i} = sprintf('%s%d%s', 'surrtuning', i, '.mat');
    plot_val{i} = sprintf('%s%d%s', 'surr_plot_values', i, '.mat');
end

%create variable to hold all response values for each surround location
%we will use this to calculate mean response difference from center only condition
%JDN 02/10/03
disp_response = zeros(length(unique_ang), length(unique_disp)-1);

%get ctr only response
disp_select = logical(disp == unique_disp(1));

plot_x = ang(disp_select & ~null_trials & ~control_trials & select_trials);
plot_y = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);

%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
[px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);

ctr_only = py;

do_anovan = 1;
for i = 1:length(unique_ang)
    ang_select = logical(ang == unique_ang(i));
    
    plot_x = disp(ang_select & ~null_trials & ~patch_off_condition & ~control_trials & select_trials);
    plot_y = spike_rates(ang_select & ~null_trials & ~patch_off_condition & ~control_trials & select_trials);
    %save(plot_val{i}, 'plot_x', 'plot_y');
    DDI(i) = Compute_DDI(plot_x, plot_y);
    
    num_reps = length(find(plot_x==unique_disp(2)));
    
    %arrange into an array for anova
    spike = [];
    spike(:,1) = plot_x';
    spike(:,2) = plot_y';
    sortedspikes = sortrows(spike, [1]);
    
    if do_anovan == 1
        angle_array = zeros(length(sortedspikes), 1);
        angle_array = angle_array + i;
        total_angles(length(total_angles)+1:length(total_angles)+length(sortedspikes),:) = angle_array;
        list_disp(length(list_disp)+1:length(list_disp)+length(sortedspikes),:) = sortedspikes;
    else
        for j = 1:length(unique_disp)-1
            spikes_arranged_by_disparity(:,j) = sortedspikes((num_reps*(j-1)+1):j*num_reps,2);
        end
        
        spikes_arranged_by_disparity_by_angle((num_reps*(i-1)+1):i*num_reps, :) = spikes_arranged_by_disparity;
    end
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    figure(sur_fig)
    subplot(2, 1, 2)
    [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
    leg{i} = sprintf('%1.3f', unique_ang(i));
    hold on;
    
    disp_response(i, :) = py';
    
    if plot_sur_tuning ==1
        figure(surr_tuning_fig);
        axes(surr_axis(i));
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
        hold on;
        null_x = [min(px) max(px)];
        null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
        null_y = [null_rate null_rate];
        
        printcurves = 0;
        if printcurves == 1
            %print out each individual tuning curve for origin
            pathsize = size(PATH,2) - 1;
            while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
                pathsize = pathsize - 1;
            end   
            PATHOUT = [PATH(1:pathsize) 'Analysis\Grad\'];
            filesize = size(FILE,2) - 1;
            while FILE(filesize) ~='.'
                filesize = filesize - 1;
            end
            FILEOUT = [FILE(1:filesize) 'surmap_curves'];
            fileid = [PATHOUT FILEOUT];
            printflag = 0;
            if (exist(fileid, 'file') == 0)    %file does not yet exist
                printflag = 1;
            end
            proffid = fopen(fileid, 'a');
            if (printflag)
                fprintf(proffid,'HDisp\tAvgResp\tStdErr\tSpon\n');
                printflag = 0;
            end
            
            for go=1:length(px)
                if (go<=2)
                    fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t%6.2f\t%6.2f\n', px(go), py(go), perr(go), null_x(go),null_y(go));
                else
                    fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t\t\n', px(go), py(go), perr(go));
                end
            end
            fclose(proffid);
        end %end printcurves;
    end %end print_sur_tuning
end %end for each surr position

%calculate response difference between avg response at each location
mean_resp_ind = mean(disp_response, 2);
mean_resp_all = mean(mean_resp_ind);

per_sur = 0;
for i=1:length(mean_resp_ind)
    per_sur = per_sur + (1-mean_resp_ind(i)/ctr_only);
end

avg_per_sur = per_sur/length(mean_resp_ind);


ang_rad = unique_ang * pi / 180;
[x_cart, y_cart] = pol2cart(ang_rad, mean_resp_ind);
x_all = mean(x_cart);
y_all = mean(y_cart);

[t, r] = cart2pol(x_all, y_all);
r = r/mean_resp_all
[x_all, y_all] = pol2cart(t, r);

print_vector = 1;
if print_vector == 1
    figure
    compass(x_cart, y_cart, 'b');
    hold on
    compass(x_all, y_all, 'r');
    
    px_temp = 0:1:359;
    ctr_temp = ctr_only;
    ind_resp = mean_resp_all;
    for j = 1:length(px_temp)-1
        ctr_temp = [ctr_temp;ctr_temp(1)];
        ind_resp = [ind_resp; ind_resp(1)];
    end
    px_temp = [px_temp'; px_temp(1)];
    ctr_temp = [ctr_temp; ctr_temp(1)];
    ind_resp = [ind_resp; ind_resp(1)];
    px_temp = px_temp * pi/180;
    polar(px_temp, ctr_temp, 'r');
    polar(px_temp, ind_resp, 'k');
end

%avg level of surround change
surchange = mean_resp_all - ctr_temp(1);

%convert t to degrees
t = t * 180/pi;

printSur = 1;
if printSur == 1
    SUR_OUT = 'Z:\Users\Jerry\Gradanalysis\AvgPerSurrI_031003.dat';
    fSUR = fopen(SUR_OUT, 'a');
    fprintf(fSUR, '%s\t%3.4f\t%3.2f\t%1.4f\t%1.4f\n', FILE, avg_per_sur, t, r, surchange);
    fclose(fSUR);
end

%calculate metric of assymetry of DDI
meanDDI = mean(DDI);
[x_ddi, y_ddi] = pol2cart(ang_rad, DDI');
xddi_mean = mean(x_ddi);
yddi_mean = mean(y_ddi);
[tddi, rddi] = cart2pol(xddi_mean, yddi_mean);
rddi = rddi/meanDDI
[x_all, y_all] = pol2cart(tddi, rddi);

print_DDIgraph = 1;
if print_DDIgraph == 1
    figure
    compass(x_ddi, y_ddi, 'b');
    hold on
    compass(x_all, y_all, 'r');
end

printDDISur = 0;
if printDDISur == 1
    SUR_OUT = 'Z:\Users\Jerry\Gradanalysis\DDIVector.dat';
    fSUR = fopen(SUR_OUT, 'a');
    fprintf(fSUR, '%s\t%3.4f\t%3.2f\t%1.4f\n', FILE, meanDDI, tddi, rddi);
    fclose(fSUR);
end

printDDI = 0;
if printDDI == 1
    DDI_OUT = 'Z:\Users\Jerry\Gradanalysis\SurrDDI.dat';
    fDDI = fopen(DDI_OUT, 'a');
    fprintf(fDDI, '%s\t%3.0f\t%1.4f\t%3.0f\t%1.4f\t%3.0f\t%1.4f\t%3.0f\t%1.4f\t%3.0f\t%1.4f\t%3.0f\t%1.4f\n', FILE, unique_ang(1), DDI(1), unique_ang(2), DDI(2), unique_ang(3), DDI(3), unique_ang(4), DDI(4), unique_ang(5), DDI(5), unique_ang(6), DDI(6));
    fclose(fDDI);
end

figure(sur_fig);
subplot(2, 1, 2)
axhandles = get(gca, 'Children');
numlines = length(unique_ang);
leghandles = [];
legindex = 1;
for legendloop = numlines:-1:1
    leghandles(legindex) = axhandles((legendloop*3));
    legindex = legindex + 1;
end
legend(leghandles, leg);

if do_anovan == 1
    list_disp = [total_angles list_disp];
    p = anovan(list_disp(:, 3), {list_disp(:, 2) list_disp(:, 1)}, 'full', 3, {'H. Disp';'Angular Position'}, 'off');
else
    p = anova2(spikes_arranged_by_disparity_by_angle, num_reps, 'off')
end

%now, get the firing rate for NULL condition trials
null_x = [min(px) max(px)];
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');

%now, get the firing rate for RF only condition trials
rf_x = [min(px) max(px)];
rf_rate = mean(data.spike_rates(SpikeChan, patch_off_condition & select_trials));
rf_y = [rf_rate rf_rate];
hold on;
plot(rf_x, rf_y, 'b--');
grid on

height = axis;
yheight = height(4);
hold on
string = sprintf('P-Values = %2.4f %2.4f %2.4f', p);
text(height(1), 0.9*yheight, string, 'FontSize', 8);

string = sprintf('Inside Patch Size: %d ', ap4_size);
text(height(1), 0.8*yheight, string, 'FontSize', 8);

string = sprintf('Outside Patch Size: %d ', ap1_size);
text(height(1), 0.7*yheight, string, 'FontSize', 8);

hold off;

%now, print out some useful information in the upper subplot
figure(sur_fig);
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

return;

