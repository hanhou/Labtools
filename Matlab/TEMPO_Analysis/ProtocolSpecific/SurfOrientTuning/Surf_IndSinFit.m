%-----------------------------------------------------------------------------------------------------------------------
%-- HGradTuningCurve.m -- Plots a horizontal disparity gradient tuning curve.  These tuning curves will plot
%--   varying angles of gradient rotation vs. responses for different mean disparities on a single graph.  multiple
%--   graphs in a single column represent different gradient magnitudes for one aperture size.  Graphs in 
%--   different columns differ by aperture size.  All graphs include	monoc and uncorrelated control conditions.
%--	JDN 2/4/00
%-----------------------------------------------------------------------------------------------------------------------
function Surf_IndSinFit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

%get the x_ctr and y_ctr to calculate eccentricity
x_ctr = data.one_time_params(RF_XCTR);
y_ctr = data.one_time_params(RF_YCTR);

eccentricity = sqrt((x_ctr^2) + (y_ctr^2));

%--------------------------------------------------------------------------
%get all variables
%--------------------------------------------------------------------------
%get entire list of slants for this experiment
slant_list = data.dots_params(DOTS_SLANT,BegTrial:EndTrial,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (slant_list == data.one_time_params(NULL_VALUE)) );

unique_slant = munique(slant_list(~null_trials)');	
num_slant = length(unique_slant);

%get entire list of tilts for this experiment
tilt_list = data.dots_params(DOTS_TILT,BegTrial:EndTrial,PATCH1);
unique_tilt = munique(tilt_list(~null_trials)');

%get list of Stimulus Types
stim_list = data.dots_params(DOTS_STIM_TYPE, BegTrial:EndTrial, PATCH1);
unique_stim = munique(stim_list(~null_trials)');

%get motion coherency value
coh_dots = data.dots_params(DOTS_COHER, BegTrial:EndTrial,PATCH1);
unique_coh = munique(coh_dots(~null_trials)');

%get the column of mean depth values
mean_depth_list = data.dots_params(DEPTH_DIST_SIM,BegTrial:EndTrial,PATCH1);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (mean_depth_list == LEYE_CONTROL) | (mean_depth_list == REYE_CONTROL) | (mean_depth_list == UNCORR_CONTROL) );

%display monoc control switch
no_display_monoc = 0;

%display monoc or not?
if no_display_monoc == 1
    unique_mean_depth = munique(mean_depth_list(~null_trials & ~control_trials)');
else
    unique_mean_depth = munique(mean_depth_list(~null_trials)');
end

num_mean_depth = length(unique_mean_depth);

%get the column of different aperture sizes
ap_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH1);
unique_ap_size = munique(ap_size(~null_trials)');

%get the average horizontal eye positions to calculate vergence
Leyex_positions = data.eye_positions(1, :);
Reyex_positions = data.eye_positions(3, :);

vergence = Leyex_positions - Reyex_positions;

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(slant_list);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stringarray = [];

%now, print out some useful information in the upper subplot
gen_data_fig = figure;
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%prepare the main graphing window where the graphs will be drawn
graph = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [500 50 500 773], 'Name', 'Tilt Tuning Curve');
ver_graph = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 50 500 773], 'Name', 'Tilt vs. Horizontal Vergence');

data_string = '';   
angle_out = '';
TDI_mdisp_out = '';
IndTDI_vs_ALLTDI = '';
%resp_data = []; verg_data=[];
p_val = zeros(length(unique_stim), length(unique_mean_depth));
pref_tilt = zeros(length(unique_stim), length(unique_mean_depth));

%for each stim type, plot out the tilt tuning curve and vergence data
%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%store one time data values (Ap Size, X ctr, Y ctr, ecc, pref dir)
line = sprintf('\t%3.1f\t%2.5f\t%2.5f\t%2.5f\t%3.2f', unique_ap_size(1), x_ctr, y_ctr, eccentricity, data.neuron_params(PREFERRED_DIRECTION, 1));
data_string = strcat(data_string, line);

curve_out = cell(1000,1);
for stim_count = 1:length(unique_stim)
    TDIdata = [];
    TDIvergdata = [];
    ancova_var = [];
    list_angles = [];
    total_mdepth = [];
    verg_ancova_var = [];
    verg_list_angles = [];
    verg_total_mdepth = [];
    for slant_count = 1:length(unique_slant)
        m_disp_max = zeros(length(unique_mean_depth), 1);
        m_disp_min = zeros(length(unique_mean_depth), 1);
        start = zeros(length(unique_mean_depth), 1);
        stop = zeros(length(unique_mean_depth), 1);
        start_verg = zeros(length(unique_mean_depth), 1);
        stop_verg = zeros(length(unique_mean_depth), 1);
        for mdepth_count=1:length(unique_mean_depth)
            spike = [];
            verg = [];
            figure(graph);
            hold on;
            subplot(length(unique_stim), num_slant, stim_count);
                
            depth_select = logical((stim_list == unique_stim(stim_count)) &(slant_list == unique_slant(slant_count)) & (mean_depth_list == unique_mean_depth(mdepth_count)));
            plot_x = tilt_list(depth_select & ~null_trials & select_trials);
            plot_y = spike_rates(depth_select & ~null_trials & select_trials);
            ver = vergence(depth_select & ~null_trials & select_trials);
            
            %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
            [px, py, perr, spk_max, spk_min] = PlotTuningCurve(plot_x', plot_y', symbols{mdepth_count}, '', 1, 1);
            p_val(stim_count,mdepth_count) = calc_anovap(plot_x, plot_y);
            [value, index_max] = max(py);
            pref_tilt(stim_count,mdepth_count) = px(index_max);
            
            %save out each curve so that we can 
            %mean shift them to calculate an avg TDI value
            start(mdepth_count) = length(TDIdata)+1;
            stop(mdepth_count) = length(plot_x)+start(mdepth_count)-1;
            TDIdata(start(mdepth_count):stop(mdepth_count), 1) = plot_x';
            TDIdata(start(mdepth_count):stop(mdepth_count), 2) = plot_y';
            
            %also calculate the TDI for the individual curve
            [single_TDI(stim_count,mdepth_count), var_term] = Compute_DDI(plot_x, plot_y);
            
            %------Spike ANOVAN code-----------------------------------------
            %save out each data point to use in ANOVAN function
            spike(:,1) = plot_x';
            spike(:,2) = plot_y';
            sortedspikes = sortrows(spike, [1]);
            ancova_var(length(ancova_var)+1:length(ancova_var)+length(sortedspikes),:) = sortedspikes;
            for temp_tilt = 1:length(unique_tilt)
                tilt_ind = find(sortedspikes(:,1) == unique_tilt(temp_tilt));
                sortedspikes(tilt_ind(1):tilt_ind(length(tilt_ind)), 1) = temp_tilt;
            end
            %to do anovan
            mdepth_array = zeros(length(sortedspikes), 1);
            mdepth_array = mdepth_array + mdepth_count;
            total_mdepth(length(total_mdepth)+1:length(total_mdepth)+length(sortedspikes),:) = mdepth_array;
            list_angles(length(list_angles)+1:length(list_angles)+length(sortedspikes),:) = sortedspikes;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            %fit the data with a distorted sin wave
            means_px = (px * pi)/180;
            raw_x = (plot_x * pi)/180;
            means = [means_px py];
            raw = [raw_x' plot_y'];
            
            %fit with a distorted sin wave
            pars{mdepth_count} = sin_exp_fit(means,raw);
%             print_pars = 0;
%             if print_pars == 1
%                 PATHOUT = 'Z:\Users\jerry\SurfAnalysis\';
%                 FILEOUT = ['SinFitParamSummary.dat'];
%                 fileid = [PATHOUT FILEOUT];
%                 proffid = fopen(fileid, 'a');
%                 fprintf(proffid, '%s\t%3.2f\t%3.2f\t%1.4f\t%1.4f\t%4.2f\t%1.2f\n', FILE, unique_mean_depth(k), pars{k}(1), pars{k}(2), pars{k}(3), pars{k}(4), pars{k}(5));
%                 fclose(proffid);
%             end
            
            x_interp = (means_px(1)): .01 : (means_px(length(means_px)));
            y_sin = sin_exp_func(x_interp, pars{mdepth_count});
            y_err = sin_exp_err(pars{mdepth_count});
            y_sin(y_sin < 0) = 0;
            hold on
            x_rad = (x_interp * 180)/pi;
            size_rad = length(x_rad);
            plot(x_rad, y_sin, lines{mdepth_count});
            
            %run chi^2 test on fit
            num_free_params = 5;
            [chi2(mdepth_count), chiP(mdepth_count)] = Chi2_Test(raw_x', plot_y', 'sin_exp_func', pars{mdepth_count}, num_free_params);
            
            %store p-values of each mdisp curve
            p_val(stim_count,mdepth_count) = calc_mdisp_anovap(depth_select, plot_x, plot_y, unique_tilt);
            chiP_list(stim_count, mdepth_count) = chiP(mdepth_count);
            [value, index_max] = max(y_sin);
            pref_tilt(stim_count,mdepth_count) = x_rad(index_max);
            max_response(stim_count,mdepth_count) = value;
            
            null_x = [min(x_rad) max(x_rad)];
            null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
            null_y = [null_rate null_rate];
            
            printcurves = 0;
            if printcurves == 1
                %print out each individual tuning curve for origin
                for go=1:length(x_interp)
                    if isempty(curve_out{go})
                        curve_out{go} = '';
                    end
                    if (go<=2)
                        curve_out{go} = sprintf('%s%1.0f\t%3.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.3f\t%6.2f\t%6.2f\t', curve_out{go}, unique_stim(stim_count), unique_mean_depth(mdepth_count), x_interp(go), y_sin(go), px(go), py(go), perr(go), null_x(go),null_y(go));
                    elseif (go<=length(px))
                        curve_out{go} = sprintf('%s%1.0f\t%3.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.3f\t\t\t', curve_out{go}, unique_stim(stim_count), unique_mean_depth(mdepth_count), x_interp(go), y_sin(go), px(go), py(go), perr(go));
                    else
                        curve_out{go} = sprintf('%s%1.0f\t%3.2f\t%6.2f\t%6.2f\t\t\t\t\t\t', curve_out{go}, unique_stim(stim_count), unique_mean_depth(mdepth_count), x_interp(go), y_sin(go));
                    end
                end
                curve_out{1}
            end %end printcurves;
            
            [max_y, ind] = max(y_sin);
            [min_y, ind_min] = min(y_sin);
            
            angle_max(mdepth_count) = x_rad(ind);
            angle_min(mdepth_count) = x_rad(ind_min);
            
            max_resp(mdepth_count) = max_y;
            min_resp(mdepth_count) = min_y;
            
            %print out statistics for this mean depth
            stat_string = sprintf('\t%1.3f\t%1.0f\t%3.2f\t%1.4f\t%1.4f\t%3.2f\t%1.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%1.3f\t%1.4f\t', unique_mean_depth(mdepth_count), stim_count, pars{mdepth_count}(1), pars{mdepth_count}(2), pars{mdepth_count}(3), pars{mdepth_count}(4), pars{mdepth_count}(5), angle_max(mdepth_count), max_resp(mdepth_count), angle_min(mdepth_count), min_resp(mdepth_count), p_val(stim_count, mdepth_count), single_TDI(stim_count, mdepth_count));
            stat_out{mdepth_count+1} = stat_string;
            %--------------------------------------------------------------
            
            null_x = [min(px) max(px)];
            null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
            null_y = [null_rate null_rate];
            hold on;
            plot(null_x, null_y, 'k--');
%            yl = YLim;
%            YLim([0 yl(2)])	% set the lower limit of the Y axis to zero
            hold off;
            hold on
            
            figure(ver_graph)
            hold on;
            subplot(length(unique_stim), num_slant, stim_count);
            [ver_px, ver_py, ver_perr, ver_max, ver_min] = PlotTuningCurve(plot_x', ver', symbols{mdepth_count}, lines{mdepth_count}, 1, 1);
            p_val_ver(stim_count,mdepth_count) = calc_anovap(plot_x, ver);
            
            %save out each vergence curve so that we can 
            %mean shift them to calculate an avg TDI value  (same way as
            %neural data)
            start_verg(mdepth_count) = length(TDIvergdata)+1;
            stop_verg(mdepth_count) = length(plot_x)+start_verg(mdepth_count)-1;
            TDIvergdata(start_verg(mdepth_count):stop_verg(mdepth_count), 1) = plot_x';
            TDIvergdata(start_verg(mdepth_count):stop_verg(mdepth_count), 2) = ver';
            
            %------Verg ANOVAN code----------------------------------------
            %save out each data point to use in ANOVAN function
            verg(:,1) = plot_x';
            verg(:,2) = ver';
            sortedverg = sortrows(verg, [1]);
            verg_ancova_var(length(verg_ancova_var)+1:length(verg_ancova_var)+length(sortedverg),:) = sortedverg;
            for temp_tilt = 1:length(unique_tilt)
                tilt_ind = find(sortedverg(:,1) == unique_tilt(temp_tilt));
                sortedverg(tilt_ind(1):tilt_ind(length(tilt_ind)), 1) = temp_tilt;
            end
            %to do anovan
            verg_mdepth_array = zeros(length(sortedspikes), 1);
            verg_mdepth_array = mdepth_array + mdepth_count;
            verg_total_mdepth(length(verg_total_mdepth)+1:length(verg_total_mdepth)+length(sortedverg),:) = verg_mdepth_array;
            verg_list_angles(length(verg_list_angles)+1:length(verg_list_angles)+length(sortedverg),:) = sortedverg;
            %--------------------------------------------------------------
            
            printcurves_nonmodel = 0;
            if printcurves_nonmodel == 1
                %print out each individual tuning curve for origin
                pathsize = size(PATH,2) - 1;
                while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
                    pathsize = pathsize - 1;
                end   
                PATHOUT = 'Z:\Users\Jerry\SurfAnalysis\Surf_curves\';
                filesize = size(FILE,2) - 1;
                while FILE(filesize) ~='.'
                    filesize = filesize - 1;
                end
                FILEOUT = [FILE(1:filesize) 'surf_curve'];
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
        end %end mean depth loop
        
        %----ANOVAN--------------------------------------------------------
        list_angles = [total_mdepth list_angles];
        [p,T,STATS,TERMS] = anovan(list_angles(:, 3), {list_angles(:, 2) list_angles(:, 1)}, 'full', 3, {'Tilt Angles';'M. Depth'}, 'off');
        MS_error = T{4, 6};
        MS_treatment = T{2, 6};
        F_index = MS_treatment/ (MS_error + MS_treatment);
        
        verg_list_angles = [verg_total_mdepth verg_list_angles];
        [verg_p,verg_T,verg_STATS,verg_TERMS] = anovan(verg_list_angles(:, 3), {verg_list_angles(:, 2) verg_list_angles(:, 1)}, 'full', 3, {'Tilt Angles';'M. Depth'}, 'off');
        verg_MS_error = verg_T{4, 6};
        verg_MS_treatment = verg_T{2, 6};
        verg_F_index = verg_MS_treatment/ (verg_MS_error + verg_MS_treatment);
        %------------------------------------------------------------------

        %now that we have all the mean depth curves, mean shift the
        %responses
        %calc average TDI
        [avgTDI(stim_count), var_term] = compute_DDI(TDIdata(:,1)', TDIdata(:,2)');
        
        %readjust mean disparity responses to fall on the same mean
        %then calc avg TDI
        total_mean = mean(TDIdata(:,2));
        for count_meandepth = 1:length(unique_mean_depth)
            depth_mean = mean(TDIdata(start(count_meandepth):stop(count_meandepth),2));
            difference = total_mean - depth_mean;
            TDIdata(start(count_meandepth):stop(count_meandepth),2) = TDIdata(start(count_meandepth):stop(count_meandepth),2) + difference;
        end

        [avgTDI_adj(stim_count), var_term] = compute_DDI(TDIdata(:,1)', TDIdata(:,2)');
        [px, py, perr, spk_max, spk_min] = PlotTuningCurve(TDIdata(:,1), TDIdata(:,2), symbols{count_meandepth+1}, lines{count_meandepth+1}, 1, 0);
        Rmax_adj = spk_max;
        Rmin_adj = spk_min;
        
        zero_data = min(TDIvergdata(:,2));
        TDIvergdata(:,2) = TDIvergdata(:,2) + abs(zero_data);
        [avgTDIverg(stim_count), var_term] = compute_DDI(TDIvergdata(:,1)', TDIvergdata(:,2)');
        
        %readjust mean disparity responses to fall on the same mean
        %then calc avg vergence TDI
        total_mean = mean(TDIvergdata(:,2));
        for count_meandepth = 1:length(unique_mean_depth)
            depth_mean = mean(TDIvergdata(start_verg(count_meandepth):stop_verg(count_meandepth),2));
            difference = total_mean - depth_mean;
            TDIvergdata(start_verg(count_meandepth):stop_verg(count_meandepth),2) = TDIvergdata(start_verg(count_meandepth):stop_verg(count_meandepth),2) + difference;
        end
        
        [avgTDIverg_adj(stim_count), var_term] = compute_DDI(TDIvergdata(:,1)', TDIvergdata(:,2)');
        [px_verg, py_verg, perr_verg, verg_max, verg_min] = PlotTuningCurve(TDIvergdata(:,1), TDIvergdata(:,2), symbols{count_meandepth+1}, lines{count_meandepth+1}, 1, 0);
        
        figure(graph);
        subplot(length(unique_stim), num_slant, stim_count);
            
        yl = YLim;
        YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
        
        YLabel('Response (spikes/sec)');
        
        height = axis;
        yheight = height(4);
        if (unique_stim(stim_count) == 0)
            XLabel('Congruent Case - Tilt Angle (deg)');
        elseif (unique_stim(stim_count) == 1)
            XLabel('Speed Only Case - Tilt Angle (deg)');
        elseif (unique_stim(stim_count) == 2)
            XLabel('Disparity Only Case - Tilt Angle (deg)');
        elseif (unique_stim(stim_count) == 3)
            XLabel('Texture Only Case - Tilt Angle (deg)');
        end
            
        string = sprintf('TDI (mean adj) = %2.4f, Vergence TDI (mean adj) = %2.4f', avgTDI_adj(stim_count), avgTDIverg_adj(stim_count));
        text(height(1)+2, 0.95*yheight, string, 'FontSize', 8);
        string = sprintf('P-Values = %1.4f %1.4f %1.4f', p);
        text(height(1)+2, 0.85*yheight, string, 'FontSize', 8);
        
        figure(ver_graph);
        subplot(length(unique_stim), num_slant, stim_count);
        
        YLabel('Response (spikes/sec)');
        
        height = axis;
        yheight = height(4);
        if (unique_stim(stim_count) == 0)
            XLabel('Congruent Case - Tilt Angle (deg)');
        elseif (unique_stim(stim_count) == 1)
            XLabel('Speed Only Case - Tilt Angle (deg)');
        elseif (unique_stim(stim_count) == 2)
            XLabel('Disparity Only Case - Tilt Angle (deg)');
        elseif (unique_stim(stim_count) == 3)
            XLabel('Texture Only Case - Tilt Angle (deg)');
        end
        
        string = sprintf('P-Values = %1.4f %1.4f %1.4f', verg_p);
        text(height(1)+2, 0.95*yheight, string, 'FontSize', 8);        
        
        %store data for postprocessing (stimulus type, TDI_neuron,
        %TDI_vergence, max tilt_neuron, max_tilt_ver, spike p-val, verg p-val, spont)
        line = sprintf('\t%3.2f\t%1.5f\t%1.5f\t%3.2f\t%1.5f\t%1.5f\t%1.5f\t%3.2f', unique_stim(stim_count), avgTDI_adj(stim_count), avgTDIverg_adj(stim_count), spk_max.x, verg_max.x, p(1), verg_p(1), null_y(1));
        data_string = strcat(data_string, line);
        
    end %end slant angle
    printme = 1;
    if (printme==1)
        %pathsize = size(PATH,2) - 1;
        %while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
        %    pathsize = pathsize - 1;
        %end   
        PATHOUT = 'Z:\Users\jerry\SurfAnalysis\';
        
        line = sprintf('%s\t', FILE);
        %     data_string = strcat(line, data_string);
        %     
        %     %print grad metrics
        %     outfile = [PATHOUT 'Surf_metrics_07.16.04.dat'];
        %     printflag = 0;
        %     if (exist(outfile, 'file') == 0)    %file does not yet exist
        %         printflag = 1;
        %     end
        %     fid = fopen(outfile, 'a');
        %     if (printflag)
        %         %(Ap Size, X ctr, Y ctr, ecc, pref dir, stimulus type, TDI_neuron, TDI_vergence, max tilt_neuron, max_tilt_ver, spike p-val, verg p-val, spont)
        %         fprintf(fid, 'File\tApSize\tX\tY\tEcc\tPrefDir\tStimType\tConTDI\tTDIverg\tPrefTiltNeuron\tPrefTiltVerg\tPSpike\tPVerg\tSpont\tStimType1\tSPDTDI\tTDIverg1\tPrefTiltNeuron1\tPrefTiltVerg1\tPSpike1\tPVerg1\tSpont1\tStimType2\tDispTDI\tTDIverg2\tPrefTiltNeuron2\tPrefTiltVerg2\tPSpike2\tPVerg2\tSpont2\tStimType3\tTxtTDI\tTDIverg3\tPrefTiltNeuron3\tPrefTiltVerg3\tPSpike3\tPVerg3\tSpont3\n');
        %         printflag = 0;
        %     end
        %     fprintf(fid, '%s', [data_string]);
        %     fprintf(fid, '\r\n');
        %     fclose(fid);
        
        %print statistics for the individual tilt tuning curves
        %print grad metrics
        outfile = [PATHOUT 'Ind_SurfMetrics_09.02.04.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            stat_out{1} = sprintf('File\tMDepth\tStimType\tAmp\tFreq\tPhase\tBaseL\tExp\tPrefTilt\tTrough\tPval\tTDI\n');
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, '%s', [stat_out{1}]);
        end      

        for i=2:length(stat_out)
            stat_out{i} = strcat(line, stat_out{i});
            fprintf(fid, '%s', [stat_out{i}]);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
    end
end %end stim type

if printcurves == 1
    %print out each individual tuning curve for origin
    
    PATHOUT = 'Z:\Users\Jerry\SurfAnalysis\Surf_curves\';
    filesize = size(FILE,2) - 1;
    while FILE(filesize) ~='.'
        filesize = filesize - 1;
    end
    FILEOUT = [FILE(1:filesize) 'model_curve'];
    fileid = [PATHOUT FILEOUT];
    printflag = 0;
    if (exist(fileid, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    proffid = fopen(fileid, 'a');
    if (printflag)
        fprintf(proffid,'StimType\tMeanDepth\tIntHDisp\tSinFit\tHDisp\tAvgResp\tStdErr\tSpon\n');
        printflag = 0;
    end
    for go = 1:size_rad
        fprintf(proffid, '%s\n', curve_out{go});
    end
    fclose(proffid);
end

do_combo = 1;
if do_combo == 1
    
    %this next piece of code checks all tilt tuning curves with a
    %significance criteria of p < 0.05
    PATHOUT = 'Z:\Users\jerry\SurfAnalysis\';
%     outfile = [PATHOUT 'delta_tilt_combinations_withinStim_p0.05.dat'];
%     fid = fopen(outfile, 'a');
%     if length(unique_mean_depth) > 1
%         %go through each combination of mean disparities and print out the relationship between pref tilts
%         for i=1:length(unique_stim)
%             for j=1:length(unique_mean_depth)
%                 pref_out = '';
%                 if (p_val(i,j) < 0.05) %if significant
%                     %if(chiP_list(i,j) > .05) %if good fit
%                         for k=j+1:length(unique_mean_depth)
%                             if (p_val(i,k) < 0.05)
%                                 %if (chiP_list(i,k) > .05)
%                                     %print out pref tilt and mean disp for this combo
%                                     pref_out = sprintf('\t%1.0f\t%1.3f\t%3.2f\t%1.4f\t%1.0f\t%1.5f\t%3.2f\t%1.4f\t%3.2f', i, unique_mean_depth(j), pref_tilt(i, j), single_TDI(i, j), i, unique_mean_depth(k), pref_tilt(i, k), single_TDI(i, k), pref_tilt(i, j)-pref_tilt(i, k));
%                                     line = sprintf('%s', FILE);
%                                     pref_out = strcat(line, pref_out);
%                                     fprintf(fid, '%s', [pref_out]);
%                                     fprintf(fid, '\r');
%                                     %end %end 2nd chiP sig test
%                             end %end 2nd mdisp sig test
%                         end %end 2nd mdisp search
%                         %end %end 1st chiP sig test
%                 end %end 1st mdisp sig test
%             end %end 1st mdisp search
%         end
%         fclose(fid);
%     end %end check for multiple mdisp

    %this next piece of code goes through each mean depth and picks out the
    %preferred tilt for each stimulus type
    outfile = [PATHOUT 'delta_tilt_combo_acrossStim_p0.05.dat'];
    fid = fopen(outfile, 'a');
    if length(unique_mean_depth) > 1
        %go through each combination of mean disparities and print out the relationship between pref tilts
        for i=1:length(unique_mean_depth)
            for j=1:length(unique_stim)
                pref_out = '';
                if (p_val(j,i) < 0.05) %if significant
                    %if(chiP_list(i,j) > .05) %if good fit
                        for k=j+1:length(unique_stim)
                            if (p_val(k,i) < 0.05)
                                %if (chiP_list(i,k) > .05)
                                    %print out pref tilt and mean disp for this combo
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.2f\t%1.4f\t%1.0f\t%3.2f\t%1.4f\t%1.4f', unique_mean_depth(i), unique_stim(j), pref_tilt(j, i), single_TDI(j, i), unique_stim(k), pref_tilt(k, i), single_TDI(k, i));
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);
                                    fprintf(fid, '%s', [pref_out]);
                                    fprintf(fid, '\r');
                                    %end %end 2nd chiP sig test
                            end %end 2nd mdisp sig test
                        end %end 2nd mdisp search
                        %end %end 1st chiP sig test
                end %end 1st mdisp sig test
            end %end 1st mdisp search
        end
        fclose(fid);
    end %end check for multiple mdisp
    
end %end check if do_combo
return;