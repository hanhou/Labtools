%-----------------------------------------------------------------------------------------------------------------------
%-- SurfTuningCurve.m -- Plots a horizontal disparity gradient tuning curve.  These tuning curves will plot
%--   varying angles of gradient rotation vs. responses for different mean disparities on a single graph.  multiple
%--   graphs in a single column represent different gradient magnitudes for one aperture size.  Graphs in 
%--   different columns differ by aperture size.  All graphs include	monoc and uncorrelated control conditions.
%--	JDN 8/07/04
%-----------------------------------------------------------------------------------------------------------------------
function SurfVectorAvg(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};
color_dots = {'b.' 'r.' 'g.' 'k.'};
color_lines = {'b*' 'r*' 'g*' 'k*'};

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
shift_negativetilt = logical(tilt_list < 0);
shift_positivetilt = logical(tilt_list > 360);

tilt_list(shift_negativetilt) = tilt_list(shift_negativetilt) + 360;
tilt_list(shift_positivetilt) = tilt_list(shift_positivetilt) - 360;
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
pol_graph = figure;
% boot_graph = figure;

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

plot_x_all = [];
plot_y_all = [];
x_len = [];
y_len = [];
for stim_count = 1:length(unique_stim)
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
        
        %extract all tilts and responses for vector average calculation
        depth_select = logical((stim_list == unique_stim(stim_count)) &(slant_list == unique_slant(slant_count)));
        x_len(stim_count) = length(tilt_list(depth_select & ~null_trials & select_trials));
        plot_x_all(stim_count, 1:x_len(stim_count)) = tilt_list(depth_select & ~null_trials & select_trials);
        y_len(stim_count) = length(spike_rates(depth_select & ~null_trials & select_trials))
        plot_y_all(stim_count, 1:y_len(stim_count)) = spike_rates(depth_select & ~null_trials & select_trials);
        
        %get null values
        null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
        null_y = null_rate;
        
        %calculate vector average
        rad_x_all(stim_count, 1:x_len(stim_count)) = plot_x_all(stim_count, 1:x_len(stim_count)) * 3.14159/180;
        [x_comp, y_comp] = pol2cart(rad_x_all(stim_count, 1:x_len), plot_y_all(stim_count, 1:y_len));
        sum_x = sum(x_comp);
        sum_y = sum(y_comp);
        %mag = sum(sqrt(x_comp.^2+y_comp.^2))
        mag = sum(plot_y_all(stim_count, 1:y_len(stim_count)))

        vect_xavg(stim_count) = sum_x/(mag);
        vect_yavg(stim_count) = sum_y/(mag);
        [th_avg(stim_count), r_avg(stim_count)] = cart2pol(vect_xavg(stim_count), vect_yavg(stim_count));
        
        for mdepth_count=1:length(unique_mean_depth)
            spike = [];
            verg = [];
                
            depth_select = logical((stim_list == unique_stim(stim_count)) &(slant_list == unique_slant(slant_count)) & (mean_depth_list == unique_mean_depth(mdepth_count)));
            plot_x = tilt_list(depth_select & ~null_trials & select_trials);
            plot_y = spike_rates(depth_select & ~null_trials & select_trials);
            ver = vergence(depth_select & ~null_trials & select_trials);
            
            %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
            [px, py, perr, spk_max, spk_min] = PlotTuningCurve(plot_x', plot_y', symbols{mdepth_count}, lines{mdepth_count}, 0, 0);
            
            null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
            null_y = null_rate;
            
            %calculate vector average
            rad_x = plot_x * 3.14159/180;
            [x_comp, y_comp] = pol2cart(rad_x, plot_y);
            sum_x = sum(x_comp);
            sum_y = sum(y_comp);
            %mag = sum(sqrt(x_comp.^2+y_comp.^2))
            mag = sum(plot_y);

            vector_xavg(stim_count,mdepth_count) = sum_x/(mag);
            vector_yavg(stim_count,mdepth_count) = sum_y/(mag);
            
            p_val(stim_count,mdepth_count) = calc_anovap(plot_x, plot_y);
           
            
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

    end %end slant angle
end %end stim type

figure(pol_graph);
for stim_count = 1:length(unique_stim)
    [th(stim_count,1:3), r(stim_count,1:3)] = cart2pol(vector_xavg(stim_count:4:12), vector_yavg(stim_count:4:12));
    polar(th(stim_count, 1:3), r(stim_count, 1:3), color_dots{stim_count});
    hold on
    polar(th_avg(stim_count), r_avg(stim_count), color_lines{stim_count});    
end

%perform permutation test
th_perm = zeros(1000, length(unique_stim));
r_perm = zeros(1000, length(unique_stim));
for rep_num = 1:1000
    for i=1:length(unique_stim)
        perm_y = randperm(length(plot_y_all(i, 1:y_len(i))));
        permuted_y = plot_y_all(i, perm_y);

        [x_comp, y_comp] = pol2cart(rad_x_all(i, 1:x_len(i)), permuted_y);
        sum_x = sum(x_comp);
        sum_y = sum(y_comp);
        mag = sum(sqrt(x_comp.^2+y_comp.^2));

        vect_xavg = sum_x/(mag);
        vect_yavg = sum_y/(mag);
        [th_perm(rep_num, i), r_perm(rep_num, i)] = cart2pol(vect_xavg, vect_yavg);
        for k=1:length(unique_mean_depth)
            depth_select = logical((stim_list == unique_stim(i)) & (mean_depth_list == unique_mean_depth(k)));
            plot_x = tilt_list(depth_select & ~null_trials & select_trials);
            plot_y = spike_rates(depth_select & ~null_trials & select_trials);
            %mean depth bootstrap loop
            list_count = 1;
            for tilt_count = 1:length(unique_tilt)
                tilt_ind = find(plot_x==unique_tilt(tilt_count));
                for tilties = 1:length(tilt_ind)
                    rand_pt = randperm(length(tilt_ind));
                    bootstrap_x(list_count) = unique_tilt(tilt_count);
                    bootstrap_y(list_count) = plot_y(tilt_ind(rand_pt(1)));
                    list_count = list_count+1;
                end
            end
            %calculate bootstrap vector average
            rad_x = bootstrap_x * 3.14159/180;
            
            [x_comp, y_comp] = pol2cart(rad_x, bootstrap_y);
            sum_x = sum(x_comp);
            sum_y = sum(y_comp);
            mag = sum(sqrt(x_comp.^2+y_comp.^2));
            boot_xavg(i,k) = sum_x/(mag);
            boot_yavg(i,k) = sum_y/(mag);
            
            
            %mean depth permutation loop
            perm_y_mdpth = randperm(length(plot_y));
            permuted_y_mdpth = plot_y(perm_y_mdpth);

            rad_x = plot_x * 3.14159/180;
            [x_comp, y_comp] = pol2cart(rad_x, plot_y(perm_y_mdpth));
            sum_x = sum(x_comp);
            sum_y = sum(y_comp);
            mag = sum(sqrt(x_comp.^2+y_comp.^2));
            permdepth_xavg(i,k) = sum_x/(mag);
            permdepth_yavg(i,k) = sum_y/(mag);
        end
        [th_md_perm(i,1:3), r_md_perm(i,1:3)] = cart2pol(permdepth_xavg(i,1:3), permdepth_yavg(i,1:3));
        [boot_th(i,1:3), boot_r(i,1:3)] = cart2pol(boot_xavg(i,1:3), boot_yavg(i,1:3));
        
        th_m1(rep_num, i) = th_md_perm(i, 1);
        th_m2(rep_num, i) = th_md_perm(i, 2);
        th_m3(rep_num, i) = th_md_perm(i, 3);
        
        r_m1(rep_num, i) = r_md_perm(i, 1);
        r_m2(rep_num, i) = r_md_perm(i, 2);
        r_m3(rep_num, i) = r_md_perm(i, 3);        
        
        boot_th_m1(rep_num, i) = boot_th(i, 1);
        boot_th_m2(rep_num, i) = boot_th(i, 2);
        boot_th_m3(rep_num, i) = boot_th(i, 3);
        
        boot_r_m1(rep_num, i) = boot_r(i, 1);
        boot_r_m2(rep_num, i) = boot_r(i, 2);
        boot_r_m3(rep_num, i) = boot_r(i, 3);
    end
end

for i=1:length(unique_stim)
    %get 95% CI for bootstraps
    sorted_r_m1 = sort(boot_r_m1(:, i));
    sorted_r_m2 = sort(boot_r_m2(:, i));
    sorted_r_m3 = sort(boot_r_m3(:, i));
    
    sorted_th_m1 = sort(boot_th_m1(:,i));
    sorted_th_m2 = sort(boot_th_m2(:,i));
    sorted_th_m3 = sort(boot_th_m3(:,i));
    
    thCI_95_m1_U(i) = sorted_th_m1(975);
    thCI_95_m1_L(i) = sorted_th_m1(25);
    thCI_95_m2_U(i) = sorted_th_m2(975);
    thCI_95_m2_L(i) = sorted_th_m2(25);
    thCI_95_m3_U(i) = sorted_th_m3(975);
    thCI_95_m3_L(i) = sorted_th_m3(25);
    
    CI_95_m1_U(i) = sorted_r_m1(975);
    CI_95_m1_L(i) = sorted_r_m1(25);
    CI_95_m2_U(i) = sorted_r_m2(975);
    CI_95_m2_L(i) = sorted_r_m2(25);
    CI_95_m3_U(i) = sorted_r_m3(975);
    CI_95_m3_L(i) = sorted_r_m3(25);
    
    list = find(r_perm(:,i) >= r_avg(i));
    sig = length(list);
    perm_p_val(i) = sig/length(r_perm(:,i));
    
    list = find(r_m1(:,i) >= r(i, 1));
    sig = length(list);
    m1_p_val(i) = sig/length(r_m1(:,i));
    
    list = find(r_m2(:,i) >= r(i, 2));
    sig = length(list);
    m2_p_val(i) = sig/length(r_m2(:,i));
    
    list = find(r_m3(:,i) >= r(i, 3));
    sig = length(list);
    m3_p_val(i) = sig/length(r_m3(:,i));
    
end

m_pval = [m1_p_val; m2_p_val; m3_p_val]';
CI_95_U = [CI_95_m1_U; CI_95_m2_U; CI_95_m3_U]';
CI_95_L = [CI_95_m1_L; CI_95_m2_L; CI_95_m3_L]';

thCI_95_U = [thCI_95_m1_U; thCI_95_m2_U; thCI_95_m3_U]';
thCI_95_L = [thCI_95_m1_L; thCI_95_m2_L; thCI_95_m3_L]';

stim_combo = 1;
if(stim_combo)
    PATHOUT = 'Z:\Users\jerry\SurfAnalysis\StimDeltaTilts\';
%     PATHOUT = 'Z:\Users\jerry\SurfAnalysis\';
    outfile = [PATHOUT 'AllSigConvsSpdCombo_12.01.04.dat'];
    fid_ConvsSpd = fopen(outfile, 'a');
    outfile = [PATHOUT 'AllSigConvsDispCombo_12.01.04.dat'];
    fid_ConvsDisp = fopen(outfile, 'a');
    outfile = [PATHOUT 'AllSigConvsTxtCombo_12.01.04.dat'];
    fid_ConvsTxt = fopen(outfile, 'a');
    
    outfile = [PATHOUT 'AllSigSpdvsDispCombo_12.01.04.dat'];
    fid_SpdvsDisp = fopen(outfile, 'a');
    outfile = [PATHOUT 'AllSigSpdvsTxtCombo_12.01.04.dat'];
    fid_SpdvsTxt = fopen(outfile, 'a');
    
    outfile = [PATHOUT 'AllSigDispvsTxtCombo_12.01.04.dat'];
    fid_DispvsTxt = fopen(outfile, 'a');
    
    for meancount = 1:length(unique_mean_depth)
        [sortedU, UI] = sort(thCI_95_U(:, meancount));
        [sortedR, RI] = sort(CI_95_U(:, meancount));
        for i=4:-1:2
            if (length(find(m_pval(1:3,meancount)<0.05)) == 3)
                if (m_pval(UI(i),meancount) < 0.05)
                    for j=i-1:-1:1
                        if (m_pval(UI(j),meancount) < 0.05)
                            %checks to see if upper confidence interval of
                            %pref tilt is greater than lower confidence
                            %interval of next highest tilt
                            if (sortedU(j)>thCI_95_L(UI(i), meancount))
                                sig = 0;
                            else
                                sig = 1;
                            end
                            
                            %checks to see if upper confidence interval of
                            %pref tilt is greater than lower confidence
                            %interval of next highest tilt                            
                            if (sortedR(j)>CI_95_L(RI(i), meancount))
                                r_sig = 0;
                            else
                                r_sig = 1;
                            end
                            %calculate adjusted angles such that the angle
                            %one falls between 180 to -180
                            pref_one = th(UI(i), meancount)*180/3.14159;
                            pref_two = th(UI(j), meancount)*180/3.14159;
                            if pref_one > 360
                                pref_adj_one = mod(pref_one, 360);
                            elseif pref_one < -360
                                pref_adj_one = mod(pref_one, -360);
                            end
                            
                            if pref_one < 0
                                if pref_one < -180
                                    pref_adj_one = mod(pref_one, 180);
                                else
                                    pref_adj_one = mod(pref_one, -180);
                                end
                            elseif pref_one > 0
                                if pref_one > 180
                                    pref_adj_one = mod(pref_one, -180);
                                else
                                    pref_adj_one = mod(pref_one, 180);
                                end
                            end
                            
                            test1 = mod(pref_two, 360);
                            test2 = mod(pref_two, -360);
                            
                            if abs(pref_adj_one-test1) < abs(pref_adj_one-test2)
                                off_adj = test1;
                            else
                                off_adj = test2;
                            end
                            
                            delta_tilt = abs(pref_adj_one - off_adj);
                            
                            %sort output based on stimulus type combination
                            if (unique_stim(UI(i)) == 0)
                                pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), sig, r_sig, pref_adj_one, off_adj, delta_tilt);
                                line = sprintf('%s', FILE);
                                pref_out = strcat(line, pref_out);
                                if (unique_stim(UI(j)) == 1)
                                    fprintf(fid_ConvsSpd, '%s', [pref_out]);
                                    fprintf(fid_ConvsSpd, '\r');
                                elseif (unique_stim(UI(j)) == 2)
                                    fprintf(fid_ConvsDisp, '%s', [pref_out]);
                                    fprintf(fid_ConvsDisp, '\r');
                                elseif (unique_stim(UI(j)) == 3)
                                    fprintf(fid_ConvsTxt, '%s', [pref_out]);
                                    fprintf(fid_ConvsTxt, '\r');
                                end
                            elseif (unique_stim(UI(i)) == 1)
                                if (unique_stim(UI(j)) == 0)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), sig, r_sig, off_adj, pref_adj_one, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);
                                    fprintf(fid_ConvsSpd, '%s', [pref_out]);
                                    fprintf(fid_ConvsSpd, '\r');
                                elseif (unique_stim(UI(j)) == 2)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), sig, r_sig, pref_adj_one, off_adj, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);                                    
                                    fprintf(fid_SpdvsDisp, '%s', [pref_out]);
                                    fprintf(fid_SpdvsDisp, '\r');
                                elseif (unique_stim(UI(j)) == 3)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), sig, r_sig, pref_adj_one, off_adj, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);                                    
                                    fprintf(fid_SpdvsTxt, '%s', [pref_out]);
                                    fprintf(fid_SpdvsTxt, '\r');
                                end
                            elseif (unique_stim(UI(i)) == 2)
                                if (unique_stim(UI(j)) == 0)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), sig, r_sig, off_adj, pref_adj_one, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);
                                    fprintf(fid_ConvsDisp, '%s', [pref_out]);
                                    fprintf(fid_ConvsDisp, '\r');
                                elseif (unique_stim(UI(j)) == 1)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), sig, r_sig, off_adj, pref_adj_one, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);
                                    fprintf(fid_SpdvsDisp, '%s', [pref_out]);
                                    fprintf(fid_SpdvsDisp, '\r');
                                elseif (unique_stim(UI(j)) == 3)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), sig, r_sig, pref_adj_one, off_adj, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);                                    
                                    fprintf(fid_DispvsTxt, '%s', [pref_out]);
                                    fprintf(fid_DispvsTxt, '\r');
                                end
                            elseif (unique_stim(UI(i)) == 3)
                                if (unique_stim(UI(j)) == 0)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), sig, r_sig, off_adj, pref_adj_one, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);                                    
                                    fprintf(fid_ConvsTxt, '%s', [pref_out]);
                                    fprintf(fid_ConvsTxt, '\r');
                                elseif (unique_stim(UI(j)) == 1)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), sig, r_sig, off_adj, pref_adj_one, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);                                    
                                    fprintf(fid_SpdvsTxt, '%s', [pref_out]);
                                    fprintf(fid_SpdvsTxt, '\r');
                                elseif (unique_stim(UI(j)) == 2)
                                    pref_out = sprintf('\t%3.2f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%3.4f\t%1.4f\t%1.0f\t%1.0f\t%3.4f\t%3.4f\t%3.4f', unique_mean_depth(meancount), unique_stim(UI(j)), th(UI(j), meancount)*180/3.14159, r(UI(j), meancount), unique_stim(UI(i)), th(UI(i), meancount)*180/3.14159, r(UI(i), meancount), sig, r_sig, off_adj, pref_adj_one, delta_tilt);
                                    line = sprintf('%s', FILE);
                                    pref_out = strcat(line, pref_out);                                    
                                    fprintf(fid_DispvsTxt, '%s', [pref_out]);
                                    fprintf(fid_DispvsTxt, '\r');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    fclose(fid_ConvsSpd);
    fclose(fid_ConvsDisp);
    fclose(fid_ConvsTxt);
    fclose(fid_SpdvsDisp);
    fclose(fid_SpdvsTxt);
    fclose(fid_DispvsTxt);
end

rho_plot = figure;
for i=1:length(unique_mean_depth)
    errorbar(unique_stim, r(:,i), CI_95_U(:,i)-r(:,i), r(:,i)-CI_95_L(:,i), symbols{i})
    hold on
end

print_CI = 1;
if(print_CI)
    PATHOUT = 'Z:\Users\jerry\SurfAnalysis\';
    outfile = [PATHOUT 'BootStrapCI_10.11.04.dat'];
    fid = fopen(outfile, 'a');
    for i=1:length(unique_stim)
        pref_out = sprintf('\t%1.0f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%3.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t', unique_stim(i), r(i,1), CI_95_U(i,1)-CI_95_L(i,1), r(i,2), CI_95_U(i,2)-CI_95_L(i,2), r(i,3), CI_95_U(i,3)-CI_95_L(i,3), th(i,1)*180/3.14159, thCI_95_U(i,1)-thCI_95_L(i,1), th(i,2)*180/3.14159, thCI_95_U(i,2)-thCI_95_L(i,2), th(i,3)*180/3.14159, thCI_95_U(i,3)-thCI_95_L(i,3));
        line = sprintf('%s', FILE);
        pref_out = strcat(line, pref_out);
        fprintf(fid, '%s', [pref_out]);
        fprintf(fid, '\r')        
    end
    fclose(fid);
end

print_vectors = 0;
if (print_vectors)
    PATHOUT = 'Z:\Users\jerry\SurfAnalysis\';
    outfile = [PATHOUT 'VectorAvgCon_MDepthSummary_12.01.04.dat'];
    fidCon = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgSpd_MDepthSummary_12.01.04.dat'];
    fidSpd = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgDisp_MDepthSummary_12.01.04.dat'];
    fidDisp = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgTxt_MDepthSummary_12.01.04.dat'];
    fidTxt = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgCon_Summary_12.01.04.dat'];
    fidConAvg = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgSpd_Summary_12.01.04.dat'];
    fidSpdAvg = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgDisp_Summary_12.01.04.dat'];
    fidDispAvg = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgTxt_Summary_12.01.04.dat'];
    fidTxtAvg = fopen(outfile, 'a'); 
    outfile = [PATHOUT 'VectorAvgSummary_12.01.04.dat'];
    fid2 = fopen(outfile, 'a');
    for i=1:length(unique_stim)
        for j=1:length(unique_mean_depth)
            pref_out = sprintf('\t%1.0f\t%3.2f\t%3.4f\t%1.4f\t%1.4f\t', unique_stim(i), unique_mean_depth(j), th(i, j)*180/3.14159, r(i,j), m_pval(i,j));
            line = sprintf('%s', FILE);
            pref_out = strcat(line, pref_out);
            if (i==1)
                fprintf(fidCon, '%s', [pref_out]);
                fprintf(fidCon, '\r');
            elseif (i==2)
                fprintf(fidSpd, '%s', [pref_out]);
                fprintf(fidSpd, '\r');
            elseif (i==3)
                fprintf(fidDisp, '%s', [pref_out]);
                fprintf(fidDisp, '\r');
            elseif (i==4)
                fprintf(fidTxt, '%s', [pref_out]);
                fprintf(fidTxt, '\r');
            end
           
        end
        pref_out = sprintf('\t%1.0f\t%3.2f\t%3.4f\t%1.4f\t%1.4\t', unique_stim(i), th_avg(i)*180/3.14159, r_avg(i), perm_p_val(i));
        line = sprintf('%s', FILE);
        pref_out = strcat(line, pref_out);
        if (i==1)
            fprintf(fidConAvg, '%s', [pref_out]);
            fprintf(fidConAvg, '\r');
        elseif (i==2)
            fprintf(fidSpdAvg, '%s', [pref_out]);
            fprintf(fidSpdAvg, '\r');
        elseif (i==3)
            fprintf(fidDispAvg, '%s', [pref_out]);
            fprintf(fidDispAvg, '\r');
        elseif (i==4)
            fprintf(fidTxtAvg, '%s', [pref_out]);
            fprintf(fidTxtAvg, '\r');
        end
    end
    fclose(fidCon);
    fclose(fidSpd);
    fclose(fidDisp);
    fclose(fidTxt);
    fclose(fidConAvg);
    fclose(fidSpdAvg);
    fclose(fidDispAvg);
    fclose(fidTxtAvg);
end
return;