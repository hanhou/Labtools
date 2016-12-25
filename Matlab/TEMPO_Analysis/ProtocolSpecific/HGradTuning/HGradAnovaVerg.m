function ver_total = HGradAnovaVerg(ver_graph, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE)

TEMPO_Defs;


%Start Data Retrieval Routines---------------------------------------------------------------------------------------------------------
%get the column of values of horiz. disparity magnitude in the dots_params matrix
mag_disp = data.dots_params(DOTS_HGRAD_MAG,BegTrial:EndTrial,PATCH1);
   
%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (mag_disp == data.one_time_params(NULL_VALUE)) );

%get the column of mean disparity values
mean_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);

   
%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get the average horizontal eye positions to calculate vergence
Leyex_positions = data.eye_positions(1, :);
Reyex_positions = data.eye_positions(3, :);
   
vergence = Leyex_positions - Reyex_positions;

%get indices of monoc. and uncorrelated controls
control_trials = logical( (mean_disp == LEYE_CONTROL) | (mean_disp == REYE_CONTROL) | (mean_disp == UNCORR_CONTROL) );
   
   
unique_mag_disp = munique(mag_disp(~null_trials & ~control_trials)');	
   
%get the column of values of horiz. disparity angle of orientation in the dots_params matrix
disp_ang = data.dots_params(DOTS_HGRAD_ANGLE,BegTrial:EndTrial,PATCH1);
unique_disp_ang = munique(disp_ang(~null_trials & ~control_trials)');
   
%display monoc control switch
no_display_monoc = 1;
   
%display monoc or not?
if no_display_monoc == 1
    unique_mean_disp = munique(mean_disp(~null_trials & ~control_trials)');
else
    unique_mean_disp = munique(mean_disp(~null_trials)');
end
%unique_mean_disp = munique(mean_disp(~null_trials & ~control_trials)');
   
%get the column of different aperture sizes
ap_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH1);
unique_ap_size = munique(ap_size(~null_trials)');

all_sizes = 1;
if all_sizes ~= 1
    unique_ap_size = unique_ap_size(length(unique_ap_size));
    num_ap_size = length(unique_ap_size);
else
    num_ap_size = length(unique_ap_size);
end

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(mag_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
%End Data Retrieval Routines---------------------------------------------------------------------------------------------------------

%since the two way anova can only test one mean disparity for one set of gradient angles and only one grad magnitude
%we need to pull out each single data set then run the anove on them seperately

verg_arranged_by_angles = [];
verg_by_angles_by_mag_disp = [];
select_ap_size = 0;
do_anovan = 1;

num_mag_disp = length(unique_mag_disp);
num_ap_size = length(unique_ap_size);
   
%p_value_fig = figure;
data_string='';
ver_total = cell(length(unique_ap_size), 1);

   %if MU channel present then graph side-by-side with SU data
   plot_MU = 0;
   if (plot_MU)
       find_chan_num = size(data.spike_data);
       number_chans = find_chan_num(1);
       if ((number_chans > 2) & (SpikeChan < 2))
           length_MUchan = 2;
           ChanArray = [SpikeChan; number_chans];
       else
           length_MUchan = 1;
           ChanArray = [SpikeChan];
       end
   else
       length_MUchan = 1;
       ChanArray = [SpikeChan];   
   end

for curr_chan =1:length_MUchan
    for j = 1:length(unique_ap_size)
        ver_ancova_var = [];
        for k = 1:length(unique_mag_disp)
            list_angles = [];
            total_mdisp = [];
            for i=1:length(unique_mean_disp)	%for each different mean disparity, plot a separate disparity tuning curve
                verg = [];
                disp_select = logical((ap_size == unique_ap_size(j)) & (mag_disp == unique_mag_disp(k)) & (mean_disp == unique_mean_disp(i)) );
                if no_display_monoc == 1
                   plot_x = disp_ang(disp_select & ~null_trials & ~control_trials & select_trials);
                   plot_y = vergence(disp_select & ~null_trials & ~control_trials & select_trials);
                else
                   plot_x = disp_ang(disp_select & ~null_trials & select_trials);
                   plot_y = vergence(disp_select & ~null_trials & select_trials);
                end
                num_reps = length(find(plot_x==unique_disp_ang(i)));
                verg(:,1) = plot_x';
                verg(:,2) = plot_y';
                sortedverg = sortrows(verg, [1]);
                ver_ancova_var(length(ver_ancova_var)+1:length(ver_ancova_var)+length(sortedverg),:) = sortedverg;
                %to do anovan
                if do_anovan == 1
                   mdisp_array = zeros(length(sortedverg), 1);
                   mdisp_array = mdisp_array + i;
                   total_mdisp(length(total_mdisp)+1:length(total_mdisp)+length(sortedverg),:) = mdisp_array;
                   
                   list_angles(length(list_angles)+1:length(list_angles)+length(sortedverg),:) = sortedverg;
                else
                   for ik = 1:length(unique_disp_ang)
                      verg_arranged_by_angles(:,ik) = sortedverg((num_reps*(ik-1)+1):ik*num_reps,2);
                   end
                   verg_by_angles_by_mean_disp((num_reps*(i-1)+1):i*num_reps, :) = verg_arranged_by_angles;
                end
            end
            if do_anovan == 1
                   list_angles = [total_mdisp list_angles];
                   if i==1
                      p = anovan(list_angles(:, 3), {list_angles(:, 2)}, 'full', 3, {'Angles'}, 'off');
                   else
                      p = anovan(list_angles(:, 3), {list_angles(:, 2) list_angles(:, 1)}, 'full', 3, {'Angles';'M. Disp'}, 'off');
                   end
            else
                if i==1
                   p = anova1(verg_by_angles_by_mean_disp, 'off');
                   close(gcf);
                else
                   p = anova2(verg_by_angles_by_mean_disp, num_reps, 'off');
                end
            end
            
            line = sprintf('     %3.1d           %1.12f         %1.12f       %1.12f', unique_ap_size(j), p');
            data_string = strcat(data_string, line);
          
            figure(ver_graph);
            if(num_ap_size >= num_mag_disp)
                subplot((length_MUchan*num_ap_size), num_mag_disp,  (curr_chan-1)*(length_MUchan*num_ap_size-1) + (k-1)*(num_mag_disp) + j);
            elseif(num_ap_size < num_mag_disp)
                subplot((length_MUchan*num_mag_disp), num_ap_size,  (curr_chan-1)*(length_MUchan*num_mag_disp-1) + (k-1)*(num_ap_size) + j);
            end
            
            height = axis;
            yheight = height(4)-height(3);
            hold on
            string = sprintf('P-Values = %2.4f %2.4f %2.4f', p);
            text(height(1)+2, 0.75*yheight+height(3), string, 'FontSize', 8);
            hold off
        end
        ver_total{j} = ver_ancova_var;
    end
end %curr_chan

    printme = 0;
    if (printme == 1)
        %pathsize = size(PATH,2) - 1;
        %while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
        %    pathsize = pathsize - 1;
        %end   
        PATHOUT = 'Z:\Users\jerry\gradanalysis\figure_data\';           
        outfile = [PATHOUT 'Grad_anova_verg_vals_8.26.02.dat'];
    
        line = sprintf('%s', FILE);
        data_string = strcat(line, data_string);
    

        fid = fopen(outfile, 'a');
        fprintf(fid, '%s', [data_string]);
        fprintf(fid, '\r\n');
        fclose(fid);
    end