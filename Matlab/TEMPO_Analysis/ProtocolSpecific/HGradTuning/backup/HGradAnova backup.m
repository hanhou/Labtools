function ancova_total = HGradAnova(graph, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE)

	TEMPO_Defs;


%Start Data Retrieval Routines---------------------------------------------------------------------------------------------------------
	%get the column of values of horiz. disparity magnitude in the dots_params matrix
	mag_disp = data.dots_params(DOTS_HGRAD_MAG,BegTrial:EndTrial,PATCH1);
   
   %get indices of any NULL conditions (for measuring spontaneous activity)
   null_trials = logical( (mag_disp == data.one_time_params(NULL_VALUE)) );

   %get the column of mean disparity values
   mean_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);


   %get indices of monoc. and uncorrelated controls
   control_trials = logical( (mean_disp == LEYE_CONTROL) | (mean_disp == REYE_CONTROL) | (mean_disp == UNCORR_CONTROL) );
   unique_mean_disp = munique(mean_disp(~null_trials & ~control_trials)');
   
   %get the column of different aperture sizes
   ap_size = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH1);
   
   all_sizes = 0
   unique_ap_size = munique(ap_size(~null_trials & ~control_trials)');
   
   if all_sizes ~= 1
       unique_ap_size = unique_ap_size(length(unique_ap_size));
       num_ap_size = length(unique_ap_size);
   else
       num_ap_size = length(unique_ap_size);
   end
         
   unique_mag_disp = munique(mag_disp(~null_trials & ~control_trials)');	
   
 	%get the column of values of horiz. disparity angle of orientation in the dots_params matrix
   disp_ang = data.dots_params(DOTS_HGRAD_ANGLE,BegTrial:EndTrial,PATCH1);
   unique_disp_ang = munique(disp_ang(~null_trials & ~control_trials)');
   

   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:length(mag_disp);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
   %find how many reps there are for seperation of angles
   

   
%End Data Retrieval Routines---------------------------------------------------------------------------------------------------------

%since the two way anova can only test one mean disparity for one set of gradient angles and only one grad magnitude
%we need to pull out each single data set then run the anova on them seperately

    spikes_arranged_by_angles = [];
    spikes_by_angles_by_mag_disp = [];
    select_ap_size = 0;
    do_anovan = 1;
    
    data_string='';
    p_out = '';
    ancova_total = cell(length(unique_ap_size), 1);
       
    for j = 1:length(unique_ap_size)
        ancova_var = [];
        for k = 1:length(unique_mag_disp)
            list_angles = [];
            total_mdisp = [];
            for i=1:length(unique_mean_disp)	%for each different mean disparity, plot a separate tilt tuning curve
                spike = [];
                disp_select = logical((ap_size == unique_ap_size(j)) & (mag_disp == unique_mag_disp(k)) & (mean_disp == unique_mean_disp(i)) );
                plot_x = disp_ang(disp_select & ~null_trials & ~control_trials & select_trials);
                plot_y = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);
                num_reps = length(find(plot_x==unique_disp_ang(i)));
                spike(:,1) = plot_x';
                spike(:,2) = plot_y';
                sortedspikes = sortrows(spike, [1]);
                ancova_var(length(ancova_var)+1:length(ancova_var)+length(sortedspikes),:) = sortedspikes;
                for temp_ang = 1:length(unique_disp_ang)
                   ang_ind = find(sortedspikes(:,1) == unique_disp_ang(temp_ang));
                   sortedspikes(ang_ind(1):ang_ind(length(ang_ind)), 1) = temp_ang;
                end
                %to do anovan
                if do_anovan == 1
                   mdisp_array = zeros(length(sortedspikes), 1);
                   mdisp_array = mdisp_array + i;
                   total_mdisp(length(total_mdisp)+1:length(total_mdisp)+length(sortedspikes),:) = mdisp_array;
                   list_angles(length(list_angles)+1:length(list_angles)+length(sortedspikes),:) = sortedspikes;
                else
                   for ik = 1:length(unique_disp_ang)
    	               spikes_arranged_by_angles(:,ik) = sortedspikes((num_reps*(ik-1)+1):ik*num_reps,2);
                   end
                   spikes_by_angles_by_mean_disp((num_reps*(i-1)+1):i*num_reps, :) = spikes_arranged_by_angles;
                end
                
                p = anovan(sortedspikes(:, 2), {sortedspikes(:, 1)}, 'full', 3, {'Angles'}, 'off');
                p_string = sprintf(' %1.3f %1.4f', unique_mean_disp(i), p);
                p_out = strcat(p_out, p_string);
               
            end %end m.disp
            if do_anovan == 1
               list_angles = [total_mdisp list_angles];
               if i==1
                  [p,T,STATS,TERMS] = anovan(list_angles(:, 3), {list_angles(:, 2)}, 'full', 3, {'Angles'}, 'off');
               else
                  [p,T,STATS,TERMS] = anovan(list_angles(:, 3), {list_angles(:, 2) list_angles(:, 1)}, 'full', 3, {'Angles';'M. Disp'}, 'off');
               end
               MS_error = T{4, 6};
               MS_treatment = T{2, 6};
               F_index = MS_treatment/ (MS_error + MS_treatment);
            else
               if i==1
                  p = anova1(spikes_by_angles_by_mean_disp, 'off');
                  %close(gcf);
               else
                  p = anova2(spikes_by_angles_by_mean_disp, num_reps, 'off')
               end
            end
            
            
            
            line = sprintf('\t%3.1d\t%1.12f\t%1.12f\t%1.12f\t%1.5f', unique_ap_size(j), p', F_index);
            data_string = strcat(data_string, line);


            %close(gcf);
            
            %save list_angle.txt list_angles -ASCII
            %fid = fopen('list_angle.txt', 'a');
            %fprintf(fid, '%d %d %d \r', list_angles(:,1), list_angles(:, 2), list_angles(:, 3));
            %fclose(fid);
         
            figure(graph)
            if(length(unique_ap_size) >= length(unique_mag_disp))
                subplot(length(unique_ap_size), length(unique_mag_disp),  (k-1)*(length(unique_mag_disp)) + j);
            elseif(length(unique_ap_size) < length(unique_mag_disp))
                subplot(length(unique_mag_disp), length(unique_ap_size),  (k-1)*(length(unique_ap_size)) + j);
            end
            height = axis;
            yheight = height(4);
            hold on
            string = sprintf('P-Values = %2.4f %2.4f %2.4f', p);
            text(height(1)+2, 0.6*yheight, string, 'FontSize', 8);
            hold off;
        end %end mag_disp
        ancova_total{j} = ancova_var;
    end %end ap_size
    
    printme = 0;
    if (printme == 1)
        %pathsize = size(PATH,2) - 1;
        %while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
        %    pathsize = pathsize - 1;
        %end   
        PATHOUT = 'Z:\Users\Jerry\GradAnalysis\figure_data\';           
        outfile = [PATHOUT 'SU_Grad_anova_vals.dat'];
    
        line = sprintf('%s', FILE);
        data_string = strcat(line, data_string);
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
           printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
           fprintf(fid,'File\tApSize1\tP1\tP2\tP3\tF1\tApSize2\tP1two\tP2two\tP3two\tF1two\n')
           printflag = 0;
        end
        
        if (length(unique_ap_size) == 2)
           fprintf(fid, '%s', [data_string]);
        else
           fprintf(fid, '%s\t\t\t\t\t', [data_string]);
        end
        fprintf(fid, '\r\n');
        fclose(fid);
        
        outfile = [PATHOUT 'SU_MDisp_pvals.dat'];
        line = sprintf('%s', FILE);
        p_out = strcat(line, p_out);
        fid= fopen(outfile, 'a');
        fprintf(fid, '%s', [p_out]);
        fprintf(fid, '\r\n');
        fclose(fid);        
    end
   
   