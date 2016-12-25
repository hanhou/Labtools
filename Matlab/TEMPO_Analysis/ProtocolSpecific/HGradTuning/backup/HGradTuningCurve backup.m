%-----------------------------------------------------------------------------------------------------------------------
%-- HGradTuningCurve.m -- Plots a horizontal disparity gradient tuning curve.  These tuning curves will plot
%--   varying angles of gradient rotation vs. responses for different mean disparities on a single graph.  multiple
%--   graphs in a single column represent different gradient magnitudes for one aperture size.  Graphs in 
%--   different columns differ by aperture size.  All graphs include	monoc and uncorrelated control conditions.
%--	JDN 2/4/00
%-----------------------------------------------------------------------------------------------------------------------
function HGradTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
   display_monoc = 1;
   
   %display monoc or not?
   if display_monoc == 1
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
       
   
   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);
   
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
   graph = figure;
   set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [500 50 1000 773], 'Name', 'Tilt Tuning Curve');
   ver_graph = figure;
   set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 50 500 773], 'Name', 'Tilt vs. Horizontal Vergence');
   
   data_string = '';   
   angle_out = '';
   GDI_mdisp_out = '';
   IndTDI_vs_ALLTDI = '';
   %resp_data = []; verg_data=[];
   p_val = zeros(length(unique_ap_size), length(unique_mean_disp));
   pref_tilt = zeros(length(unique_ap_size), length(unique_mean_disp));
   for j=1:length(unique_ap_size)
      GDIdata = [];
      GDIvergdata = [];
      for k=1:length(unique_mag_disp)
         Gmean=zeros(length(unique_disp_ang), 2);
         Vmean=zeros(length(unique_disp_ang), 2);
         m_disp_max = zeros(length(unique_mean_disp), 1);
         m_disp_min = zeros(length(unique_mean_disp), 1);
         start = zeros(length(unique_mean_disp), 1);
         stop = zeros(length(unique_mean_disp), 1);
         start_verg = zeros(length(unique_mean_disp), 1);
         stop_verg = zeros(length(unique_mean_disp), 1);
         for l=1:length(unique_mean_disp)
            figure(graph);
            hold on;
            if(num_ap_size >= num_mag_disp)
               subplot(num_ap_size, num_mag_disp,  (k-1)*(num_mag_disp) + j);
            elseif(num_ap_size < num_mag_disp)
               subplot(num_mag_disp, num_ap_size,  (k-1)*(num_ap_size) + j);
            end
            disp_select = logical((ap_size == unique_ap_size(j)) & (mag_disp == unique_mag_disp(k)) & (mean_disp == unique_mean_disp(l)) );
      
            if display_monoc == 1
                plot_x = disp_ang(disp_select & ~null_trials & ~control_trials & select_trials);
                plot_y = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);
                ver = vergence(disp_select & ~null_trials & ~control_trials & select_trials);
            else
                plot_x = disp_ang(disp_select & ~null_trials & select_trials);
                plot_y = spike_rates(disp_select & ~null_trials & select_trials);
                ver = vergence(disp_select & ~null_trials & select_trials);
            end
            
            %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
            [px, py, perr, spk_max, spk_min] = PlotTuningCurve(plot_x', plot_y', symbols{l}, lines{l}, 1, 1);
            m_disp_max(l) = spk_max.x;
            m_disp_min(l) = spk_min.x;
            null_x = [min(px) max(px)];
            null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
            null_y = [null_rate null_rate];

            p_val(j,l) = calc_mdisp_anovap(disp_select, plot_x, plot_y, unique_disp_ang);
            [value, index_max] = max(py);
            pref_tilt(j,l) = px(index_max);
            
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
               FILEOUT = [FILE(1:filesize) 'grad_curve2'];
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

             
                
            Gmean(:, 1) = Gmean(:,1) + px;
            Gmean(:, 2) = Gmean(:,2) + py;
                
            start(l) = length(GDIdata)+1;
            stop(l) = length(plot_x)+start(l)-1;
            GDIdata(start(l):stop(l), 1) = plot_x';
            GDIdata(start(l):stop(l), 2) = plot_y';
                
            %compute a GDI using plot_x and plot_y
            GDI(l,2) = unique_mean_disp(l);
            [GDI(l,1), var_term] = Compute_DDI(plot_x, plot_y);
                
            hold off;
                
            figure(ver_graph)
            hold on
            if(num_ap_size >= num_mag_disp)
               subplot(num_ap_size, num_mag_disp,  (k-1)*(num_mag_disp) + j);
            elseif(num_ap_size < num_mag_disp)
               subplot(num_mag_disp, num_ap_size,  (k-1)*(num_ap_size) + j);
            end

            hold on;
            [vx, vy, verr, temp_max, temp_min] = PlotTuningCurve(plot_x', ver', symbols{l}, lines{l}, 0, 1);
            Vmean(:, 1) = Vmean(:,1) + vx;
            Vmean(:, 2) = Vmean(:,2) + vy;
            hold off;
                
            start_verg(l) = length(GDIvergdata)+1;
            stop_verg(l) = length(plot_x)+start_verg(l)-1;
            GDIvergdata(start_verg(l):stop_verg(l), 1) = plot_x';
            GDIvergdata(start_verg(l):stop_verg(l), 2) = ver';
                
            leg{l} = sprintf('%1.3f', unique_mean_disp(l));
            angle_string = sprintf('\t%1.3f\t%3.2f', unique_mean_disp(l), m_disp_max(l));
            angle_out = strcat(angle_out, angle_string);
            
            GDI_mdisp_string = sprintf('\t%1.3f\t%1.4f', unique_mean_disp(l), GDI(l,1));
            GDI_mdisp_out = strcat(GDI_mdisp_out, GDI_mdisp_string);
         end %end m_disp loop
          
         figure(graph);
         axhandles = get(gca, 'Children');
         numlines = length(unique_mean_disp);
         leghandles = [];
         legindex = 1;
         for legendloop = numlines:-1:1
            leghandles(legindex) = axhandles((legendloop*3));
            legindex = legindex + 1;
         end
           
         legend(leghandles, leg, 4);
           
         figure(ver_graph);
         axhandles = get(gca, 'Children');
         numlines = length(unique_mean_disp);
         leghandles = [];
         legindex = 1;
         for legendloop = numlines:-1:1
            leghandles(legindex) = axhandles((legendloop*3));
            legindex = legindex + 1;
         end
           
         legend(leghandles, leg, 4);
           
         %calculate metric of gradient tuning selectivity  by first averaging together the mean disparities
         %then subtract off spontaneous from the peak and trough.  then divide the peak by the trough
         Gmean(:,1) = Gmean(:,1)/length(unique_mean_disp);
         Gmean(:,2) = Gmean(:,2)/length(unique_mean_disp);
           
         %calc average GDI
         [avgGDI(j), var_term] = compute_DDI(GDIdata(:,1)', GDIdata(:,2)');
           
         %readjust mean disparity responses to fall on the same mean
         %then calc avg GDI
         %shifted_graphs = figure;
         %set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [750 50 500 773], 'Name', 'Mean Adjusted Tilt Tuning Curves');
         total_mean = mean(GDIdata(:,2));
         for count_meandisp = 1:length(unique_mean_disp)
            disp_mean = mean(GDIdata(start(count_meandisp):stop(count_meandisp),2));
            difference = total_mean - disp_mean;
            GDIdata(start(count_meandisp):stop(count_meandisp),2) = GDIdata(start(count_meandisp):stop(count_meandisp),2) + difference;
            %figure(shifted_graphs);
            %hold on
            %PlotTuningCurve(GDIdata(start(count_meandisp):stop(count_meandisp),1), GDIdata(start(count_meandisp):stop(count_meandisp),2), symbols{count_meandisp}, lines{count_meandisp}, 1, 1);
         end
         %hold on
         [avgGDI_adj(j), var_term] = compute_DDI(GDIdata(:,1)', GDIdata(:,2)');
         [px, py, perr, spk_max, spk_min] = PlotTuningCurve(GDIdata(:,1), GDIdata(:,2), symbols{count_meandisp+1}, lines{count_meandisp+1}, 1, 0);
         
         print_adj = 0;
         if print_adj == 1
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
            FILEOUT = [FILE(1:filesize) 'adj_grad_curve'];
            fileid = [PATHOUT FILEOUT];
            printflag = 0;
            if (exist(fileid, 'file') == 0)    %file does not yet exist
               printflag = 1;
            end
            proffid = fopen(fileid, 'a');
            if (printflag)
               fprintf(proffid,'HDisp\tAvgResp\tStdErr\tAngle\tSpon\n');
               printflag = 0;
            end
            for go =1:length(px)
               if (go<=2)
                  fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t%6.2f\t%6.2f\n', px(go), py(go), perr(go), null_x(go),null_y(go));
               else
                  fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t\t\n', px(go), py(go), perr(go));
               end
            end
            
            for count_meandisp = 1:length(unique_mean_disp)
               [px, py, perr, spk_max, spk_min] = PlotTuningCurve(GDIdata(start(count_meandisp):stop(count_meandisp),1), GDIdata(start(count_meandisp):stop(count_meandisp),2), symbols{count_meandisp}, lines{count_meandisp}, 1, 0);
               for go =1:length(px)
                  fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t\t\n', px(go), py(go), perr(go));
               end
            end
            
            fclose(proffid);
         end
       
                  
         %calculate average VDI
         zero_data = min(GDIvergdata(:,2));
         GDIvergdata(:,2) = GDIvergdata(:,2) + abs(zero_data);
         [avgGDIverg(j), var_term] = compute_DDI(GDIvergdata(:,1)', GDIvergdata(:,2)');
         
         %readjust mean disparity responses to fall on the same mean
         %then calc avg VDI
         total_mean = mean(GDIvergdata(:,2));
         for count_meandisp = 1:length(unique_mean_disp)
            disp_mean = mean(GDIvergdata(start_verg(count_meandisp):stop_verg(count_meandisp),2));
            difference = total_mean - disp_mean;
            GDIvergdata(start_verg(count_meandisp):stop_verg(count_meandisp),2) = GDIvergdata(start_verg(count_meandisp):stop_verg(count_meandisp),2) + difference;
         end
         
         [avgGDIverg_adj(j), var_term] = compute_DDI(GDIvergdata(:,1)', GDIvergdata(:,2)');
           
         figure(graph);
         if(num_ap_size >= num_mag_disp)
            subplot(num_ap_size, num_mag_disp,  (k-1)*(num_mag_disp) + j);
         elseif(num_ap_size < num_mag_disp)
            subplot(num_mag_disp, num_ap_size,  (k-1)*(num_ap_size) + j);
         end
         hold on;
         [px, py, perr, spk_max, spk_min] = PlotTuningCurve(Gmean(:,1), Gmean(:,2), symbols{l+1}, lines{l+1}, 1, 1);
         mean_max = spk_max.x;
         mean_min = spk_min.x;
         
         Vmean(:,1) = Vmean(:,1)/length(unique_mean_disp);
         Vmean(:,2) = Vmean(:,2)/length(unique_mean_disp);
                                
         figure(ver_graph);
         if(num_ap_size >= num_mag_disp)
            subplot(num_ap_size, num_mag_disp,  (k-1)*(num_mag_disp) + j);
         elseif(num_ap_size < num_mag_disp)
            subplot(num_mag_disp, num_ap_size,  (k-1)*(num_ap_size) + j);
         end
         hold on;
         [px, py, perr, verg_max, verg_min] = PlotTuningCurve(Vmean(:,1), Vmean(:,2), symbols{l+1}, lines{l+1}, 0, 1);
         ax_han = gca;
         set(ax_han, 'box', 'on');
           
         hold on;
         figure(graph);
         null_x = [min(px) max(px)];
         null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
         null_y = [null_rate null_rate];
         hold on;
    	   plot(null_x, null_y, 'k--');
    	   yl = YLim;
         YLim([0 yl(2)])	% set the lower limit of the Y axis to zero
         hold off;
         hold on
           
         %now that we have the spontaneous rate, we need to subtract this from the mean peak and mean trough.
         %mean_peak = spk_max.y;
         %mean_trough = spk_min.y;
         GTI(j) = (spk_max.y-spk_min.y)/(spk_max.y-null_rate);
           
         hold off;
         figure(graph);
            
	      yl = YLim;
         YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
         XLabel('Tilt Angle (deg)');
         YLabel('Response (spikes/sec)');
            
         height = axis;
		   yheight = height(4);
		   string = sprintf('Ap Size = %4.1f', unique_ap_size(j));
		   text(height(1)+2, 0.95*yheight, string, 'FontSize', 8);
		   string = sprintf('Disp Mag = %4.2f', unique_mag_disp(k));
		   text(height(1)+2, 0.9*yheight, string, 'FontSize', 8);  
         string = sprintf('avgGDI = %2.4f', avgGDI(j));
         text(height(1)+2, 0.85*yheight, string, 'FontSize', 8);  
         string = sprintf('avgGDI mean adj = %2.4f, avgVDI mean adj = %2.4f', avgGDI_adj(j), avgGDIverg_adj(j));
         text(height(1)+2, 0.8*yheight, string, 'FontSize', 8);  
         string = sprintf('GTI = %2.4f', GTI(j));
		   text(height(1)+2, 0.75*yheight, string, 'FontSize', 8);  
         string = sprintf('Mean Peak = %2.4f', spk_max.x);
		   text(height(1)+2, 0.7*yheight, string, 'FontSize', 8);  
         string = sprintf('Mean Trough = %2.4f', spk_min.x);
		   text(height(1)+2, 0.65*yheight, string, 'FontSize', 8);  
           
         hold off;
         figure(ver_graph);
            
         XLabel('Tilt Angle (deg)');
         YLabel('Visual Angle (deg)');
            
         height = axis;
			yheight = height(4)-height(3);
			string = sprintf('Ap Size = %4.1f', unique_ap_size(j));
			text(height(1)+2, 0.95*yheight+height(3), string, 'FontSize', 8);
			string = sprintf('Disp Mag = %4.2f', unique_mag_disp(k));
			text(height(1)+2, 0.85*yheight+height(3), string, 'FontSize', 8);  
              
         line = sprintf('\t%3.1d\t%3.5f\t%1.12f\t%1.12f\t%3.2f\t%3.2f\t%6.1f\t%2.4f\t%2.4f\t%3.2f\t%3.2f\t%3.2f', unique_ap_size(j), eccentricity, avgGDI_adj(j), avgGDIverg_adj(j), spk_max.x, verg_max.x, data.neuron_params(PREFERRED_DIRECTION, 1), x_ctr, y_ctr, mean_max, mean_min, null_y(1));
         data_string = strcat(data_string, line);

         if j == length(unique_ap_size)
            for temp_counter = 1:length(unique_mean_disp)
                 if p_val(j, temp_counter) < .01
                     IndTDI_vs_ALLTDI_temp = sprintf('%s\t%1.12f\t%1.12f\n', FILE, avgGDI_adj(j), GDI(temp_counter, 1));
                     IndTDI_vs_ALLTDI = sprintf('%s%s', IndTDI_vs_ALLTDI, IndTDI_vs_ALLTDI_temp)
                 end
            end
         end
      end %end mag disp
      hold off;
   end % end ap_size
    
   printme = 0;
   if (printme==1)
      %pathsize = size(PATH,2) - 1;
      %while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
      %    pathsize = pathsize - 1;
      %end   
      PATHOUT = 'Z:\Users\Jerry\GradAnalysis\figure_data\';
            
      line = sprintf('%s\t%d\t', FILE, SpikeChan);
      data_string = strcat(line, data_string);
      angle_out = strcat(line, angle_out);
      GDI_mdisp_out = strcat(line, GDI_mdisp_out);
      
      %print grad metrics
      outfile = [PATHOUT 'SU_Grad_metrics.dat'];
      printflag = 0;
      if (exist(outfile, 'file') == 0)    %file does not yet exist
         printflag = 1;
      end
      fid = fopen(outfile, 'a');
      if (printflag)
         fprintf(fid,'File\tChannel\tApSize1\tEcc1\tTDIadj1\tTDIvergadg1\tmax1\tvergmax1\tPrefDir1\txctr1\tyctr1\tmeanmax1\tmeanmin1\tSpont1\tApSize2\tEcc2\tTDIadj2\tTDIvergadg2\tmax2\tvergmax2\tPrefDir2\txctr2\tyctr2\tmeanmax2\tmeanmin2\tSpont2\n')
         printflag = 0;
      end

      if (length(unique_ap_size) == 2)
         fprintf(fid, '%s', [data_string]);
      else
         fprintf(fid, '%s\t\t\t\t\t\t\t\t\t\t\t', [data_string]);
      end
      fprintf(fid, '\r\n');
      fclose(fid);
 
      %print mean disp angles
      outfile = [PATHOUT 'SU_Mean_Disp_angles.dat'];
      fid = fopen(outfile, 'a');
      fprintf(fid, '%s', [angle_out]);
      fprintf(fid, '\r\n');
      fclose(fid);        
      
      %print TDI for each mean disp curve
      outfile = [PATHOUT 'SU_Mean_Disp_TDI.dat'];
      fid = fopen(outfile, 'a');
      fprintf(fid, '%s', [GDI_mdisp_out]);
      fprintf(fid, '\r\n');
      fclose(fid);        
      
      %print out overall TDI vs Individual TDIs
      %outfile = [PATHOUT 'AllTDI_vs_IndTDI.dat'];
      %printflag = 0;
      %if (exist(outfile, 'file') == 0)    %file does not yet exist
      %   printflag = 1;
      %end
      %fid = fopen(outfile, 'a');
      %if (printflag)
      %   fprintf(fid,'File\tAllTDI\tIndTDI\n');
      %   printflag = 0;
      %end

      %fprintf(fid, '%s', [IndTDI_vs_ALLTDI]);

      %fclose(fid);
    end
    
   ancova_var = HGradAnova(graph, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
   ver_ancova_var = HGradAnovaVerg(ver_graph, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    
    
   for i  = 1:length(unique_ap_size)
      [h,atab,ctab,stats] = aoctool(ver_ancova_var{i}(:,2), ancova_var{i}(:,2), ancova_var{i}(:,1), .05, 'Vergence', 'Response', 'Angle', 'off')
    
      figure(gen_data_fig);
      %print out the values associated with each parameter graphed    
      subplot(2, 1, 2);
      %get out the column name for the ancova
      s = atab(1,6);
      se = atab(2:4, 6);
      for j=1:length(se)
         f_num(i,j) = se{j};
      end
   end
   
   if i>1
      variable_names = {'Ap Size:' 'Disp Mag:' 'M.Disp:' 'Disp Ang:' 'Avg GDI' char(s) char(s) 'Speed' 'Coh'};
      variable_data = {unique_ap_size unique_mag_disp unique_mean_disp unique_disp_ang avgGDI f_num(1,:)' f_num(2,:)' unique_speed unique_coh};
      PrintDataArray(variable_data, variable_names);
   else
      variable_names = {'Ap Size:' 'Disp Mag:' 'M.Disp:' 'Disp Ang:' 'Avg GDI' char(s) 'Speed' 'Coh'};
      variable_data = {unique_ap_size unique_mag_disp unique_mean_disp unique_disp_ang avgGDI f_num(1,:)' unique_speed unique_coh};
      PrintDataArray(variable_data, variable_names);
   end

return;


