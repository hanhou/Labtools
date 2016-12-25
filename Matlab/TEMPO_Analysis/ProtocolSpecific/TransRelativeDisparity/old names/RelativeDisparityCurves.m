%-----------------------------------------------------------------------------------------------------------------------
% RelativeDisparityCurves.m -- Module to display center disparity tuning for different surround
%	disparities, or vice-versa.  Starting 6/20/00...
%-----------------------------------------------------------------------------------------------------------------------

function TransRelativeDisparityCurves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	TEMPO_Defs;

	symbols = {'ko' 'r*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
   lines = {'k-' 'r--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

	%get the column of values of horiz. disparity of center in the dots_params matrix
	hor_disp_ctr = data.dots_params(DOTS_HDISP,:,PATCH1);
   
	%get the column of values of horiz. disparity of surround in the dots_params matrix
   hor_disp_surr = data.dots_params(DOTS_HDISP,:,PATCH4);
   
   %get indices of any NULL conditions (for measuring spontaneous activity)
   null_trials = logical( (hor_disp_surr == data.one_time_params(NULL_VALUE)) );
   
   % Surr_Control = no dots in surround, ctr disp varies
   surr_control = logical(hor_disp_surr == data.one_time_params(PATCH_OFF));
   
   % Ctr_Control = no dots in ctr, surround varies
   ctr_control = logical(hor_disp_ctr == data.one_time_params(PATCH_OFF));

   
   unique_disp_surr = munique(hor_disp_surr(~null_trials)');
   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);
   
	%get indices of monoc. and uncorrelated controls
   control_trials = logical( (hor_disp_surr == LEYE_CONTROL) | (hor_disp_surr == REYE_CONTROL) | (hor_disp_surr == UNCORR_CONTROL) );
   
	%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:length(hor_disp_ctr);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );


	figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 150 500 473], 'Name', 'Relative Disparity Tuning Curve');
   subplot(2, 1, 2);
   
   all_resp_data = [];
   
   for i=1:length(unique_disp_surr)	%for each different surround disparity value, plot a separate disparity tuning curve
      surr_disp_select = logical( (hor_disp_surr == unique_disp_surr(i)) );
      
		plot_x = hor_disp_ctr(surr_disp_select & ~null_trials & ~ctr_control & ~control_trials & select_trials);
      plot_y = spike_rates(surr_disp_select & ~null_trials & ~ctr_control & ~control_trials & select_trials); 
      
    	%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
	   hold on;
      [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
      
      all_resp_data{i}.cdisp = px;
      all_resp_data{i}.resp = py;
      all_resp_data{i}.resp_err = perr;
      
      if unique_disp_surr(i) ~= data.one_time_params(PATCH_OFF)
      	surr_control_trials = logical(ctr_control & surr_disp_select);
         if ( sum(surr_control_trials) > 0 )
	         surr_resp = spike_rates(surr_control_trials);
   	   	hold on;
      		errorbar(max(px)*1.07, mean(surr_resp), std(surr_resp)/sqrt(sum(surr_control)), std(surr_resp)/sqrt(sum(surr_control)), symbols{i});
            text(max(px)*1.12, mean(surr_resp), num2str(unique_disp_surr(i)));
         end
      end
         
      p_value(i) = spk_anova(plot_y, plot_x, px);
      avg_resp(i) = mean(plot_y); 
      hold on;
      H(i) = plot(px, py, symbols{i});
	   legend_string{i} = sprintf('%5.2f', unique_disp_surr(i));
	end
   
   %write out data in a form for 2-way ANOVA
   unique_disp_ctr = munique(hor_disp_ctr(~null_trials & ~ctr_control)');
   for i=1:length(unique_disp_surr)	%for each different surround disparity value, plot a separate disparity tuning curve
      for j=1:length(unique_disp_ctr)
         select1 = ( logical(hor_disp_surr == unique_disp_surr(i)) & ~surr_control);
         select2 = ( logical(hor_disp_ctr == unique_disp_ctr(j)) & ~ctr_control);
         %[hor_disp_ctr(select1&select2)' hor_disp_surr(select1&select2)' spike_rates(select1&select2)']
         temp = [];
         temp = spike_rates(select1&select2);
         for k=1:length(temp)
            buff = sprintf('%d %d %6.3f', j, i, temp(k));
            %disp(buff);
         end
      end
   end   
   
   
   %write out data in a form for plotting with Origin, etc.
   for j=1:length(unique_disp_ctr)
      line = '';
      for i=1:length(unique_disp_surr)	
         buff = sprintf('%7.3f %7.2f %7.3f ', all_resp_data{i}.cdisp(j), all_resp_data{i}.resp(j), all_resp_data{i}.resp_err(j) );
         line = [line buff];
      end
      %disp(line);
   end
   
   yl = YLim;
   YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
   XLabel('Horizontal Disparity of Center(deg)');
   YLabel('Response (spikes/sec)');
   legend(H, legend_string, -1);

   %now, get the firing rate for NULL condition trials and add spontaneous rate to plot
   null_x = [min(px) max(px)];
   null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
   null_y = [null_rate null_rate];
   hold on;
   plot(null_x, null_y, 'k--');
   hold off;
   
     %now, print out some useful information in the upper subplot
   subplot(2, 1, 1);
   PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
 
   %now print out useful values for H.Disp specific 
   % pmax, pmin, py, 
   PrintRelDispData(data, p_value, avg_resp, pmax, pmin, px, null_rate, unique_disp_surr, PATH, FILE);
         
   return;
   
   