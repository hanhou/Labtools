%-----------------------------------------------------------------------------------------------------------------------
%-- HDispTuningCurve.m -- Plots a horizontal disparity tuning curve, possibly for multiple sizes, and including
%--	monoc and uncorrelated control conditions.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function HDispTuningCurveSize(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	TEMPO_Defs;
   
   symbols = {'r*' 'bo' 'go' 'ko' 'b*' 'r*' 'g*' 'k*'};
   lines = {'r-' 'b-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--'};
   
	%get the column of values of horiz. disparity in the dots_params matrix
	hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
   
   %get indices of any NULL conditions (for measuring spontaneous activity)
   null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );
   
   %get the column of size values
   size = data.dots_params(DOTS_AP_XSIZ,:,PATCH1);
   unique_size = munique(size(~null_trials)');
   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);

   %get indices of monoc. and uncorrelated controls
   control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );
   
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:length(hor_disp);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
   figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 250 500 573], 'Name', 'Horizontal Disparity Tuning Curve');
   subplot(2, 1, 2);
   
   for i=1:length(unique_size)	%for each different size value, plot a separate disparity tuning curve
	   size_select = logical( (size == unique_size(i)) );
        
      plot_x = hor_disp(size_select & ~null_trials & ~control_trials & select_trials);
   	plot_y = spike_rates(size_select & ~null_trials & ~control_trials & select_trials); 
      
   	%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
	   hold on;
      [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
      hold off;
      
      hold on;
      errorbar(px, py, perr, perr, symbols{i});
   	hold off;
      
      p_value(i) = spk_anova(plot_y, plot_x, px);
      avg_resp(i) = mean(plot_y);      
      
   end
   
   %now, get the firing rate for NULL condition trials and add spontaneous rate to plot
   null_x = [min(px) max(px)];
   null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
   null_y = [null_rate null_rate];
   hold on;
   plot(null_x, null_y, 'k--');
   hold off;
   
   yl = YLim;
   YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
   XLabel('Horizontal Disparity(deg)');
   YLabel('Response (spikes/sec)');
   
   %now, print out some useful information in the upper subplot
   subplot(2, 1, 1);
   PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
 
return;