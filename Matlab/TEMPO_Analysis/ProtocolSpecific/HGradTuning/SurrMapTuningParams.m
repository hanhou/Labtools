%-----------------------------------------------------------------------------------------------------------------------
%-- Surround Mapping curve will plot a polar plot for each set of surrounds -JDN 3/28/00
%-----------------------------------------------------------------------------------------------------------------------
function [ctr_only, ctrdisp, surrsize, ctrsize, unique_ang, px, py, plot_x, plot_y] = SurrMapTuningParams(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, UseSyncPulses);

	TEMPO_Defs;
   
   %following function checks offset times and outputs starting and ending analysis (in spike bins)
   num_trials = size(data.event_data, 3);
   [StartOffsetBin StopOffsetBin] = CheckTimeOffset(data, num_trials, StartCode, StopCode, StartOffset, StopOffset, UseSyncPulses);

   if (~isempty(data.spike_data))
      %compute the firing rate over all trials during the period between StartCode and StopCode
      data.spike_rates = ComputeSpikeRates(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin);
   end
   
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
   
   %this is alright since we only use one size to map
   surrsize = unique_ap_size(2);
   
   ctrsize = data.dots_params(DOTS_AP_XSIZ,BegTrial:EndTrial,PATCH4);
   unique_ctr_size = munique(ctrsize(~null_trials)');
   ctrsize = unique_ctr_size;
   
   ctrdisp = data.dots_params(DOTS_HDISP, BegTrial:EndTrial, PATCH4);
   unique_ctr_disp = munique(ctrdisp(~null_trials)');
   ctrdisp = unique_ctr_disp;
   
   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);
   
   %get indices of monoc. and uncorrelated controls
   control_trials = logical( (disp == LEYE_CONTROL) | (disp == REYE_CONTROL) | (disp == UNCORR_CONTROL) );
   
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:length(disp);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
   %now plot a horizontal disparity tuning curve for each point
   unique_ang = unique_ang(2:length(unique_ang));
   
   ctr_resp = spike_rates(~null_trials & patch_off_condition & ~control_trials & select_trials);
   ctr_only = mean(ctr_resp);

   for i = 1:length(unique_ang)
      ang_select = logical(ang == unique_ang(i));

      plot_x(i, :) = disp(ang_select & ~null_trials & ~patch_off_condition & ~control_trials & select_trials);
      plot_y(i, :) = spike_rates(ang_select & ~null_trials & ~patch_off_condition & ~control_trials & select_trials);
      
      %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
      [px(:, i), py(:, i), perr(:, i)] = PlotTuningCurve(plot_x(i, :)', plot_y(i, :)', symbols{i}, lines{i}, 1, 0);
      null_x = [min(px) max(px)];
      null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
      null_y = [null_rate null_rate];
   end
   
return;
