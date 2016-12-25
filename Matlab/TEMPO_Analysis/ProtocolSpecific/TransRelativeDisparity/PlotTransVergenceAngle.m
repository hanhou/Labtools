%Vergence Angle during Relative Disparity Fixations = PlotVergenceAngle
function PlotTransVergenceAngle(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE)

TEMPO_Defs;

symbols = {'ko' 'r*' 'go' 'mo' 'b*' 'c*' 'ro' 'g*' };
lines = {'k-' 'r--' 'g-' 'm-' 'b--' 'c--' 'r-' 'g--' };

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
   
%get the average horizontal eye positions to calculate vergence
Leyex_positions = data.eye_positions(1, :);
Reyex_positions = data.eye_positions(3, :);
   
vergence = Leyex_positions - Reyex_positions;


%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp_surr == LEYE_CONTROL) | (hor_disp_surr == REYE_CONTROL) | (hor_disp_surr == UNCORR_CONTROL) );
   
%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp_ctr);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%p_value_fig = figure;
   
figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 250 500 573], 'Name', 'Vergence Angle');
   subplot(2, 1, 2);
   
   all_verg_data = [];
   
   for i=1:length(unique_disp_surr)	%for each different surround disparity value, plot a separate disparity tuning curve
	   surr_disp_select = logical( (hor_disp_surr == unique_disp_surr(i)) );

		plot_x = hor_disp_ctr(surr_disp_select & ~null_trials & ~ctr_control & ~control_trials & select_trials);
      plot_y = vergence(surr_disp_select & ~null_trials & ~ctr_control & ~control_trials & select_trials); 
      
    	%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
	   hold on;
      [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
      
      all_verg_data{i}.cdisp = px;
      all_verg_data{i}.verg = py;
      all_verg_data{i}.verg_err = perr;
            
      if unique_disp_surr(i) ~= data.one_time_params(PATCH_OFF)
      	surr_control_trials = logical(ctr_control & surr_disp_select);
         if ( sum(surr_control_trials) > 0 )
	     		surr_resp = vergence(surr_control_trials);
   	   	hold on;
      		errorbar(max(px)*1.07, mean(surr_resp), std(surr_resp)/sqrt(sum(surr_control)), std(surr_resp)/sqrt(sum(surr_control)), symbols{i});
            text(max(px)*1.12, mean(surr_resp), num2str(unique_disp_surr(i)));
         end
      end
     
    p_value(i) = spk_anova(plot_y, plot_x, px);
    avg_resp(i) = mean(plot_y);      
  end


   %write out data in a form for 2-way ANOVA
   unique_disp_ctr = munique(hor_disp_ctr(~null_trials & ~ctr_control)');
   for i=1:length(unique_disp_surr)	%for each different surround disparity value, plot a separate disparity tuning curve
      for j=1:length(unique_disp_ctr)
         select1 = ( logical(hor_disp_surr == unique_disp_surr(i)) & ~surr_control);
         select2 = ( logical(hor_disp_ctr == unique_disp_ctr(j)) & ~ctr_control);
         %[hor_disp_ctr(select1&select2)' hor_disp_surr(select1&select2)' vergence(select1&select2)']
         temp = [];
         temp = vergence(select1&select2);
         for k=1:length(temp)
            buff = sprintf('%d %d %6.3f', j, i, temp(k));
            disp(buff);
         end
      end
   end   

   %write out data in a form for plotting with Origin, etc.
   for j=1:length(unique_disp_ctr)
      line = '';
      for i=1:length(unique_disp_surr)	
         buff = sprintf('%7.3f %7.5f %7.5f ', all_verg_data{i}.cdisp(j), all_verg_data{i}.verg(j), all_verg_data{i}.verg_err(j) );
         line = [line buff];
      end
      disp(line);
   end


	XLabel('Horizontal Disparity of Center(deg)');
   YLabel('Vergence Angle');
   

    %now, print out some useful information in the upper subplot
   subplot(2, 1, 1);
   PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   
   %now print out useful values for Vergence Angle 
   % pmax, pmin, py, 
   %


PrintTransVergenceData(p_value, pmax, pmin, px, unique_disp_surr, PATH, FILE);

return