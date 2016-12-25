%-----------------------------------------------------------------------------------------------------------------------
%-- ContrastResponseCurve.m -- Plots contrast-response curves.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function ContrastResponseCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
   
   symbols = {'r*' 'bo' 'go' 'ko' 'b*' 'r*' 'g*' 'k*'};
   lines = {'r-' 'b-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--'};
   
	%get the column of values of speeds in the dots_params matrix
	speed = data.dots_params(DOTS_SPEED,:,PATCH1);
   
  	%get the column of values of directions in the dots_params matrix
	direction = data.dots_params(DOTS_DIREC,:,PATCH1);
   
	%get the column of values of speeds in the dots_params matrix
	contrast = data.dots_params(DOTS_CONTRAST,:,PATCH1);
   
   %get the contrast values for PATCH2
	p2_status = data.dots_params(DOTS_AP_STATUS,:,PATCH2);
   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);
   
   %[contrast' direction' speed' spike_rates']
   
   %get indices of any NULL conditions (for measuring spontaneous activity
   null_trials = logical( (contrast == data.one_time_params(NULL_VALUE)) );
   
   % get unique values of vars
   unique_speed = munique(speed(~null_trials)');
   unique_direc = munique(direction(~null_trials)');
   unique_contrast = munique(contrast(~null_trials)');
   
   %now, remove trials from speed and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:length(contrast);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
%   [speed' null_trials' select_trials']

   figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 250 500 573], 'Name', 'Contrast Response Curve');
   subplot(2, 1, 2);

	if (length(unique_direc) > 1)	
      unique_vals = unique_direc;
      values = direction;
	elseif (length(unique_speed) >= 1)	
      unique_vals = unique_speed;
      values = speed;
   end
   
   for i=1:length(unique_vals)
	   value_select = logical( (values == unique_vals(i)) );
      
	   plot_x = contrast(~null_trials & select_trials & value_select);
      plot_y = spike_rates(~null_trials & select_trials & value_select);
   
	   %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
      hold on;
      [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 0, 1);    
	   %errorbar(px, py, perr, perr, 'bo');
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
   XLabel('Contrast');
   YLabel('Response (spikes/sec)');
   
   %now, print out some useful information in the upper subplot
   subplot(2, 1, 1);
   PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   
   % first calculate some metrics and stats
   p_value = spk_anova(plot_y, plot_x, px)
   avg_resp = mean(plot_y);
   
   %write raw data out to a file
   i = size(PATH,2) - 1;
	while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
   	i = i - 1;
   end   
   PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
   i = size(FILE,2) - 1;
   while FILE(i) ~='.'
     	i = i - 1;
   end
   FILEOUT = [FILE(1:i) 'ct'];
     
   fileid = [PATHOUT FILEOUT];
   fwriteid = eval(['fopen(fileid, ''w'')']);
   % header for data fields
   fprintf(fwriteid, 'Contrast	Speed		Direction	Response  PATCH2OnOFF\n');
   
   for j = 1:length(contrast)
     	fprintf(fwriteid, '%8.4f   %6.2f    %8.1f       %7.4f      %8.4f\n', contrast(j), speed(j), direction(j), spike_rates(j), p2_status(j));	
   end    
   fclose(fwriteid);
   
return;