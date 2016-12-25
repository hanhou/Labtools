%-----------------------------------------------------------------------------------------------------------------------
%-- NeuroPsychoPlot.m -- Plots neurometric and psychometric functions on the same plot; used ROC analysis for neurometrics
%--	GCD, 5/26/00
%-----------------------------------------------------------------------------------------------------------------------
function NeuroPsychoPlot_bnd(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	TEMPO_Defs;		
   
   [conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

    
   
   symbols = {'bo' 'r*' 'gd' 'kx' 'bs' 'rv' 'g*' 'k*' 'c*'};
   lines = {'b-' 'r-.' 'g:' 'k--' 'b--' 'r--' 'g--' 'k--' 'c--'};
   
   var_values = conditions(6,:);
   unique_var_values = munique(var_values');
   
   % use more intelligent method for determining parameters varied across trials
   %now, select trials that fall between BegTrial and EndTrial
   trials = 1:size(conditions,2);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

   figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Neurometric/Psychometric Comparison');
   subplot(2, 1, 2);
   
   dx = 0.1;
   
   %% ********* NEUROMETRIC ANALYSIS ********************
   %loop through each binocular correlation levels, and do ROC analysis for each
 %  fit_x = min(unique_var_values) :dx: max(unique_var_values);
 fit_x = 0 :dx: 360;

      
   %% *********** PSYCHOMETRIC ANALYSIS ****************************
   pct_correct = []; N_obs = [];
   for i=1:length(unique_var_values)
        trials = (var_values == unique_var_values(i) ) & select_trials;
   	    correct_trials = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
        pct_correct(i) = sum(correct_trials)/sum(trials);
      	N_obs(i) = sum(trials);
      	% data for Weibull fit
   	    fit_data(i,1) = unique_var_values(i);
        fit_data(i,2) = pct_correct(i);
        fit_data(i,3) = N_obs(i);
   end
   hold on;
   Handl(2) = plot(unique_var_values, pct_correct, 'ko');
   hold off;
   
   [monkey_alpha monkey_beta] = weibull_fit_bnd(fit_data);
   monkey_fit_y = weibull_curve_bnd(fit_x, [monkey_alpha monkey_beta]);
   hold on;
   plot(fit_x, monkey_fit_y, 'k--');
   hold off;
   xlabel('Obj Phase Jitter (degrees)');
   ylabel('Fraction Correct');
   
   YLim([0 1]);
   %comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
   % set(gca, 'XScale', 'log');
   XLim([0 360]);
   
   %now, print out some useful information in the upper subplot
   subplot(2, 1, 1);
   PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   
   %now, print out some specific useful info.
   xpos = -10; ypos = 10;
   font_size = 11;
   bump_size = 8;
    line = sprintf( 'Behavioral Threshold = %6.3f, slope = %6.3f\n', monkey_alpha, monkey_beta );
   	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   
return;