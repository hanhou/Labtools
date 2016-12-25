%-----------------------------------------------------------------------------------------------------------------------
%-- NeuroPsychoPlot.m -- Plots neurometric and psychometric functions on the same plot; used ROC analysis for neurometrics
%--	GCD, 5/26/00
%-----------------------------------------------------------------------------------------------------------------------
function NeuroPsychoPlot_bnd(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	TEMPO_Defs;		
   
   [conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

    
   BKGND = 0;
   OBJ = 1;
   BAR = 2;
   
   %currently only handles obj 1
   obj_num = 1;
   active_bars = data.bar_params(BAR_STATUS, 1, obj_num, :);
   active_bars = squeeze(active_bars);
   bar_num = active_bars(1);
   
   symbols = {'bo' 'r*' 'gd' 'kx' 'bs' 'rv' 'g*' 'k*' 'c*'};
   lines = {'b-' 'r-.' 'g:' 'k--' 'b--' 'r--' 'g--' 'k--' 'c--'};
   
   %% use object, background, and bar keywords
   ParamList = {'Bkgnd Color ' 'Obj Y Position' 'Bkgnd H disp '};
   
   % use more intelligent method for determining parameters varied across trials
   
   var_num = 1;
	%need to make analysis backwards compatible
   if (size(data.obj_params,1) < NUM_OBJ_PARAMS)
      num_obj_vars = size(data.obj_params,1);      
   else 
      num_obj_vars = NUM_OBJ_PARAMS;
   end
   
   if (size(data.bkgnd_params,1) < NUM_BKGND_PARAMS)
      num_bkgnd_vars = size(data.bkgnd_params,1);      
   else 
      num_bkgnd_vars = NUM_BKGND_PARAMS;
   end

   
   % skip obj status and obj_dir for now 

   for param = OBJ_SKEW_AMPL: num_obj_vars      
  		unique_vals = munique(data.obj_params(param ,:, obj_num)');
        if (length(unique_vals) > 1  )
        		if (param == OBJ_JITTER)%% || (param == OBJ_LUM_MULT) 
            	keyed_param = var_num;     
            end 
            
         	var_index(var_num) = param;
         	var_values(var_num,:) = data.obj_params(param ,:,obj_num);
         	var_type(var_num) = OBJ;
            unique_var_values{var_num} = unique_vals;
         	var_num = var_num + 1;
			end
  	end
   
   for param = 2: num_bkgnd_vars      
         unique_vals = munique(data.bkgnd_params(param ,:)'); 
         if  (length(unique_vals) > 1 )
            var_index(var_num) = param;
         	var_values(var_num,:) = data.bkgnd_params(param ,:);
            var_type(var_num) = BKGND;
            unique_var_values{var_num} = unique_vals;
            var_num = var_num + 1;
         end  
   end
   
   % decrement counter by 1 to find final number of nested variables used in experiment.
	var_num = var_num - 1;
   
   % & ~(param == BKGND_DOT_COLOR) & (length(unique_vals) == length(var_values(var_num - 1, :) ) ) 

	if (keyed_param == 1)
   	anal_order = [2 1];
	else 
		anal_order = [1 2];
   end   
      
   %now, select trials that fall between BegTrial and EndTrial
   trials = 1:length(var_values(1,:));		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

   figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 250 500 573], 'Name', 'Neurometric/Psychometric Comparison');
   subplot(2, 1, 2);
   
   dx = 0.1;
   
   %% ********* NEUROMETRIC ANALYSIS ********************
   %loop through each binocular correlation levels, and do ROC analysis for each
   fit_x = unique_var_values{keyed_param}(1):dx: unique_var_values{keyed_param}(end);

      
   %% *********** PSYCHOMETRIC ANALYSIS ****************************
   pct_correct = []; N_obs = [];
   for i=1:length(unique_var_values{anal_order(1)})
      for j = 1:length(unique_var_values{anal_order(2)})
	      trials = (var_values(anal_order(1),:) == unique_var_values{anal_order(1)}(i) ) & (var_values(anal_order(2),:) == unique_var_values{anal_order(2)}(j) ) & select_trials;
   	   correct_trials = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
         pct_correct(i,j) = sum(correct_trials)/sum(trials);
         
      	N_obs(j) = sum(trials);
      	% data for Weibull fit
   	   fit_data(j,1) = unique_var_values{anal_order(2)}(j);
    	   fit_data(j,2) = pct_correct(i,j);
    	   %fit_data(j,2) = pct_incorrect(j);
         fit_data(j,3) = N_obs(j);
     	end
      hold on;
      Handl(2) = plot(unique_var_values{anal_order(2)}, pct_correct(i,:), symbols{i});
      [monkey_alpha(i) monkey_beta(i)] = inv_weibull_fit(fit_data);
   	monkey_fit_y(i,:) = inv_weibull_curve(fit_x, [monkey_alpha(i) monkey_beta(i)]);
      hold on;
      plot(fit_x', monkey_fit_y(i,:), lines{i} );
      hold off; 
		leg{2*i-1} = [ParamList{keyed_param} num2str(unique_var_values{anal_order(1)}(i)) ''];      
      leg{2*i} = [ParamList{keyed_param} num2str(unique_var_values{anal_order(1)}(i)) ''];  
   end
   
   pct_correct
   unique_var_values{1}
   unique_var_values{2}
   
   xlabel('Obj Phase Jitter (degrees)');
   ylabel('Fraction Correct');
   legend(leg);
   
   YLim([0.4 1]);
   %comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
   % set(gca, 'XScale', 'log');
   % XLim([1 100]);
   
   %now, print out some useful information in the upper subplot
   subplot(2, 1, 1);
   PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   
   %now, print out some specific useful info.
   xpos = -10; ypos = 10;
   font_size = 11;
   bump_size = 8;
   for i=1:length(unique_var_values{anal_order(1)})
	   line = sprintf( '%s %5.3f: Monkey: threshold = %6.3f, slope = %6.3f\n', ParamList{keyed_param}, unique_var_values{anal_order(1)}(i), monkey_alpha(i), monkey_beta(i) );
   	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   end
   
   %print some text to console, for batch purposes
%   buff = sprintf('%s  %4.1f %5.2f %5.3f %4.2f %4.2f %5.2f %6.3f %7.4f %6.3f %7.4f %5.3f %5.3f %4d', FILE, data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED),...
%      data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER), ...
%      neuron_alpha, neuron_beta, monkey_alpha, monkey_beta, unique_hdisp(1), unique_hdisp(2), (1+ EndTrial - BegTrial) );
%   disp(buff);
   
return;