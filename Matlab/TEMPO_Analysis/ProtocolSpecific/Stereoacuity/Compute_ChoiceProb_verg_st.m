%-----------------------------------------------------------------------------------------------------------------------
%-- Compute_ChoiceProb.m -- Uses ROC analysis to compute a choice probability for each different stimulus level
%--	GCD, 5/30/00
%-----------------------------------------------------------------------------------------------------------------------
function Compute_ChoiceProb_verg_st(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	TEMPO_Defs;		%defns like IN_T1_WIN_CD
	ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
   
   %disp('computing choice probabilities...');
   
   Pref_HDisp = data.one_time_params(PREFERRED_HDISP);
   
	%get the column of values of horiz. disparities in the dots_params matrix
   h_disp_p1 = data.dots_params(DOTS_HDISP,:,PATCH1);
   unique_hdisp_p1 = munique(h_disp_p1');
   
   %get the column of values of horiz. disparities in the dots_params matrix
   h_disp_p4 = data.dots_params(DOTS_HDISP,:,PATCH4);
   unique_hdisp_p4 = munique(h_disp_p4');
   
   %compute the unsigned horizontal disparity
   unsigned_hdisp = abs(h_disp_p1-h_disp_p4);
   unsigned_hdisp = (round(10000 * unsigned_hdisp)) / 10000;
   unique_unsigned_hdisp = munique(unsigned_hdisp');
   
   %get the average eye positions to calculate vergence
   Leyex_positions = data.eye_positions(1, :);
   Leyey_positions = data.eye_positions(2, :);
   Reyex_positions = data.eye_positions(3, :);
   Reyey_positions = data.eye_positions(4, :);
   
   vergence_h = Leyex_positions - Reyex_positions;
   vergence_v = Leyey_positions - Reyey_positions;
   
   if (data.eye_calib_done == 1)
        Leyex_positions = data.eye_positions_calibrated(1, :);
        Leyey_positions = data.eye_positions_calibrated(2, :);
        Reyex_positions = data.eye_positions_calibrated(3, :);
        Reyey_positions = data.eye_positions_calibrated(4, :);
         
        vergence_h = Leyex_positions - Reyex_positions;
        vergence_v = Leyey_positions - Reyey_positions;
   end
    
   %get the firing rates for all the trials 
    spike_rates = data.spike_rates(SpikeChan, :);

    %get indices of any NULL conditions (for measuring spontaneous activity
   null_trials = logical( (h_disp_p1 == data.one_time_params(NULL_VALUE)) );
   
   %now, select trials that fall between BegTrial and EndTrial
   trials = 1:length(h_disp_p1);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
   %now, determine the choice that was made for each trial, PREFERRED or NULL
   %by definition, a preferred choice will be made to Target1 and a null choice to Target 2
   %thus, look for the events IN_T1_WIN_CD and IN_T2_WIN_CD.  GCD, 5/30/2000
   num_trials = length(h_disp_p1);
   PREFERRED = 1;
   NULL = 2;
   for i=1:num_trials
      temp = data.event_data(1,:,i);
      events = temp(temp>0);  % all non-zero entries
      if (sum(events == IN_T1_WIN_CD) > 0)
         choice(i) = PREFERRED;
      elseif (sum(events == IN_T2_WIN_CD) > 0)
         choice(i) = NULL;
      else
         disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
      end        
   end
   
   %now, plot the spike distributions, sorted by choice, for each disparity level
   figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [50 150 500 473], 'Name', 'Choice Probabilities');
   num_disp = length(unique_hdisp_p1);
   num_unsigned_disp = length(unique_unsigned_hdisp);
   choice_prob = [];
   
   for i=1:num_disp
      subplot(num_unsigned_disp, 2, i);
      if (Pref_HDisp - unique_hdisp_p4(1) < 0)
	      near_choices = ( (choice == PREFERRED) & (h_disp_p1 == unique_hdisp_p1(i)) );
         near_dist_h{i} = vergence_h(near_choices & select_trials);
	      far_choices = ( (choice == NULL) & (h_disp_p1 == unique_hdisp_p1(i)) );
         far_dist_h{i} = vergence_h(far_choices & select_trials);
      elseif (Pref_HDisp - unique_hdisp_p4(1) > 0)
         far_choices = ( (choice == PREFERRED) & (h_disp_p1 == unique_hdisp_p1(i)) );
         far_dist_h{i} = vergence_h(far_choices & select_trials);
	      near_choices = ( (choice == NULL) & (h_disp_p1 == unique_hdisp_p1(i)) );
         near_dist_h{i} = vergence_h(near_choices & select_trials);
      end
      
      %plot the distributions.  This uses a function (in CommonTools) that I wrote.  GCD
      PlotTwoHists(far_dist_h{i}, near_dist_h{i});
         
	   lbl = sprintf('%5.3f deg', unique_hdisp_p1(i)-unique_hdisp_p4(1) );
      ylabel(lbl);
         
      if ( (length(far_dist_h{i}) > 0) & (length(near_dist_h{i}) > 0) )
       	choice_prob(i) = rocN(far_dist_h{i}, near_dist_h{i}, 100);
   	   cp = sprintf('%5.2f', choice_prob(i));
         xl = XLim; yl = YLim;
         text(xl(2), yl(2)/2, cp);
      end
   end    
   
   str = sprintf('%s vergence', FILE );
   xlabel(str);
      
   %far_dist{i} and near_dist{i} are cell arrays that hold the far and near choice
   %distributions for each disparity level.
   %NOW, we want to Z-score the distributions (far and near choices together) and combine across
   %disparities.  GCD, 8/10/00
   for i=1:num_disp
         %for each condition, combine the preferred and null choices into one dist., then find mean and std
         all_choices = []; mean_val = []; std_val = [];
         all_choices = [far_dist_h{i}  near_dist_h{i}];
         mean_val = mean(all_choices);
         std_val = std(all_choices);
         %now use the mean_val and std_val to Z-score the original distributions and store separately
         Z_far_dist{i} = (far_dist_h{i} - mean_val)/std_val;
         Z_near_dist{i} = (near_dist_h{i} - mean_val)/std_val;      
   end      
   
   %Now, combine data across correlation to get a grand choice probability, and plot distributions again
   figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [600 200 400 300], 'Name', 'Grand Choice Probability');
   Zfar_grand = []; Znear_grand = [];
   %combine data across correlations into grand distributions
   for i=1:num_disp
      %include in grand distribution only if the monkey made choices in both direction: at least 25% of the trials in either direction.
      if (min(length(Z_far_dist{i}),length(Z_near_dist{i})) / max(length(Z_far_dist{i}),length(Z_near_dist{i})) > 1/3)
      	Zfar_grand = [Zfar_grand Z_far_dist{i}];   
         Znear_grand = [Znear_grand Z_near_dist{i}];   
         
      end
   end
   
   PlotTwoHists(Zfar_grand, Znear_grand);
   
   %do permutation test to get P value for grand CP
   [grandCP, grandPval] = ROC_signif_test(Zfar_grand, Znear_grand);
   titl = sprintf('grand CP = %5.2f, P = %6.4f', grandCP, grandPval);
   title(titl);
   
   str = sprintf('%s vergence', FILE );
   xlabel(str);

   %Now, make a scatterplot of vergence angles sorted by choice for calibrated and uncalibrated data TU, 01/10/01
	near_dist_h = []; near_dist_v = []; near_dist_h = []; near_dist_v = [];
   if (Pref_HDisp - unique_hdisp_p4(1) < 0)
		near_choices = (choice == PREFERRED);
         near_dist_h = vergence_h(near_choices & select_trials);
         near_dist_v = vergence_v(near_choices & select_trials);
         near_dist_h = vergence_h(near_choices & select_trials);
         near_dist_v = vergence_v(near_choices & select_trials);
	      far_choices = (choice == NULL);
         far_dist_h = vergence_h(far_choices & select_trials);
         far_dist_v = vergence_v(far_choices & select_trials);
         far_dist_h = vergence_h(far_choices & select_trials);
         far_dist_v = vergence_v(far_choices & select_trials);
   elseif (Pref_HDisp - unique_hdisp_p4(1) > 0)
         far_choices = (choice == PREFERRED);
         far_dist_h = vergence_h(far_choices & select_trials);
         far_dist_v = vergence_v(far_choices & select_trials);
         far_dist_h = vergence_h(far_choices & select_trials);
         far_dist_v = vergence_v(far_choices & select_trials);
	      near_choices = (choice == NULL);
         near_dist_h = vergence_h(near_choices & select_trials);   
         near_dist_v = vergence_v(near_choices & select_trials);
         near_dist_h = vergence_h(near_choices & select_trials);
         near_dist_v = vergence_v(near_choices & select_trials);
   end
      
   figure; %Scatter plot using uncalibrated data
   set(gcf, 'PaperPosition', [.2 .2 8 10.7], 'Position',[600 100 400 300], 'Name', 'Uncalibrated Vergence Plot');  
   hold on;
   Handl(1) = plot(near_dist_h, near_dist_v, 'ko', 'MarkerFaceColor', 'k');
   Handl(2) = plot(far_dist_h, far_dist_v, 'ko');
   hold off;
   
   	xlabel('Horizontal Disparity (deg)');
	ylabel('Vertical Disparity (deg)');
   legend(Handl, 'Near', 'Far', 2);
      
   figure; %Scatter plot using calibrated data
   set(gcf, 'PaperPosition', [.2 .2 8 10.7], 'Position',[600 100 400 300], 'Name', 'Calibrated Vergence Plot');  
   hold on;
   Handl(1) = plot(near_dist_h, near_dist_v, 'ko', 'MarkerFaceColor', 'k');
   Handl(2) = plot(far_dist_h, far_dist_v, 'ko');
   hold off;
   
   	xlabel('Horizontal Disparity (deg)');
	ylabel('Vertical Disparity (deg)');
    legend(Handl, 'Near', 'Far', 2);      
    
    %now, Z-score the spike rates for each bin_corr and disparity condition
    Z_Spikes = spike_rates;
    for i=1:length(unique_hdisp_p1)
        select = (h_disp_p1 == unique_hdisp_p1(i));
        z_dist = spike_rates(select);
        z_dist = (z_dist - mean(z_dist))/std(z_dist);
        Z_Spikes(select) = z_dist;
    end

    %[binoc_corr_ranks' hdisp_ranks' choice' spike_rates' h_verg']

    %run the ANCOVA analysis tool
    aoctool(vergence_h', Z_Spikes', choice);

return;
   