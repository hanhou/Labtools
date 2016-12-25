%-----------------------------------------------------------------------------------------------------------------------
%-- PlotResponseDistributions.m -- Plots neural response histograms, sorted by disparity, at each binoc. correlation level
%--	GCD, 6/5/00
%-----------------------------------------------------------------------------------------------------------------------
function PlotResponseDistributions(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
   
	%get the column of values of horiz. disparities in the dots_params matrix
   h_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
   unique_hdisp = munique(h_disp');
   
   %get the binocular correlations
   binoc_corr = data.dots_params(DOTS_BIN_CORR, :, PATCH1);
   unique_bin_corr = munique(binoc_corr');
   
   %now, get the firing rates for all the trials 
   spike_rates = data.spike_rates(SpikeChan, :);

   %get indices of any NULL conditions (for measuring spontaneous activity
   null_trials = logical( (binoc_corr == data.one_time_params(NULL_VALUE)) );
   
   %now, select trials that fall between BegTrial and EndTrial
   trials = 1:length(binoc_corr);		% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
   
   %[h_disp' binoc_corr' spike_rates' null_trials' select_trials']
   
   Pref_HDisp = data.one_time_params(PREFERRED_HDISP);
   
   %now, plot the spike distributions, sorted by disparity, for each correlation level
   figure;
	set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 250 500 573], 'Name', 'Responses Sorted by Disparity');
   PREFERRED = 1; NULL = 2;
   num_corrs = length(unique_bin_corr);
   ROC_values = []; pref_dist = []; null_dist = [];
   for i=1:num_corrs
      subplot(num_corrs, 1, i);
	   MAX_REPS = 100;
   	Mtemp = ones(MAX_REPS, 2)*NaN;
      pref_trials = ( (h_disp == Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );
      pref_dist{i} = spike_rates(pref_trials & select_trials);
      Mtemp(1:length(pref_dist{i}), PREFERRED) = pref_dist{i}';
      null_trials = ( (h_disp ~= Pref_HDisp) & (binoc_corr == unique_bin_corr(i)) );
      null_dist{i} = spike_rates(null_trials & select_trials);
      Mtemp(1:length(null_dist{i}), NULL) = null_dist{i}'; 
      
      ROC_values(i) = rocN(pref_dist{i}, null_dist{i}, 100)
      
      hist(Mtemp);
   end   
   
   %now, print out spike rates for all trials at each correlation, sorted by disparity
  	i = size(PATH,2) - 1;
   while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
     	i = i - 1;
  	end   
  	PATHOUT = [PATH(1:i) 'Analysis\NeuroPsychoCurves\'];
  	i = size(FILE,2) - 1;
  	while FILE(i) ~='.'
     	i = i - 1;
  	end
  	FILEOUT = [FILE(1:i) 'resp_dists'];
     
  	fileid = [PATHOUT FILEOUT];
  	fwriteid = fopen(fileid, 'w');
     
   len = [];  
  	for i=1:num_corrs
      temp = [length(pref_dist{i}) length(null_dist{i})];
      len(i) = max(temp);
   end
   len
   max_vals = max(len);
   min_vals = min(len);
   
  	for i=1:num_corrs
   	fprintf(fwriteid,'P%.1f\tN%.1f\t', unique_bin_corr(i), unique_bin_corr(i));  
   end
   fprintf(fwriteid, '\n');
   
   for i=1:min_vals
      for j=1:num_corrs
         fprintf(fwriteid, '%6.3f\t%6.3f\t', pref_dist{j}(i), null_dist{j}(i));
      end
	   fprintf(fwriteid, '\n');
   end
   
   for i=min_vals+1:max_vals
      for j=1:num_corrs
	      if (length(pref_dist{j}) >= i)
   	      fprintf(fwriteid, '%6.3f\t', pref_dist{j}(i));
      	else
         	fprintf(fwriteid, '\t');
         end
         
	      if (length(null_dist{j}) >= i)
   	      fprintf(fwriteid, '%6.3f\t', null_dist{j}(i));
      	else
         	fprintf(fwriteid, '\t');
         end
         
	   fprintf(fwriteid, '\n');         
      end
   end
   
   fclose(fwriteid);
   
   
   
return;