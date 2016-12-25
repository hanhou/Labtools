function Stereoacuity_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	switch(Analysis{1})
	   case 'Plot Neurometric/Psychometric'
         NeuroPsychoPlot_st(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      case 'Response Distributions'
         PlotResponseDistributions_st(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      case 'Plot Rasters/Histograms'
         PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, Protocol);
      case 'Compute Choice Probabilities'
         Compute_ChoiceProb_st(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      case 'Compute Choice Probabilities from Vergence'
         Compute_ChoiceProb_verg_st(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      case 'Plot Microstim'
         Microstim_st(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      case 'Behavioral Bootstrap'
         BehavBootstrap(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      
    end

return;