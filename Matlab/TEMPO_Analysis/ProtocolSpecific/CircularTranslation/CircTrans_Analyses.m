function CircTrans_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
	switch(Analysis{1})
    case 'Response Distributions'
        PlotResponseDistributions(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot PSTH'
        PSTH(data, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Compute Choice Probabilities'
        Compute_ChoiceProb(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Cross Correlograms'
        CrossCorr(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    case 'Plot Spike-Triggered Averages'
        PlotSpikeTA(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    case 'Compare Correlation Levels',
        BindingPermutationTest(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    case 'Plot Auto Correlograms'
        AutoCorr(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    case 'Plot Event Times'
        EventTimes(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Eye Traces'
        PlotEyeTraces(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot LFPs'
        PlotLFPs(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rasters'
        PlotRasters(data, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Circ Trans Tuning Curve'
        CircTransTuningCurve(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    
    end
return;