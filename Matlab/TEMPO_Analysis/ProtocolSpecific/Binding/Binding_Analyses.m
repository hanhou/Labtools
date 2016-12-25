function Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

	switch(Analysis{1})
    case 'Plot Direction Tuning Curve'
        BindingDirectionTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);  
    case 'Plot Spatial Frequency Tuning Curve'
        BindingSpatFreqTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);  
    case 'Plot Temporal Frequency Tuning Curve'
        BindingTempFreqTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);  
    case 'Plot Neurometric/Psychometric'
        NeuroPsychoPlot_bnd(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
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
    case 'Plot Joint PSTH'
        JointPSTH(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    case 'Plot Event Times'
        EventTimes(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot ISI Histogram'
        ISIHist(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Spike-Triggered Averages' 
        PlotSTSA(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);    
    case 'Plot Eye Traces'
        PlotEyeTraces(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot LFPs'
        PlotLFPsJustin(data, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rasters'
        PlotRasters(data, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Fit with 2D Gaussian'
        BindingGauss2D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Correlogram Timecourse'
        CrossCorrTimeCourse(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    case 'Plot Coherence Spectrogram'
        Compute_Coherence(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    case 'Simulate Spike Trains'
        SimulateSpikeTrains(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   
    end

return;