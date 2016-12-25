function HGradTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

switch(Analysis{1})
    case 'Plot Grad Angle Tuning Curve'
        HGradTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Simulate Gradient Response'
        HGradModel(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    case 'Fit with Sine Wave'
        HGradSinFit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    case 'Seq F-Test for Fixed/Unfixed Freq in Sin Fit'
        HGradSinFit_FreqTest(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    case 'Fit with Multiple Sine Waves'
        HGradMultiSinFit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    case 'Simultaneous Fit with Phase Slop'
        HGradFitPhaseSlop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    case 'Simulate Heterogeneous Surround'
        Surround_Model(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    case 'Sorted Spike Correlogram'
        HGradSpikeCorr(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, PATH, FILE);
    case 'Plot Slant Tuning Curve'
        SlantTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Vergence Data' 
        disp('plot vergence here');   
end

return;