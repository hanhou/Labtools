function SurfTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

switch(Analysis{1})
    case 'Fit with Ind. Sin Fxn'
        Surf_IndSinFit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'Plot Tuning Curve'
        SurfOrientTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Vector Average'
        SurfVectorAvg(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Vector Average normalized by STD'
        SurfVectorAvgSTD(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Vector Average with Tuning Bias (STD)'
        SurfVectorAvgSTD_TB(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tilt Slant Tuning Curves'
        SurfSlantTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
end

return;