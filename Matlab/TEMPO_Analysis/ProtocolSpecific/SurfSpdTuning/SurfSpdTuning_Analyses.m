function SurfTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

switch(Analysis{1})
    case 'Plot Tuning Curve'
        SurfSpdTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
end

return;