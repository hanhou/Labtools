function CuedDirectionDiscrim_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

switch(Analysis{1})
    case 'Plot Psychometric'
        PsychCuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1d Psychometric'
        Psych1dCuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Bootstrap Curve Fits'
        bootstrap_curve_fits(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Neurometric/Psychometric'
        NeuroPsychoCuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot N/P sorted by CueDir'
        NPCuedDirec_CueDir(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Cue Direction Effects'
        CueDirecEffects(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'PSTH'
        PSTH_CuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'LFP PSTH'
        LFP_PSTH_CuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Compute CP'
        CP_CuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Psych SepDir Curves'
        Psych1d_sepDir_curves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Neural Response Functions'
        NeuralResponse_CuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
end