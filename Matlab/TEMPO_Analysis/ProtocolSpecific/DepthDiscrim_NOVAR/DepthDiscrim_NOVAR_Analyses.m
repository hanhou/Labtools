function DepthDiscrim_NOVAR_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

switch(Analysis{1})
case 'Plot Psychometric Var/Novar'
    NOVAR_PsychoPlot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
end

return;