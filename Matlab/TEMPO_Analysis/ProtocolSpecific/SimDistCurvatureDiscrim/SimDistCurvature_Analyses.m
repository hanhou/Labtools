function SimDistCurvature_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

switch(Analysis{1})
case 'Plot Psychometric'
    CurvaturePsychometric(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
end

return;