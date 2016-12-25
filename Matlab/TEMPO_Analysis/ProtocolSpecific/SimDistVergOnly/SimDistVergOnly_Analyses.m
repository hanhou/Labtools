function SimDistVergOnly_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

switch(Analysis{1})
case 'Plot Vergence Data'
    SimDistVergence(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
end

return;