function MovingTarget_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

switch(Analysis{1})
case 'Fit with Gaussian'
  MovTargFit(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
end

return;