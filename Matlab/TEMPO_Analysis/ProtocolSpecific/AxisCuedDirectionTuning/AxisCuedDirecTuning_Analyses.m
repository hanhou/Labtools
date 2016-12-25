function AxisCuedDirecTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

%PATH = 'Z:\Data\Tempo\Robbins\Raw\'
%FILE = 'm3c1165r112.htb'
switch(Analysis{1})
case 'Plot Psychometric'
  PsychAxisCuedDirec(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
end

return;