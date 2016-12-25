function MunitAnalysis(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

addpath('Z:\Users\Xiaodong\data_ana\PSTHAnalysis\');
OutputPath='Z:\Users\Xiaodong\data_ana\mu_PSTH\';

[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
MOOG_PSTH_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);

return