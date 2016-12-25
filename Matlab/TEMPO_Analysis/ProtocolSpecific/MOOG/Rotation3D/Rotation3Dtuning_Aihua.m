function Rotation3Dtuning_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% addpath('Z:\Users\Aihua\CurveFittingAnalysis')

% OutputPath=['Z:\Users\Aihua\z_tempOutputs\'];
% OutputPath=['Z:\Users\Tanya\z_tempOutputs\'];
OutputPath=['C:\Aihua\z_TempOutputs\'];

% SaveTrials(data, Protocol, Analysis, SpikeChan, StartCode, StopCode,BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath)
[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);

%*****************************************************************%
%spatial-temporal tuning analysis
% addpath('Z:\Users\Aihua\PSTHAnalysis');
% MOOG_PSTH_Excit_Inhibit(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% Rotation_PSTH_save(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_PSTH_cah_LoadVersion(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_PSTH_cah(data, 1, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_PSTH_cah_LoadVersion(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% PSTHEnlargeWhole(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol,OutputPath);
% TuningPeakTimeEnlarge(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_ModeStep_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);

%*****************************************************************%
%CurveFitting Analysis
% addpath('Z:\Users\Aihua\CurveFittingAnalysis')
addpath('Z:\Users\Aihua\STCurveFitting')
% CurveFitting_3Plane(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
CurveFitting_3Plane_Sti(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
% CurveFitting_3Plane_Sti_Astart(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
% CurveFitting_3Plane_FEF(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
% CurveFitting_3Plane_NDC2(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
% CurveFitting_3Plane_DC2_Compare(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
%  CurveFitting_Hor(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);



return;