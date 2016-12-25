%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_Fix.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-- Modified for Vary_Fixation protocol: Can now handle multiple stim types and multiple gaze angles.  CRF + Yong, 12/03
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_Fix(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

addpath('Z:\Users\Aihua\PSTHAnalysis');
addpath('Z:\Users\Aihua\CurveFittingAnalysis')

% OutputPath=['Z:\Users\Aihua\Temporal_OutputData\'];
OutputPath=['C:\Aihua\z_TempOutputs\'];

% SaveTrials(data, Protocol, Analysis, SpikeChan, StartCode, StopCode,BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath)
[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
% VarMean=MOOG_MeanVar_PeakTime(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% CurveFitting_3Plane(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
% MOOG_PSTH_Excit_Inhibit(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MeanVar=MOOG_MeanVar_cah1(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MeanVar=MOOG_MeanVar_cah2(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MeanVar=MOOG_MeanVar_cah3(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MeanVar=MOOG_MeanVar_cah4(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_PSTH_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_PSTH_cah_LoadVersion(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% PSTHEnlargeWhole(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol,OutputPath);
% TuningPeakTimeEnlarge(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
%--------------------------------------------------------------------%
%CurveFitting Analysis
% addpath('Z:\Users\Aihua\CurveFittingAnalysis')
addpath('Z:\Users\Aihua\STCurveFitting')
% CurveFitting_3Plane(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
CurveFitting_3Plane_Sti(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
