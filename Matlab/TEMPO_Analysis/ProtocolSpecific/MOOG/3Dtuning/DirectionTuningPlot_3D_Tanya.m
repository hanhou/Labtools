% DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	YONG, 6/27/03  
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_3D_Tanya(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% addpath('Z:\Users\Aihua\PSTHAnalysis');
addpath('Z:\Users\Aihua\CurveFittingAnalysis');

% OutputPath=['Z:\Users\Aihua\z_tempOutputs\'];
% OutputPath=['Z:\Users\Tanya\z_tempOutputs\'];
OutputPath=['C:\z_tempoutputs\'];
% OutputPath=['C:\Fisher\'];

 added4RecoveringSpikeData;

% SaveTrials(data, Protocol, Analysis, SpikeChan, StartCode, StopCode,BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath)
[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
% MOOG_PSTH_Excit_Inhibit(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_PSTH_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MeanVar=MOOG_MeanVar_PeakTime(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% CurveFitting_3Plane(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
%   CurveFitting_3Plane_FEF(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)% 
% CurveFitting_Hor(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% CurveFitting_3Plane_Sti(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% CurveFitting_3Plane_Sti_Tania(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
 CurveFitting_3Plane_Sti_Astart(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% PSTHEnlargeWhole(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol,OutputPath);
%%%%%%%%%%%%%%%%%%%
%  CurveFitting_3Plane_Tanya(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)
% MOOG_TuningStep_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_ModeStep_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);

% MOOG_PSTH_cah_LoadVersion(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_PSTH_Tanya_LoadVersion(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% MOOG_TuningStep_temp_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);
% PSTHEnlargeWhole(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol,OutputPath);
% TuningPeakTimeEnlarge(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);

% MOOG_Translation_PSTH_LoadVersion(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Translation_PSTH_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Translation_PSTH_cah_try(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Translation_TuningStep_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Translation_TuningStep_cah_try(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Translation_corr_cah(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% MOOG_Translation_sep(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

return;
