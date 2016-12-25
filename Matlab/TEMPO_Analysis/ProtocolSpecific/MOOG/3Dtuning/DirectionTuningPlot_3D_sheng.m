% DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	YONG, 12/10/08  
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_3D_sheng(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

% %%%% for variance to mean
% [StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
%   addpath 'Z:\Users\sheng\program\variancemean';
%  OutputPath =  'Z:\Users\sheng\program\variancemean\'
%  getmoogvarmeandata(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath);



%%%for heading temporal
addpath 'Z:\Users\sheng\program\3dsimple';

% OutputPath =  'Z:\Users\Aihua\z_tempOutputs\';   
OutputPath = 'Z:\Users\sheng\choice\input\';

[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
    
if Protocol == 100 | Protocol == 107  | Protocol == 104  
    addpath('Z:\Users\Aihua\HeadingDiscriAnalysis');
%     addpath 'Z:\Users\sheng\program\3dsimple';
%     SysDelay=115;    StartEventBin=StartEventBin+SysDelay;    StopEventBin=StopEventBin+SysDelay;
%     Search_StartPoint=400
%     [SpikeArray PSTH_Smooth Step WindowInterval SmoothBin]=PSTH_TemporalModulation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Search_StartPoint,OutputPath);        
%     PeakTime_hor(FILE,SpikeArray,PSTH_Smooth,OutputPath,StartEventBin, StopEventBin,Step,WindowInterval,SmoothBin,Search_StartPoint,OutputPath, SpikeChan);            
WindowInterval = 400;             
[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);
            f3dsimple_temporal(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath,  WindowInterval);
else
    addpath 'Z:\Users\sheng\program\heading';
    HeadingDis_cum_origin_temporal(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath);

    %     HeadingDis_cum_origin_afferent(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, OutputPath); 
end

return;
