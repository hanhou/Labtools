function CurveFitting_3Plane(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)

warning off MATLAB:divideByZero;
warning off MATLAB:singularMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;clc
% load time_vel_acc.mat;%load('C:\MATLAB6p5\work\time_vel_acc.mat');
% load m14c136r1.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
switch Protocol
    case 100 % DIRECTION_TUNING_3D 
        temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
        temp_elevation = data.moog_params(ELEVATION,:,MOOG);
        CurrentProtocol=['Translation'];
    case 112 %ROTATION_TUNING_3D
        temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
        temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
        CurrentProtocol=['Rotation'];
    case 104 %DIR3D_VARY_FIXATION 
        temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
        temp_elevation = data.moog_params(ELEVATION,:,MOOG);
        CurrentProtocol=['DIR3D VARY FIXATION '];
end

temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_spike_data = data.spike_data(SpikeChan,:);
temp_spike_rates = data.spike_rates(SpikeChan, :); 

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_trials ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

% find spontaneous trials which azimuth,elevation,stim_type=-9999
azimuth = temp_azimuth(~null_trials & select_trials);unique_azimuth = munique(azimuth');
elevation = temp_elevation(~null_trials & select_trials);unique_elevation = munique(elevation');
stim_type = temp_stim_type(~null_trials & select_trials);unique_stim_type = munique(stim_type');
spike_rates= temp_spike_rates(~null_trials & select_trials);
% spike_rates_step=temp_SpikeRates(:,~null_trials & select_trials);

condition_num = stim_type;unique_condition_num = munique(condition_num');
h_title{1}='Vestibular';h_title{2}='Visual';h_title{3}='Combined';
stim_duration = length(temp_spike_data)/length(temp_azimuth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove null trials, bad trials, and trials outside Begtrial~Endtrial for spike_data
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    %temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :Discard_trials(i)*stim_duration ) = 99;
end
spike_data = temp_spike_data( temp_spike_data~=99 );
spike_data( find(spike_data>89) ) = 1; % something is absolutely wrong 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the SpikeArray
trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)) * length(unique_stim_type)+1;
repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep); 

unique_azimuth0=[0:45:315]';
for k=unique_condition_num(1):unique_condition_num(end)%1: length(unique_condition_num)
    for j=1:length(unique_elevation)      
        for i=1: length(unique_azimuth0)
            clear select; select = logical( (azimuth==unique_azimuth0(i)) & (elevation==unique_elevation(j)) & (condition_num==k) );
%             clear select; select = logical( (azimuth==unique_azimuth0(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );
            act_found = find( select==1 );
            if length(act_found)>repetitions
                NumRepetition=repetitions;
            else
                NumRepetition=length(act_found);
            end
            for repeat=1:NumRepetition
                spikeTracker(repeat,:)=spike_data(1,stim_duration*(act_found(repeat)-1)+1:stim_duration*(act_found(repeat)));
            end
            SpikeArray{k,j,i}(:,:)=spikeTracker;
            SpikeHistArray{k,j,i}=mean(spikeTracker);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Step = 100;%Step = 20;
WindowInterval=100;%WindowInterval=50;
WindowInterval_Big=400;
%AftSti=100;%
AftSti=500;
StepEventBin1=[StartEventBin(1,1):-Step:StartEventBin(1,1)-500];StepEventBin1=fliplr(StepEventBin1);StepEventBin1=StepEventBin1(1:end-1);
StepEventBin2=[StartEventBin(1,1):Step:StopEventBin(1,1)+AftSti];
StepEventBin=[StepEventBin1 StepEventBin2];
StartIndex=find(StepEventBin==StartEventBin(1,1));
[min_diff,EndIndex]=min(abs(StepEventBin-StopEventBin(1,1)));

for k=unique_condition_num(1):unique_condition_num(end)
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth0)      
            clear StartBinReset;StartBinReset=StepEventBin(1,1);%StartBinReset=StartEventBin(1,1)-500;
            clear Index; Index=1;
            clear tempDC; tempDC=SpikeArray{k,j,i}(:,StartEventBin(1,1)-100:StartEventBin(1,1)+300);
            clear DCMean; DCMean=sum(tempDC')*1000/400;
            
            tempPSTH=[];            
            while StartBinReset<StopEventBin(1,1)+AftSti                
                clear CurrentValue;CurrentValue=SpikeArray{k,j,i}(:,StartBinReset-0.5*WindowInterval+1:StartBinReset+0.5*WindowInterval); 
                clear CurrentValueMean; CurrentValueMean=sum(CurrentValue')*1000/WindowInterval;                 %CurrentValueMean=mean(CurrentValue');                 
                StepMatrix{k}(j,i,Index)=mean(CurrentValueMean);
                StdMatrix{k}(j,i,Index)=std(CurrentValueMean);

                clear CurrentValue;CurrentValue=SpikeArray{k,j,i}(:,StartBinReset-0.5*WindowInterval+1:StartBinReset+0.5*WindowInterval);
                clear CurrentSpikeCount; CurrentSpikeCount=sum(CurrentValue')*1000/WindowInterval;
%                 SpikeCount_Trial{k,Index}(j,i,:)=CurrentSpikeCount;
                SpikeCount_Trial{k,j,i}(Index,:)=CurrentSpikeCount;                
                clear CurrentValue;CurrentValue=SpikeArray{k,j,i}(:,StartBinReset-0.5*WindowInterval_Big+1:StartBinReset+0.5*WindowInterval_Big);
                clear CurrentValueMean; CurrentValueMean=sum(CurrentValue')*1000/WindowInterval;                 
                tempPSTH(1,Index)=mean(CurrentValueMean);
                
                Index=Index+1;                
                StartBinReset=StartBinReset+Step;
            end                
            tempPSTH(1:StartIndex-1)=NaN;tempPSTH(EndIndex:end)=NaN;%tempPSTH(Index:end)=NaN;
            clear MaxIndex; MaxIndex=find(tempPSTH==max(tempPSTH));
            clear PeakTimeIndex; PeakTimeIndex=StepEventBin(MaxIndex(1,1));            
            clear PeakValues; PeakValues=SpikeArray{k,j,i}(:,PeakTimeIndex-0.5*WindowInterval_Big+1:PeakTimeIndex+0.5*WindowInterval_Big);
            clear PeakValueMean; PeakValueMean=sum(PeakValues')*1000/WindowInterval;  
            Value_peak{k}(j,i)=mean(PeakValueMean); 
            Time_peak{k}(j,i)=PeakTimeIndex*Step/1000-0.5;
            TimeIndex_peak{k}(j,i)=MaxIndex(1,1);
            
            [temp_p_peak(1),h_peak] = ranksum(PeakValueMean,DCMean,0.01); clear h_peak;            
            clear PeakValues; PeakValues=SpikeArray{k,j,i}(:,PeakTimeIndex-0.5*WindowInterval_Big+1-Step:PeakTimeIndex+0.5*WindowInterval_Big-Step);
            clear PeakValueMean; PeakValueMean=sum(PeakValues')*1000/WindowInterval_Big;            
            p_peak{k}(j,i)=max(temp_p_peak);
            
            clear MinIndex; MinIndex=find(tempPSTH==min(tempPSTH)); 
            clear TroughTimeIndex; TroughTimeIndex=StepEventBin(MinIndex(1,1));            
            clear TroughValues; TroughValues=SpikeArray{k,j,i}(:,TroughTimeIndex-0.5*WindowInterval_Big+1:TroughTimeIndex+0.5*WindowInterval_Big);
            clear TroughValueMean; TroughValueMean=sum(TroughValues')*1000/WindowInterval_Big;   
            Value_trough{k}(j,i)=mean(TroughValueMean); 
            Time_trough{k}(j,i)=TroughTimeIndex*Step/1000-0.5;
            TimeIndex_trough{k}(j,i)=MinIndex(1,1);
            
            [temp_p_trough(1),h_Trough] = ranksum(TroughValueMean,DCMean,0.01); clear h_Trough;            
            clear TroughValues; TroughValues=SpikeArray{k,j,i}(:,TroughTimeIndex-0.5*WindowInterval_Big+1-Step:TroughTimeIndex+0.5*WindowInterval_Big-Step);
            clear TroughValueMean; TroughValueMean=sum(TroughValues')*1000/WindowInterval_Big;                      
            p_trough{k}(j,i)=max(temp_p_trough);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Azi_3D=[0:45:315];Ele_3D=[-90:45:90];%[AziGrid,EleGrid]=meshgrid(UniAzi,UniEle);  
% OutputPath=['C:\Aihua\z_TempOutputs\'];% OutputPath=['Z:\Users\Aihua\z_tempOutputs\'];
%%%******************************%%%
%xy: Horizontal plane
Azi_hor=[0:45:315];Ele_hor=0*ones(1,8);
CurveFit_Hor=CuveFitting_Plot(FILE,'Horizontal Plane',Azi_3D,Ele_3D,Azi_hor,Ele_hor,Step,StepMatrix,StdMatrix,SpikeCount_Trial,p_peak,Value_peak,TimeIndex_peak,p_trough,Value_trough,TimeIndex_trough);
%Save the figure 
FigureIndex=2;
if ~isempty(CurveFit_Hor)
    figure(FigureIndex); set(gcf, 'PaperOrientation', 'portrait');
    saveas(gcf,[OutputPath FILE(1:end-4) '_CurveFitting_Hor.png'],'png');close(FigureIndex);
    %Save the Data
    SaveFileName=[OutputPath FILE(1:end-4) '_CurveFitHor'];
    save(SaveFileName,'CurveFit_Hor'); clear SaveFileName;
else
end

%%%******************************%%%
%xz: Frontal plane
Azi_Front=[0 0 0 180 180 180 0 0];Ele_Front=[0 45 90 45 0 -45 -90 -45];
CurveFit_Front=CuveFitting_Plot(FILE,'Frontal Plane',Azi_3D,Ele_3D,Azi_Front,Ele_Front,Step,StepMatrix,StdMatrix,SpikeCount_Trial,p_peak,Value_peak,TimeIndex_peak,p_trough,Value_trough,TimeIndex_trough);
if ~isempty(CurveFit_Front)
    %Save the figure 
    figure(FigureIndex); set(gcf, 'PaperOrientation', 'portrait');
    saveas(gcf,[OutputPath FILE(1:end-4) '_CurveFitting_Front.png'],'png');close(FigureIndex);
    %Save the Data
    SaveFileName=[OutputPath FILE(1:end-4) '_CurveFitFront'];
    save(SaveFileName,'CurveFit_Front'); clear SaveFileName;
else
end

%%%******************************%%%
%yz: Median plane (or mid-sagital plane)
Azi_Med=[90 90 0 270 270 270 0 90];Ele_Med=[0 45 90 45 0 -45 -90 -45];
CurveFit_Med=CuveFitting_Plot(FILE,'Median Plane',Azi_3D,Ele_3D,Azi_Med,Ele_Med,Step,StepMatrix,StdMatrix,SpikeCount_Trial,p_peak,Value_peak,TimeIndex_peak,p_trough,Value_trough,TimeIndex_trough);
if ~isempty(CurveFit_Med)
    %Save the figure 
    figure(FigureIndex); set(gcf, 'PaperOrientation', 'portrait');
    saveas(gcf,[OutputPath FILE(1:end-4) '_CurveFitting_Med.png'],'png');close(FigureIndex);
    %Save the Data
    SaveFileName=[OutputPath FILE(1:end-4) '_CurveFitMed'];
    save(SaveFileName,'CurveFit_Med'); clear SaveFileName;
else
end








