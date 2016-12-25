function test()
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
%get the relevant spike_rates 
Step=500;%running window
WindowInterval=500;% ms (timebin for plot PSTH)
clear StartBinReset;StartBinReset=StartEventBin(1,1)-500;

clear Index;Index=1;
while StartBinReset<StopEventBin(1,1)+500 %StartBinReset<StopEventBin(1,1)-0.5*WindowInterval
    StartOffsetBin=StartBinReset+1-StartEventBin(1,1)-0.5*WindowInterval;
    StopOffsetBin=(StartBinReset+0.5*WindowInterval)-StopEventBin(1,1);%StopOffsetBin=StopEventBin(1,1)-(StartBinReset+WindowInterval);    
    clear TempSpikeRates;TempSpikeRates = ComputeSpikeRates(data, size(data.event_data, 3), StartCode, StopCode, StartOffsetBin, StopOffsetBin);
    temp_SpikeRates(Index,:)=TempSpikeRates(SpikeChan,:);
    Index=Index+1;
    StartBinReset=StartBinReset+Step;
end
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
spike_rates_step=temp_SpikeRates(:,~null_trials & select_trials);

condition_num = stim_type;unique_condition_num = munique(condition_num');
h_title{1}='Vestibular';h_title{2}='Visual';h_title{3}='Combined';
stim_duration = length(temp_spike_data)/length(temp_azimuth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeSeq=[-500:Step:2500];
StartIndex=find(timeSeq==0);
EndIndex=find(timeSeq==2500);
current_sti=1;
%xy: Horizontal plane
temp_Azi=[0:45:315];temp_Ele=0*ones(1,8);
% %xz: Frontal plane
% temp_Azi=[0 0 0 180 180 180 0 0];temp_Ele=[0 45 90 45 0 -45 -90 -45];
%yz: Median plane (or mid-sagital plane)
% temp_Azi=[90 90 0 270 270 270 0 90];temp_Ele=[0 45 90 45 0 -45 -90 -45];

spacetime_data=[];
for i=1:length(temp_Azi)    
    clear select; select=find(stim_type==current_sti & azimuth==temp_Azi(i) & elevation==temp_Ele(i));
    clear tempdata; tempdata=spike_rates_step(StartIndex:EndIndex,select);tempdata=tempdata';
    spacetime_data(i,:)=nanmean(tempdata); 
    spacetime_data_s(i,:)=nanmean(sqrt(tempdata)); 
    spacetime_std(i,:)=std(sqrt(tempdata));
    spacetime_err(i,:)=std(sqrt(tempdata))/sqrt(length(select));
end
spacetime_data_s(i+1,:)=spacetime_data_s(1,:); 
spacetime_data(i+1,:)=spacetime_data(1,:); 
spacetime_std(i+1,:)=spacetime_std(1,:);
spacetime_err(i+1,:)=spacetime_err(1,:);
figure;contourf(spacetime_data');colorbar
title(['bin = ' num2str(Step) 'ms']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = ([0 45 90 135 180 225 270 315 360]*pi/180)'; % for fitting
for i = 1:size(spacetime_data,2)%101
    xdata(:, i) = x;    
end
% t=0:Step*0.001:0.02*(size(spacetime_data,2)-1);%t = 0:0.02:2;  %% 190 bins -- ingore the last bins with zeros
tdata=[];
t=timeSeq(StartIndex:EndIndex)*0.001;
for i = 1:size(spacetime_data,1)
    tdata(i,:) = t;
end
xtdata = [xdata; tdata]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FigureIndex=2; 
figure(FigureIndex);set(FigureIndex,'Position', [50,100 1200,800], 'Name', 'CurveFitting');orient landscape; 
% text(-0.1,1.06,[FILE '  ' Plane]); axis off;    
clear XAzi Ytime;[XAzi,Ytime] = meshgrid([0:45:360],[min(t):Step*0.001:max(t)]);
figure(FigureIndex);axes('position',[0.05 0.77 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,spacetime_data');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');title('raw data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model1: Velocity only
allow_negative=1;
global model_use
model_use=1;
% vect_Vel=[47.2 92.3 0.8 340.2*pi/180 0.904 0.166 0.489];
% spacefit_Vel= funccosnlin(vect_Vel,xtdata);
[spacefit_Vel,vect_Vel,r_squared_Vel,CI_Vel] = MSF_Vel_fit(spacetime_data,Step*0.001,allow_negative);
% clear model_Vel; model_Vel=zeros(size(PSTH_tempPlane));
% model_Vel(:,1:25) = NaN; model_Vel(:,26:151)=spacefit_Vel(1:8,:);

figure(FigureIndex);axes('position',[0.05 0.53 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_Vel',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
%xlabel('Azimuth, X (deg)');  
ylabel('Time (sec)');title('model: Vel  ');

error_surf = spacetime_data - spacefit_Vel;
% current_err = cosnlin_err(vect_Vel);
err_Vel = cosnlin_err(vect_Vel);
figure(FigureIndex);axes('position',[0.26 0.53 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,error_surf'); colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');
axis off;
title({[ 'Err: ' num2str(err_Vel, '%0.2f') ]},  'FontSize', 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do the chi-square goodness fit;
spike_rates_step_n=sqrt(spike_rates_step);
for i=1:length(temp_Azi)    
    clear select; select=find(stim_type==current_sti & azimuth==temp_Azi(i) & elevation==temp_Ele(i));
    clear tempdata; tempdata=spike_rates_step_n(StartIndex:EndIndex,select);tempdata=tempdata';
    spacetime_data_n(i,:)=nanmean(tempdata);    
    spacetime_std(i,:)=std(tempdata);
    spacetime_err(i,:)=std(tempdata)/sqrt(length(select));
end
spacetime_data_n(i+1,:)=spacetime_data(1,:); 
spacetime_std(i+1,:)=spacetime_std(1,:); 

% Do chi-square goodness of fit test
[spacefit_Vel,newvect] = funccosnlin(vect_Vel,xtdata); 
spacefit_Vel(spacefit_Vel<0)=0;% err = sqrt(sum( sum(( spacefit_Vel- spacetime_data) .^2) ));
clear DataDiff;DataDiff=spacetime_data_n-sqrt(spacefit_Vel);%DataDiff=sqrt(spacetime_data)-sqrt(spacefit_Vel);
Data_std=spacetime_std;%Data_std=sqrt(spacetime_std);
clear sd_zero;sd_zero=find(Data_std==0);Data_std(sd_zero)=NaN;
dof_Vel=(size(spacetime_data,1)-1)*(size(spacetime_data,2)-1)-length(sd_zero)-length(vect_Vel);%degree of freedom
chi2_Vel=sum(nansum( (DataDiff./Data_std).^2));
RChi2_Vel = chi2_Vel/dof_Vel;%Reduced Chi-Squared value
chi2P_Vel=1-chi2cdf(chi2_Vel,dof_Vel);

[chi2_Vels, chi2P_Vels] = Chi2_Test_3D(azimuth, spike_rates_step, 'funccosnlin', vect_Vel, length(vect_Vel))];
% [chi2, chiP] = Chi2_Test_3D(datax, datay, funcname, params, num_free_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_use=2;
[spacefit_VelAcc,vect_VelAcc,r_squared_VelAcc,CI_VelAcc] = MSF_VelAcc_fit(spacetime_data, vect_Vel,Step*0.001,allow_negative);
% vect_VelAcc=[48.1 101.1 0.4 175.2*pi/180 1.238 0.259 0.160 0 166.9*pi/180];
% spacefit_VelAcc= funccosnlin(vect_VelAcc,xtdata);
% clear model_VelAcc; model_VelAcc=zeros(size(PSTH_tempPlane));
% model_VelAcc(:,1:25) = NaN; model_VelAcc(:,26:151)=spacefit_VelAcc(1:8,:);

figure(FigureIndex);axes('position',[0.05 0.31 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_VelAcc');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
%xlabel('Azimuth, X (deg)');  
ylabel('Time (sec)');title('model: Vel + Acc ');

error_surf = spacetime_data - spacefit_VelAcc;
%current_err = cosnlin_err(vect_VelAcc);
err_VelAcc = cosnlin_err(vect_VelAcc);

figure(FigureIndex);axes('position',[0.26 0.31 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); colorbar;   
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');
axis off;
title({[ 'Err: ' num2str(err_VelAcc, '%0.2f') ]},  'FontSize', 10);% axis off;

% Do chi-square goodness of fit test
clear DataDiff;DataDiff=spacetime_data_n-sqrt(spacefit_VelAcc);%DataDiff=sqrt(spacetime_data)-sqrt(spacefit_VelAcc);
dof_VelAcc=(size(spacetime_data,1)-1)*(size(spacetime_data,2)-1)-length(sd_zero)-length(vect_VelAcc);%degree of freedom
chi2_VelAcc=sum(nansum( (DataDiff./Data_std).^2));
RChi2_VelAcc = chi2_VelAcc/dof_VelAcc;%Reduced Chi-Squared value
chi2P_VelAcc=1-chi2cdf(chi2_VelAcc,dof_VelAcc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_use=3;
% vect_VelAccPos=[47.7 105.2 0.2 171.9*pi/180 1.239 0.264 0.108 0 169.6*pi/180 0.986 145.2*pi/180 0.089];
% spacefit_VelAccPos= funccosnlin(vect_VelAccPos,xtdata);
[spacefit_VelAccPos,vect_VelAccPos,r_squared_VelAccPos,CI_VelAccPos] = MSF_VelAccPos_fit(spacetime_data, vect_VelAcc,Step*0.001,allow_negative);
% clear model_VelAccPos; model_VelAccPos=zeros(size(PSTH_tempPlane));
% model_VelAccPos(:,1:25) = NaN; model_VelAccPos(:,26:151)=spacefit_VelAccPos(1:8,:);

figure(FigureIndex);axes('position',[0.05 0.09 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,spacefit_VelAccPos',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)'); 
ylabel('Time (sec)');title('model: Vel + Acc + Pos');

error_surf = spacetime_data - spacefit_VelAccPos;
%current_err = cosnlin_err(vect_VelAccPos);
err_VelAccPos = cosnlin_err(vect_VelAccPos);
figure(FigureIndex);axes('position',[0.26 0.09 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); colorbar;   
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:45:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
axis off; %xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');
title({[ 'Err: ' num2str(err_VelAccPos, '%0.2f') ]},  'FontSize', 10);% axis off;

% Do chi-square goodness of fit test
clear DataDiff;DataDiff=sqrt(spacetime_data)-sqrt(spacefit_VelAccPos);
dof_VelAccPos=(size(spacetime_data,1)-1)*(size(spacetime_data,2)-1)-length(sd_zero)-length(vect_VelAccPos);%degree of freedom
chi2_VelAccPos=sum(nansum((DataDiff./Data_std).^2));%chi2_VelAccPos=sum(nansum((DataDiff./Data_std).^2));
RChi2_VelAccPos = chi2_VelAccPos/dof_VelAccPos;%Reduced Chi-Squared value
% chi2P_VelAccPos=1-chi2cdf(chi2_VelAccPos,dof_VelAccPos);
chi2P_VelAccPos=1-chi2cdf(err_VelAccPos,dof_VelAccPos);
chi2P_VelAcc=1-chi2cdf(err_VelAcc,dof_VelAcc);
chi2P_Vel=1-chi2cdf(err_Vel,dof_Vel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[];
out = 'CurveFitting:    R0 |  Amp  |  n  | muAzi |  muT  | sigmaT |  DC2  | wVel  |ThetaAcc|  wAcc | ThetaPos | wPos'; 
out = strvcat(out, sprintf('------------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=vect_Vel;
OutputValue(4)=OutputValue(4)*180/pi;
out=strvcat(out, sprintf('Vel Model:     %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f', OutputValue));

clear OutputValue;OutputValue=vect_VelAcc;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('Vel+Acc Model: %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f  ', OutputValue));

clear OutputValue; OutputValue=vect_VelAccPos;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
OutputValue(11)=OutputValue(11)*180/pi;
out=strvcat(out, sprintf('Vel+Acc+Pos  : %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f |%6.3f | %6.1f   |%6.3f', OutputValue)); 
figure(FigureIndex);axes('position',[0.26 0.74 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

out=[];
out = 'Goodness of fitting: RChi2(Vel) | ChiP(Vel) | RChi2(VelAcc) | ChiP(VelAcc) | RChi2(VelAccPos)  | ChiP(VelAccPos) '; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[RChi2_Vel chi2P_Vel RChi2_VelAcc chi2P_VelAcc RChi2_VelAccPos chi2P_VelAccPos];
out=strvcat(out, sprintf('                        %6.3f  |   %6.3f  |      %6.3f   |   %6.3f     |  %6.3f           |   %6.3f    ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.59 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

x_start = [0, 0];
x_stop =  [2, 2];
y_marker=[0,1.1*max(max(spacetime_data))];

xscale = [0.8 0.77 0.65 0.52 0.5 0.52 0.65 0.77];   
yscale = [0.35 0.5 0.65 0.5 0.35 0.2 0.05 0.2];

% x_time=[1:size(PSTH_tempPlane,2)]*0.02-0.5;
x_time=t;
for i=1:length(Azi_temp)
    i
    axes('position',[xscale(i) yscale(i) 0.16 0.10]);
    bar(x_time,PSTH_tempPlane(i,:));hold on;%bar(x_time, count_y{i,j,k});    hold on;
   
    plot(x_time,model_Vel(i,:),'r','LineWidth',2); hold on;   
    plot(x_time,model_VelAcc(i,:), 'g', 'LineWidth', 2);hold on;
    plot(x_time,model_VelAccPos(i,:),'c','LineWidth',2);hold on;
    
    plot( x_start, y_marker, 'k-');hold on;
    plot( x_stop,  y_marker, 'k-');hold on;
    xlim([-0.5,2.5]);    ylim([0,1.2*max(max(spacetime_data))]);
    
    set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
    set(gca, 'xtick',[-0.5:0.5:2.5]);
    set( gca, 'xticklabel', '-0.5|0|0.5|1|1.5|2|2.5');
    set( gca, 'yticklabel', ' ' );
    title({['azi=' num2str(Azi_temp(i)) '; Ele=' num2str(Ele_temp(i))]}, 'FontSize', 8);
end  














