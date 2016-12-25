function CuveFitting_Hor(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, OutputPath)

warning off MATLAB:divideByZero;
warning off MATLAB:singularMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc
load time_vel_acc.mat;%load('C:\MATLAB6p5\work\time_vel_acc.mat');
load m14c136r1.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); time = clock; disp(['Start Time = ' num2str([time(4) time(5) time(6)])]);

Path_Defs;
ProtocolDefs;  %contains protocol specific keywords - 1/4/01 BJP

plot_psth = 1;
plot_spatial = 1;

dftcutoff = 2.5;    % Select the DFTR cut-off 
polswitch = 1;      % Select 1 for switching signs on the component phases
use_max_phase = 0;

load('Z:\Users\Aihua\CurveFittingAnalysis\outputs\time_vel_acc.mat');

% load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\bpr_lookup_table.mat');
% load('Z:\Data\Tempo\Batch Files\Suhrud\bpr_lookup_table.mat');

if Protocol == 100 || Protocol == 104
	%get the column of values for azimuth and elevation and stim_type
	temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
	temp_elevation = data.moog_params(ELEVATION,:,MOOG);
	temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
	temp_spike_data = data.spike_data(SpikeChan,:);
	% fixation data added by CRF, 3/2007
	temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
	temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
	temp_fix_x(isnan(temp_fix_x)) = 0;
	temp_fix_y(isnan(temp_fix_y)) = 0;
elseif Protocol == 112
	%get the column of values for azimuth and elevation and stim_type
	temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
	temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
	temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
	temp_spike_data = data.spike_data(SpikeChan,:);
	% fixation data added by CRF, 3/2007
	temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
	temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
	temp_fix_x(isnan(temp_fix_x)) = 0;
	temp_fix_y(isnan(temp_fix_y)) = 0;    
end
    
    
%mean firing rate of each trial depending on the start and stop offsets
temp_spike_rates = data.spike_rates(SpikeChan, :);

%--------------------------------------------------------------
% Some cells have more than the usual 26 directions, so we need to remove them
extra_directions =[];
if length(munique(temp_azimuth')) > 9  % it's 9 because the null trials (-9999) are included
    unique_azimuth_withnull = [0;45;90;135;180;225;270;315;-9999];
    unique_elevation_withnull = [-90;-45;0;45;90;-9999];
    extra_azis = ones(1,length(temp_azimuth));
    extra_eles = ones(1,length(temp_elevation));
    for n = 1:length(unique_azimuth_withnull)
        extra_azis(find(temp_azimuth == unique_azimuth_withnull(n))) = 0;
    end
    for m = 1:length(unique_elevation_withnull)
        extra_eles(find(temp_elevation == unique_elevation_withnull(m))) = 0;
    end
    extra_directions = (extra_azis | extra_eles);
    temp_elevation(extra_directions) = [];
    temp_stim_type(extra_directions) = [];
    temp_fix_x(extra_directions) = [];
    temp_fix_y(extra_directions) = [];
    temp_spike_rates(extra_directions) = [];    
    temp_azimuth(extra_directions) = [];
    
    EndTrial = EndTrial - sum(extra_directions(1:EndTrial-(BegTrial-1)));
end
%--------------------------------------------------------------

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth); % a vector of trial indices
%at any given point, cell response cannot be greater than 1
abnormal = find(temp_spike_data > 1);
temp_spike_data(1,abnormal) = 1;
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response

if ( bad_trials ~= NaN)
    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else
    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
fix_x = temp_fix_x(~null_trials & select_trials);
fix_y = temp_fix_y(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

fix_x_withnull = temp_fix_x(select_trials);
fix_y_withnull = temp_fix_y(select_trials);
spike_rates_withnull = temp_spike_rates(select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_fix_x    =  munique(fix_x');
unique_fix_y    =  munique(fix_y');

temp_condition_num = temp_stim_type;
condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_psth = 1;
timebin = 20; time = time20; vel = vel20; acc = acc20;

% sample frequency depends on test duration
% frequency=length(temp_spike_data)/length(select_trials); % actually its 5000 now for 5 sec 
frequency = 5000;

% length of x-axis
x_length = round(frequency/timebin);
% x-axis for plot PSTH
x_time=1:x_length;
x_time_bincenter = timebin/2000;            % using bin-centers now: start at time = 0, plus a half-bin...
for nn = 1 : round(5/(timebin/1000)) - 1    % then fill it out with a full 2-seconds (i.e., 40 bins with timebin = 50 ms)
    x_time_bincenter(end+1) = x_time_bincenter(end) + (timebin/1000);
end

% remove 'extra direction' trials from spikedata
discard_trials = find(extra_directions);
for i = 1 : length(discard_trials)
    temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) : discard_trials(i)*frequency ) = 9999;
end
temp_spike_data = temp_spike_data( temp_spike_data~=9999 );

% remove trials outside Begtrial~Endtrial
discard_trials = [find(trials <BegTrial | trials >EndTrial)];
for i = 1 : length(discard_trials)
    temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) : discard_trials(i)*frequency ) = 9999;
end
spike_data_withnull = temp_spike_data( temp_spike_data~=9999 );

% remove null trials and trials outside Begtrial~Endtrial
discard_trials = [find(null_trials==1 | trials <BegTrial | trials >EndTrial)];
for i = 1 : length(discard_trials)
    temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) :  discard_trials(i)*frequency ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );

% count spikes from raster data (spike_data)
max_count = 1;
time_step = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_stim = 1;  % vestibular=1, visual=2
actual_stim = current_stim;
unique_condition_num = 1;
% % For visual-only blocks, the index 'k' (below) will only reach 1, and no data
% % will be created at index k = 2.  Therefore must change current_stim to 1. -CRF
% if unique_condition_num == 2 & current_stim == 2
%     current_stim = 1;
%     actual_stim = 2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(unique_condition_num)  % it's the stim_type
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)         % SELECT now restricted to only 0 deg fixation trials -- CRF 3/2007
            select = logical(azimuth==unique_azimuth(i) & elevation==unique_elevation(j) & condition_num==unique_condition_num(k) & fix_x==0 & fix_y==0 );             
            if sum(select) > 0
                resp{k}(j,i) = mean(spike_rates(select));
                act_found = find( select==1 );
                % count spikes per timebin on every same condition trials
                clear temp_count dummy_count;
                for repeat=1:length(act_found) 
                    for n=1:(x_length)
                        temp_count(repeat,n)=sum( spike_data( 1, frequency*(act_found(repeat)-1) + time_step : round(frequency*(act_found(repeat)-1)+n*timebin) ) );
                        dummy_count{repeat}(n) = sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):round(frequency*(act_found(repeat)-1)+n*timebin)));
                        
                        if timebin == 16.6667
                            if floor((n+1)/3) == (n+1)/3
                                time_step = floor(time_step + timebin);
                            else                               % kluge for 16.6667 ms timebin (just rounding was missing bins in spike_data)
                                time_step = round(time_step + timebin);
                            end
                        else
                            time_step = time_step + timebin;   % this may need to be modified for timebins other than 50
                        end
                    end
                    time_step=1;
                    
                    if k == current_stim
                        count_condition{i,j,repeat} = dummy_count{repeat};
                    end
                end

                count_y_trial{i,j,k}(:,:) = temp_count;  % each trial's PSTH in vestibular condition
                if k == 1
                    count_y_trial1{i,j}(:,:) = temp_count;
                end

                dim=size(temp_count);
                if dim(1) > 1;
                    count_y{i,j,k} = mean(temp_count);
                    std_y{i,j,k}=std(temp_count);
                else
                    count_y{i,j,k}= temp_count;     % for only one repetition cases
                    std_y{i,j,k}=NaN;
                end
            else
                resp{k}(j,i) = 0;
                count_y{i,j,k}=0;
                std_y{i,j,k}=NaN;
            end
            % normalize count_y
            if max(count_y{i,j,k})~=0;
                count_y_norm{i,j,k}=count_y{i,j,k} / max(count_y{i,j,k});
            else
                count_y_norm{i,j,k}=0;
            end
        end
    end
    % now find the peak
    [row_max, col_max] = find( resp{k}(:,:)==max(max(resp{k}(:,:))) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k}=row_max(1);
    col_m{k}=col_max(1);
    if max(count_y{col_max(1), row_max(1), k})~=0;
        count_y_max{k} = count_y{col_max(1), row_max(1), k} / max(count_y{col_max(1), row_max(1), k});
    else
        count_y_max{k} =0;
    end
    % find the largest y to set scale later
    if max(count_y{col_max(1), row_max(1), k}) > max_count
        max_count = max(count_y{col_max(1), row_max(1), k});
    end
end

% % Do the same, but for null trials (for new significance test on DFTR below) -- CRF 10/2007
% select = null_trials(select_trials) & temp_fix_x(select_trials)==0 & temp_fix_y(select_trials)==0;
% if sum(select) > 0
%     resp_null = mean(spike_rates_withnull(select));
%     act_found = find( select==1 );
%     % count spikes per timebin on every same condition trials
%     clear temp_count dummy_count;
%     for repeat=1:length(act_found) 
%         for n=1:(x_length)
%             temp_count(repeat,n)=sum( spike_data_withnull( 1, frequency*(act_found(repeat)-1) + time_step : round(frequency*(act_found(repeat)-1)+n*timebin) ) );
%             dummy_count{repeat}(n) = sum(spike_data_withnull(1,(frequency*(act_found(repeat)-1)+time_step):round(frequency*(act_found(repeat)-1)+n*timebin)));
%             if timebin == 16.6667
%                 if floor((n+1)/3) == (n+1)/3
%                     time_step = floor(time_step + timebin);
%                 else                               % kluge for 16.6667 ms timebin (just rounding was missing bins in spike_data)
%                     time_step = round(time_step + timebin);
%                 end
%             else
%                 time_step = time_step + timebin;   % this may need to be modified for timebins other than 50
%             end
%         end
%         time_step=1;
%     end
%     count_y_trial_null = temp_count;
%     dim=size(temp_count);
%     if dim(1) > 1;
%         count_y_null = mean(temp_count);
%     else
%         count_y_null = temp_count;     % for only one repetition cases
%     end
% else
%     resp_null=0;
%     count_y_null=0;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xy: Horizontal plane
Azi_hor=[0:45:315];
Ele_hor=0*ones(1,8);

%xz: Frontal plane
Azi_Front=[0 0 0 0 0 180 180 180];
Ele_Front=[-90:45:90 -45:45:45];

%yz: Median plane (or mid-sagital plane)
Azi_Med=[0 0 90 90 90 -90 -90 -90];
Ele_Med=[-90 90 -45:45:45 -45:45:45];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the space time data for the horizontal plane







%%% --------------------TAKE PSTHs in one elevation, plot as direction-time maps-------------------------
%%% Smooth the data using Box-Car filter (same as a running
%%% average), remove the last 60 ms which contains no data and then
%%% save into Elevation1 (ele = -90 degrees), Elevation2 (-45), Elevation3 (0),
%%% Elevation4 (+45) and Elevation5 (+90)

Ele3(:,1) = count_y(1,3);
Ele3(:,2) = count_y(2,3);
Ele3(:,3) = count_y(3,3);
Ele3(:,4) = count_y(4,3);
Ele3(:,5) = count_y(5,3);
Ele3(:,6) = count_y(6,3);
Ele3(:,7) = count_y(7,3);
Ele3(:,8) = count_y(8,3);
Elevation3 = flipud(squeeze(Ele3(:,:)));

%%% Second term in boxcarfilter is number of bins to use for averaging
A3 = zeros(8, 250);
for i = 1:length(Elevation3)
    A3(i,:) = Elevation3{i};
end

for i =1:8
    A3Smooth2(i,:) = BoxcarFilter(A3(i,:), 2);
    A3Smooth4(i,:) = BoxcarFilter(A3(i,:), 4);
    %         A3Smooth4(i,:) = A3(i,:); 
end
% A3junk = A3Smooth4(:, 1:190);
% A3Smooth4 = A3junk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------------DIRECTION TIME FITTING----------------------------------
% global xdata tdata spacetime_data;
xt_data = A3Smooth4(:,50:175);%xt_data = A3Smooth4(:,50:150);
%xt_data = [A2Smooth4(:,50:150); A3Smooth4(:,50:150); A4Smooth4(:,50:150)];  %% Data to use for curvefitting -- (bins 50:150 = 2s stim)
xscaling = 1:8;
for i = 1:size(xt_data,2)%101
    xs(i,:) = xscaling;
end
yscaling = zeros([1 size(xt_data,2)]);%yscaling = zeros([1 101]);
for i = 1:length(yscaling)
    yscaling(i) = i*20/1000 - 0.02;
end
for i = 1:8
    ys(i,:) = yscaling;
end
ys = ys';

%%% This is the main loop for curvefitting. It uses lsqcurvefit (can be
%%% changed to fmincon or any other algorithm). The loop calls two
%%% functions viz. funccosnlin and cosnlin_err. Two models are
%%% implemented: first one contains only the velocity term whereas the
%%% second one also has the acceleration term.
xt_data=xt_data;%/timebin*1000;
spacetime_data = xt_data;%spacetime_data = xt_data((ele_2_use-1)*8+1: 8*ele_2_use, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the original space-time figure; 
FigureIndex=2; 
figure(FigureIndex);set(FigureIndex,'Position', [50,100 1200,800], 'Name', 'CurveFitting');orient landscape; 
% text(-0.1,1.06,[FILE '       Elevevation = 0 deg']); axis off;    
clear XAzi Ytime;[XAzi,Ytime] = meshgrid([0:45:315],[0:0.02:max(yscaling)]);
figure(FigureIndex);axes('position',[0.05 0.77 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,spacetime_data');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(yscaling)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');title('raw data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model1: Velocity only
global model_use
model_use=1;
[spacefit_Vel,vect_Vel,r_squared_Vel,CI_Vel] = MSF_Vel_fit(spacetime_data);
model_Vel = zeros(250, 8);model_Vel(1:49, :) = NaN; 
model_Vel(50:175, :) = spacefit_Vel';model_Vel(176:250, :) = NaN;

figure(FigureIndex);axes('position',[0.05 0.53 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_Vel',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(yscaling)]);
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
set(gca, 'ytick',[0:0.5:max(yscaling)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');
axis off;
title({[ 'Err: ' num2str(err_Vel, '%0.2f') ]},  'FontSize', 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model 2: Velocity + Acceleration
model_use=2;
[spacefit_VelAcc,vect_VelAcc,r_squared_VelAcc,CI_VelAcc] = MSF_VelAcc_fit(spacetime_data, vect_Vel);
model_VelAcc = zeros(250, 8);model_VelAcc(1:49, :) = nan; 
model_VelAcc(50:175, :) = spacefit_VelAcc';model_VelAcc(176:250, :) = nan;

% Velocity + Acceleration model:
figure(FigureIndex);axes('position',[0.05 0.31 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_VelAcc');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(yscaling)]);
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
set(gca, 'ytick',[0:0.5:max(yscaling)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');
axis off;
title({[ 'Err1: ' num2str(err_VelAcc, '%0.2f') ]},  'FontSize', 10);% axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model 3: Velocity + Acceleration + Position
model_use=3;
[spacefit_VelAccPos,vect_VelAccPos,r_squared_VelAccPos,CI_VelAccPos] = MSF_VelAccPos_fit(spacetime_data, vect_VelAcc);
model_VelAccPos = zeros(250, 8);model_VelAccPos(1:49, :) = nan; 
model_VelAccPos(50:175, :) = spacefit_VelAccPos';model_VelAccPos(176:250, :) = nan;

figure(FigureIndex);axes('position',[0.05 0.09 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,spacefit_VelAccPos',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(yscaling)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)'); 
ylabel('Time (sec)');title('model: Vel + Acc + Pos');

error_surf = spacetime_data - spacefit_VelAccPos;
%current_err = cosnlin_err(vect_VelAccPos);
err_VelAccPos = cosnlin_err(vect_VelAccPos);
figure(FigureIndex);axes('position',[0.26 0.09 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); colorbar;   
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:45:315]);set(gca, 'xticklabel','0|45|90|135|180|215|270|315');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(yscaling)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
axis off; %xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');
title({[ 'Err: ' num2str(err_VelAccPos, '%0.2f') ]},  'FontSize', 10);% axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[];
out = 'CurveFitting:      R0 |  Amp  |   n   | muAzi | muTime | sigmaTime |  DC2  | wVel  |ThetaAcc|  wAcc | ThetaPos'; 
out = strvcat(out, sprintf('-----------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=vect_Vel;
OutputValue(4)=OutputValue(4)*180/pi;
out=strvcat(out, sprintf('Vel Model:     %6.3f |%6.3f |%6.3f |%6.1f | %6.3f |  %6.3f   |%6.3f', OutputValue));

clear OutputValue;OutputValue=vect_VelAcc;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('Vel+Acc Model: %6.3f |%6.3f |%6.3f |%6.1f | %6.3f |  %6.3f   |%6.3f |%6.3f | %6.1f  ', OutputValue));

clear OutputValue; OutputValue=vect_VelAccPos;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
OutputValue(11)=OutputValue(11)*180/pi;
out=strvcat(out, sprintf('Vel+Acc+Pos  : %6.3f |%6.3f |%6.3f |%6.1f | %6.3f |  %6.3f   |%6.3f |%6.3f | %6.1f |%6.3f |%6.1f', OutputValue)); 
figure(FigureIndex);axes('position',[0.26 0.74 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing regression to get R^2 and p-values for F-test
Data_Raw=reshape(spacetime_data,size(spacetime_data,1)*size(spacetime_data,2),1);
Data_Vel=reshape(spacefit_Vel,size(spacefit_Vel,1)*size(spacefit_Vel,2),1);
Data_VelAcc=reshape(spacefit_VelAcc,size(spacefit_VelAcc,1)*size(spacefit_VelAcc,2),1);
Data_VelAccPos=reshape(spacefit_VelAccPos,size(spacefit_VelAccPos,1)*size(spacefit_VelAccPos,2),1);

clear tempstd; %tempstd=std_y(3,50:175);
tempstd(1,:) = std_y(1,3);
tempstd(2,:) = std_y(2,3);
tempstd(3,:) = std_y(3,3);
tempstd(4,:) = std_y(4,3);
tempstd(5,:) = std_y(5,3);
tempstd(6,:) = std_y(6,3);
tempstd(7,:) = std_y(7,3);
tempstd(8,:) = std_y(8,3);
std_hor0 = cell2mat(tempstd);
std_hor=std_hor0(:,50:175);
Data_std=reshape(std_hor,size(std_hor,1)*size(std_hor,2),1);
clear sd_zero_count; sd_zero_count=find(Data_std==0);Data_std(sd_zero_count)=NaN;

clear X1;X1 = [ones(size(Data_Vel,1),1) Data_Vel];%X1 = [ones(808,1) Data_Vel];% y_fit = [ones(length(y_fit),1) y_fit];
[b_Vel,bint,r_Vel,rint,stats_Vel] = regress(Data_Raw,X1,0.05);

% Do chi-square goodness of fit test
df_mean=(length(Data_Raw)-length(sd_zero_count))-length(vect_Vel);
chi2_Vel = nansum( (sqrt(Data_Raw)-sqrt(Data_Vel)).^2 ./ Data_std.^2 );
chiP_Vel = 1 - chi2cdf(chi2_Vel, df_mean);

clear X1;X1 = [ones(size(Data_VelAcc,1),1) Data_VelAcc];
[b_VelAcc,bint,r_VelAcc,rint,stats_VelAcc] = regress(Data_Raw,X1,0.05);
% Do chi-square goodness of fit test
df_mean=(length(Data_Raw)-length(sd_zero_count))-length(vect_VelAcc);
chi2_VelAcc = nansum( (sqrt(Data_Raw)-sqrt(Data_VelAcc)).^2 ./ Data_std.^2 );
chiP_VelAcc = 1 - chi2cdf(chi2_VelAcc, df_mean);

clear X1;X1 = [ones(size(Data_VelAccPos,1),1) Data_VelAccPos];
[b_VelAccPos,bint,r_VelAccPos,rint,stats_VelAccPos] = regress(Data_Raw,X1,0.05);
df_mean=(length(Data_Raw)-length(sd_zero_count))-length(vect_VelAccPos);
chi2_VelAccPos = nansum( (sqrt(Data_Raw)-sqrt(Data_VelAccPos)).^2 ./ Data_std.^2 );
chiP_VelAccPos = 1 - chi2cdf(chi2_VelAccPos, df_mean);

Ftest_1vs2=[(err_Vel-err_VelAcc)/(length(vect_VelAcc)-length(vect_Vel))]/[err_VelAcc/(size(spacetime_data,1)*size(spacetime_data,2)-length(vect_Vel))];
p_1vs2=1-fcdf(Ftest_1vs2,length(vect_VelAcc)-length(vect_Vel),size(spacetime_data,1)*size(spacetime_data,2)-length(vect_Vel));

Ftest_2vs3=[(err_VelAcc-err_VelAccPos)/(length(vect_VelAccPos)-length(vect_VelAcc))]/[err_VelAccPos/(size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc))];
p_2vs3=1-fcdf(Ftest_2vs3,length(vect_VelAccPos)-length(vect_VelAcc),size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc));

out=[];
out = 'Goodness of fitting:  r2 (Vel) | r2 (VelAcc)| r2 (VelAccPos) | Ftest(1vs2) | p(1vs2)  | Ftest(2vs3) | p(2vs3)'; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[r_squared_Vel r_squared_VelAcc r_squared_VelAccPos Ftest_1vs2 p_1vs2 Ftest_2vs3 p_2vs3];
out=strvcat(out, sprintf('                        %6.3f |    %6.3f  |      %6.3f    |   %6.3f   |  %6.3f  |   %6.3f    |  %6.3f   ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.65 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

out=[];
out = 'Goodness of fitting:  Chi2 (Vel) | ChiP (Vel)| Chi2 (VelAcc) | ChiP(VelAcc) | Chi2(VelAccPos)  | ChiP(VelAccPos) '; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[chi2_Vel chiP_Vel chi2_VelAcc chiP_VelAcc chi2_VelAccPos chiP_VelAccPos];
out=strvcat(out, sprintf('                        %6.3f |    %6.3f  |      %6.3f    |   %6.3f   |  %6.3f  |   %6.3f    ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.59 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

% find when in the 5-sec trial period that the spike data cuts off
sum_dat = zeros(1,x_length);
for j = 1:26
    sum_dat = sum_dat + count_y{plot_col(j), plot_row(j), current_stim};
end
% and accordingly set length of data to read from PSTH
for ii = 1:x_length
    if sum(sum_dat(ii:end)) == 0
        span = ii;
        break
    end
end

if span < .80*x_length
    span = round(.68*x_length); % usually this is how much data there is
else % this is one of katsu's cells
    span = round(.80*x_length);
end

current_elevation = 0;
ax = axes('position',[0,0,1,1],'visible','off');
% tx = text(0.01,0.98, [FILE ' ' 'Elevation = '  num2str(current_elevation)  ' ' 'F3=' num2str(ftest3, '%0.2f') ' ' 'p3=' num2str(p3, '%0.2f')]);
tx = text(0.01,0.98, [FILE ' ' 'Elevation = '  num2str(current_elevation)]);
set(tx,'fontweight','bold');

j = find(unique_elevation==current_elevation);
k = current_stim;
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
xscale = [0.8 0.77 0.65 0.52 0.5 0.52 0.65 0.77];   
yscale = [0.35 0.5 0.65 0.5 0.35 0.2 0.05 0.2];
for i=1:length(unique_azimuth)
    axes('position',[xscale(i) yscale(i) 0.16 0.10]);
    bar(x_time, count_y{i,j,k});    hold on;
    plot(x_time, model_Vel(:, i), 'r', 'LineWidth', 2);
    plot(x_time,model_VelAcc(:, i), 'g', 'LineWidth', 2);
    plot(x_time, model_VelAccPos(:,i),'c','LineWidth',2);
    plot( x_start, y_marker, 'k-');
    plot( x_stop,  y_marker, 'k-');
    xlim([0,span]);    ylim([0,max_count]);
    set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
    set(gca, 'xtick',[0:50:span]);
    set( gca, 'xticklabel', '-1|0|1|2|3');
    set( gca, 'yticklabel', ' ' );
    title({['azi = ' num2str(unique_azimuth(i))]}, 'FontSize', 8);
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the figure
OutputPath=['C:\Aihua\z_TempOutputs\'];
figure(FigureIndex); 
set(gcf, 'PaperOrientation', 'portrait');
saveas(gcf,[OutputPath FILE(1:end-4) '_CurveFitting.png'],'png');
close(FigureIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output some data for further analysis
CurveFitHor.spacetime_data=spacetime_data;
CurveFitHor.spacefit_Vel=spacefit_Vel;
CurveFitHor.spacefit_VelAcc=spacefit_VelAcc;
CurveFitHor.spacetime_VelAccPos=spacefit_VelAccPos;
CurveFitHor.vect_Vel=vect_Vel;
CurveFitHor.vect_VelAcc=vect_VelAcc;
CurveFitHor.vect_VelAccPos=vect_VelAccPos;
CurveFitHor.CI_Vel=CI_Vel;
CurveFitHor.CI_VelAcc=CI_VelAcc;
CurveFitHor.CI_VelAccPos=CI_VelAccPos;
CurveFitHor.stats_Vel=stats_Vel;
CurveFitHor.stats_VelAcc=stats_VelAcc;
CurveFitHor.stats_VelAccPos=stats_VelAccPos;
CurveFitHor.err_Vel=err_Vel;
CurveFitHor.err_VelAcc=err_VelAcc;
CurveFitHor.err_VelAccPos=err_VelAccPos;
CurveFitHor.r_squared_Vel=r_squared_Vel;
CurveFitHor.r_squared_VelAcc=r_squared_VelAcc;
CurveFitHor.r_squared_VelAccPos=r_squared_VelAccPos;
CurveFitHor.Ftest_1vs2=Ftest_1vs2;
CurveFitHor.p_1vs2=p_1vs2;
CurveFitHor.Ftest_2vs3=Ftest_2vs3;
CurveFitHor.p_2vs3=p_2vs3;

CurveFitHor.chi2_Vel=chi2_Vel;
CurveFitHor.chiP_Vel=chiP_Vel;
CurveFitHor.chi2_VelAcc=chi2_VelAcc;
CurveFitHor.chiP_VelAcc=chiP_VelAcc;
CurveFitHor.chi2_VelAccPos=chi2_VelAccPos;
CurveFitHor.chiP_VelAccPos=chiP_VelAccPos;

SaveFileName=[OutputPath FILE(1:end-4) '_CurveFitHor'];
save(SaveFileName,'CurveFitHor'); clear SaveFileName;


