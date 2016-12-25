function Pursuit_Visual_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,2);

temp_pursuit_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
temp_pursuit_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
temp_pursuit_amplitude = data.moog_params(ROT_AMPLITUDE,:,2);
temp_pursuit_duration = data.moog_params(ROT_DURATION,:,MOOG);
% temp_rot_type = data.moog_params(25,:,MOOG);

temp_spike_data = data.spike_data(SpikeChan,:);
temp_event_data = data.spike_data(2,:); %2

% remove null trials, bad trials, and trials outside Begtrial~Endtrial
data_duration = length(temp_spike_data)/length(temp_azimuth); %usually 5000
spike_data_good = temp_spike_data;

% Each column is a new trial
trials = 1:length(temp_azimuth);
spike_data = reshape(spike_data_good, data_duration, length(trials));
event_data = reshape(temp_event_data, data_duration, length(trials));


%Find trials with no visual pursuit
null_trials = logical(temp_pursuit_amplitude == -9999);

[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, size(data.event_data, 3), 4, 5, 500, -500, data.UseSyncPulses);

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_pursuit_azimuth);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);

pursuit_azimuth = temp_pursuit_azimuth(~null_trials & select_trials);
pursuit_elevation = temp_pursuit_elevation(~null_trials & select_trials);
pursuit_amplitude = temp_pursuit_amplitude(~null_trials & select_trials);
pursuit_duration = temp_pursuit_duration(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
% rot_type = temp_rot_type(~null_trials & select_trials);

% List of unique conditions
unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

unique_pursuit_azimuth = munique(pursuit_azimuth');
unique_pursuit_elevation = munique(pursuit_elevation');
unique_pursuit_amplitude = munique(pursuit_amplitude');
unique_pursuit_duration = munique(pursuit_duration');
% unique_rot_type = munique(rot_type');

% spont_trials = find(pursuit_amplitude ~= -9999);
% pursuit_heading_trials = find(rot_type == 8 & stim_type ==2);
% pursuit_only_trials = find(rot_type == 8 & stim_type ==0);
% fix_heading_trials = find(rot_type == -1);

eye = 1; % 1 left eye; 2 right eye;

if eye == 1
    eye_x=1;
    eye_y=2;
else
    eye_x=3;
    eye_y=4;
end

select_left_trial =  logical(pursuit_amplitude > 0); %Leftward pursuit
left_trials = find(select_left_trial ==1);
select_right_trial =  logical(pursuit_amplitude < 0); %Rightward pursuit
right_trials = find(select_right_trial ==1);

temp_eyex(:,:) = data.eye_data(eye_x,:,:);
temp_eyey(:,:) = data.eye_data(eye_y,:,:);

eyex = temp_eyex(:,~null_trials & select_trials);
eyey = temp_eyey(:,~null_trials & select_trials);

for i= 1:length(pursuit_amplitude)
    pursuit_start(i) = find(data.event_data(1,:,i) == 4 );
    fix_on(i) = find(data.event_data(1,:,i) == 3 );
end

temp_event_data = data.spike_data(2,:);
event_data_all = reshape(temp_event_data, 5000, length(temp_amplitude));
event_data = event_data_all(:,~null_trials & select_trials);

for i = 1:length(pursuit_amplitude)
    for j = 1:5000
    if (event_data(j,i) == 1)
        stim_on_pulse(i) = floor(j/5); %using the sync pulse
        break
    end
    end
end

for i = 1:length(pursuit_amplitude)
    if (stim_on_pulse(i) == 0)
        stim_on_pulse(i) = 213; %Average stim_on_pulse value
    end
end

for i = 1 : length(left_trials)
    stim_on_pulse_l(i) = stim_on_pulse(left_trials(i));
    eyex_l(:,i) = eyex((stim_on_pulse_l(i)-30:stim_on_pulse_l(i)+330), left_trials(i)); %starting from fix pt on...
    eyey_l(:,i) = eyey((stim_on_pulse_l(i)-30:stim_on_pulse_l(i)+330), left_trials(i));
end

for i = 1 : length(right_trials)
    stim_on_pulse_r(i) = stim_on_pulse(right_trials(i));
    eyex_r(:,i) = eyex((stim_on_pulse_r(i)-30:stim_on_pulse_r(i)+330), right_trials(i));
    eyey_r(:,i) = eyey((stim_on_pulse_r(i)-30:stim_on_pulse_r(i)+330), right_trials(i));
end

m = 10; %bin size
n = 3; %step size

% Use Gaussian smoothing on position before calculating velocity (like Yong does)
for i=1:length(left_trials)
    gauss2 = normpdf(1:100,50,10); % smooth data at a SD of 50ms (10 points. Probably too much)
    eye_temp = conv(eyex_l(:,i), gauss2);
    eyex_l_smooth(:,i) = eye_temp(50:end-50);
    eyex_vel_l(:,i) = (diff(eyex_l_smooth(:,i)))*200;

%     for j = 1:35
%         eyex_vel_l_diff(j,i) = (mean(eyex_l_smooth(1+(m*(j-1)+1:(m*j+n)),i)) - mean(eyex_l_smooth(1+(m*(j-1)+1:(m*j)),i)))/(.005*n/2);
%     end
%     eyex_vel_l_fderiv(:,i) =fderiv(eyex_l(1:40,i), 5, 200);
end

for i=1:length(right_trials)
    gauss2 = normpdf(1:100,50,10); % smooth data at a SD of ___
    eye_temp = conv(eyex_r(:,i), gauss2);
    eyex_r_smooth(:,i) = eye_temp(50:end-50);
    eyex_vel_r(:,i) = (diff(eyex_r_smooth(:,i)))*200;
%     for j = 1:35
%         eyex_vel_r_diff(j,i) = (mean(eyex_r_smooth(1+(m*(j-1)+1:(m*j+n)),i)) - mean(eyex_r_smooth(1+(m*(j-1)+1:(m*j)),i)))/(.005*n/2);
%         %     eyex_vel_r_diff(j,i) = (eyex_r(stim_on_pulse_r(i)+20*(j-1)+1,i) - eyex_r(stim_on_pulse_r(i)+20*j,i))/(.005*20);
%     end
%     eyex_vel_r_fderiv(:,i) =fderiv(eyex_r(1:40,i), 5, 200);
end

peak_velocity = abs(unique_pursuit_amplitude(1))*2/(.75+1.5);
position_l = (-peak_velocity*(0:750)/1000)+(10-3.333);
position_r = (peak_velocity*(0:750)/1000)+(-10+3.333);

mean_eyex_l = mean(eyex_l,2);
mean_eyex_r = mean(eyex_r,2);
mean_eyey_l = mean(eyey_l,2);
mean_eyey_r = mean(eyey_r,2);
mean_eyex_vel_r = mean(eyex_vel_r,2);
mean_eyex_vel_l = mean(eyex_vel_l,2);

error_l = std(eyex_vel_l,0,2);
error_r = std(eyex_vel_r,0,2);

figure(2)
for i = 1:length(left_trials)
    subplot(2,2,1), plot(5*(-30:330),eyex_l(:,i));
    hold on
    subplot(2,2,3), plot(5*(-29:330),eyex_vel_l(:,i));
    hold on
end
subplot(2,2,1), plot(5*(-30:330),mean_eyex_l,'r','LineWidth',2);
hold on;
subplot(2,2,1), plot((375+65:375+750+65),position_l+1.75,'k','LineWidth',1);
hold on
subplot(2,2,1), plot((375+65:375+750+65),position_l-1.75,'k','LineWidth',1);
title('Position - Left Pursuit');
ylabel('X position');
set(gca,'xlim',[0 1500]);
subplot(2,2,3), plot(5*(-29:330),mean_eyex_vel_l,'r','Linewidth',2);
hold on
subplot(2,2,3), plot(5*(-29:330),mean_eyex_vel_l+error_l, 'r-',5*(-29:330),mean_eyex_vel_l-error_l,'r-');
subplot(2,2,3), line([375+65;375+750+65],[-peak_velocity, -peak_velocity],'Color','g', 'Linewidth',2);
subplot(2,2,3), line([0+65;375+65],[0, -peak_velocity],'Color','g', 'Linewidth',2);
subplot(2,2,3), line([375+750+65;1500+65],[-peak_velocity, 0],'Color','g', 'Linewidth',2);
title('Velocity - Left Pursuit');
xlabel('Time (ms)');
set(gca,'ylim',[-40 20]);
set(gca,'xlim',[0 1500]);

for i = 1:length(right_trials)
    subplot(2,2,2), plot(5*(-30:330),eyex_r(:,i));
    hold on
    subplot(2,2,4), plot(5*(-29:330),eyex_vel_r(:,i));
    hold on
end
subplot(2,2,2), plot(5*(-30:330),mean_eyex_r, 'r','LineWidth',2);
hold on;
subplot(2,2,2), plot((375+65:375+750+65),position_r+1.75,'k','LineWidth',1);
hold on
subplot(2,2,2), plot((375+65:375+750+65),position_r-1.75,'k','LineWidth',1);
title('Position - Right Pursuit');
ylabel('X position');
set(gca,'xlim',[0 1500]);
subplot(2,2,4), plot(5*(-29:330),mean_eyex_vel_r,'r','Linewidth',2);
hold on
subplot(2,2,4), plot(5*(-29:330),mean_eyex_vel_r+error_r,'r-',5*(-29:330),mean_eyex_vel_r-error_r,'r-');
subplot(2,2,4), line([375+65;375+750+65],[peak_velocity, peak_velocity],'Color','g', 'Linewidth',2);
subplot(2,2,4), line([0+65;375+65],[0, peak_velocity],'Color','g', 'Linewidth',2);
subplot(2,2,4), line([375+750+65;1500+65],[peak_velocity, 0],'Color','g', 'Linewidth',2);
title('Velocity - Right Pursuit');
xlabel('Time (ms)');
set(gca,'ylim',[-20 40]);
set(gca,'xlim',[0 1500]);

axes('position',[0,0,1,1],'visible','off');
text(0.88,0.98,FILE);
hold on



% m_path = 'Z:\Users\Adhira\';
% saveas(gcf,[m_path, FILE(1:end-4)],'fig');
% saveas(gcf,[m_path, FILE(1:end-4)]);

return;