% Accelerometer_study.m -- figures out whether accelerometer data is
% comparable to that of the sled

%SURESFIT1h.m curtesy Jimmy
%
%  fitting interface for single unit (su)/eyemov routines
%
% function that fits a single sinusoid to data
% to average response and stimulus cycles
% outputs: gain, phase, dc
% estimation of gain (expressed relative to ampl=peak stimulus velocity)
% phase (in deg) and dc (in deg/s) for each cycle


%--	AYANNA 09/08/2006 based on YONG, 07/12/04
%-----------------------------------------------------------------------------------------------------------------------
function Sinusoid_analysis(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_duration = data.moog_params(SIN_DURATION,:,MOOG);
temp_frequency = data.moog_params(SIN_FREQUENCY,:,MOOG);
temp_trans_amplitude = data.moog_params(SIN_TRANS_AMPLITUDE,:,MOOG);
temp_rot_amplitude = data.moog_params(SIN_ROT_AMPLITUDE,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_interocular_dist = data.moog_params(INTEROCULAR_DIST,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);   

%get indices of any NULL conditions (for measuring spontaneous activity)
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);
motion_coherence = temp_motion_coherence(~null_trials & select_trials);
interocular_dist = temp_interocular_dist(~null_trials & select_trials);
num_sigmas = temp_num_sigmas(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_motion_coherence = munique(motion_coherence');
unique_interocular_dist = munique(interocular_dist');
unique_num_sigmas = munique(num_sigmas');

condition_num = stim_type;
h_title{1}='Sinusoid';
unique_condition_num = munique(condition_num');

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

%%%%%%%%%%%%%%%%%%%%ACCELEROMETER DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data.eye_data(6,1:200,:);%this is to look at the accelerometer data
% (either channel 5 or 6, where channel 5 is for IA and 6 is for NO).

% % NEED THE CONVERSION FACTOR FOR DEGREES TO VOLTS!
% This is what loadtempodata does to eye position signal data.
%%%temp1(DA_H, :, :) = temp1(DA_H, :, :).*(one_time_params(X_DEG_FULL_SCALE)/one_time_params(AD_RANGE));
% Accelerometer data is erroneously being converted like eye data. We need to get rid of this.

% position2 = data.eye_data(6,:,1); 
% position1 = position2(601:5600); 
% trace = position1 .* 2.5/(mean(position2(1:600))); %2.5volts over the mean baseline position.
% figure
% hold on
% title('accelerometer')
% xlabel('time (d/a units)')
% ylabel('voltage (1V = 1G)')
% plot(trace)
% hold off

% NOaccel = data.eye_data(6,:,1);
% IAaccel = data.eye_data(5,:,1);

% Make a loop that runs through all the trials for analisys

spiketrain = data.spike_data(SpikeChan,:,trials); %vectors for the spikechan data for all the trials

%Sinusoid fitting of the spikes
signal = data.spike_data(SpikeChan,:,1);
time = length(signal); %doesn't matter what units.
frqf = 0.5; %in Hz

output = suresfit1h_noah(frqf,time,signal);

sines = out(1)*sin(frqf*2*pi*signal + (out(2)/180)*pi) + out(3);
figure
plot(sines)
xlabel = ('time (d/a units)');
title('sine func.m plot');

NOaccel = data.eye_data(6,:,1); %this is the data in degrees that I get out of the computer.
position2 = NOaccel .* (65536/100); %NOW the data is in raw d/a units... not degrees.
position1 = position2(601:5600); %this is the data in d/a units, minus baseline at the start and end.
trace = position1 .* 2.5/(mean(position2(1:600))); %this turns data into volts. 1V=1G for accelerometer.
figure
hold on
title('accelerometer IA')
xlabel('time (d/a units)')
ylabel('voltage (1V = 1G)')
plot(trace)
hold off
gval1 = max(trace(1:1500))-min(trace(1:1500));
gval2 = max(trace(1501:3500))-min(trace(1501:3500));
gval3 = max(trace(3501:5000))-min(trace(3501:5000));
gavg1 = (gval1 + gval2 + gval3)/3

IAaccel = data.eye_data(5,:,2); %this is the data in degrees that I get out of the computer.
position5 = NOaccel*(65536/100); %NOW the data is in raw d/a units... not degrees.
position4 = position5(601:5600); %this is the data in d/a units, minus baseline at the start and end.
trace2 = position4 .* 2.5/(mean(position5(1:600))); %this turns data into volts. 1V=1G for accelerometer.
figure
hold on
title('accelerometer NO')
xlabel('time (d/a units)')
ylabel('voltage (1V = 1G)')
plot(trace2)
hold off

gval4 = max(trace2(1:1500))-min(trace2(1:1500));
gval5 = max(trace2(1501:3500))-min(trace2(1501:3500));
gval6 = max(trace2(3501:5000))-min(trace2(3501:5000));
gavg2 = (gval4 + gval5 + gval6)/3


%sinusoid fitting of the actual accelerometer data.
signal = trace; %accelerometer trace.
time = length(trace); %doesn't matter what units.
frqf = 1; %in Hz

out = suresfit1h_noah(frqf,time,signal);

% % Sine fitting part = data.eye_data(1,:,2); %spike channel
% x = (length(data.eye_data(1,:,2))/0.005); %time
% frqf = 0.5; %sinusoid frequency
% out = suresfit1h_noah(frqf,x,y); %gain (expressed relative to ampl=peak stimulus velocity),phase (deg) and dc (deg/sec)
% results = out;