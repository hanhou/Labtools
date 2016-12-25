%----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Yong, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------

function Azimuth1D_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_spike_data = data.spike_data(SpikeChan,:);
temp_spike_rates = data.spike_rates(SpikeChan, :);  
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_trials ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
% stim_duration = length(temp_spike_data)/length(temp_azimuth);
% spike_data = data.spike_data(1, ((BegTrial-1)*stim_duration+1):EndTrial*stim_duration);
spike_rates= temp_spike_rates(~null_trials & select_trials);
num_sigmas = temp_num_sigmas(~null_trials & select_trials);
% notice that this bad_trials is the number without spon trials 

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');

if length(unique_num_sigmas)>1
    condition_num = num_sigmas;
else
    condition_num = stim_type;
end
unique_condition_num = munique(condition_num');

% add parameters here
% timebin for plot PSTH
timebin=50;   % in ms
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% find spontaneous trials which azimuth,elevation,stim_type=-9999
spon_found = find(null_trials==1);     

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(null_trials==1 | temp_spike_rates > 3000 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
for k=1: length(unique_condition_num)
    for i=1: length(unique_azimuth)
        select = logical( (azimuth==unique_azimuth(i)) & (condition_num==unique_condition_num(k)) );            
        resp{k}(i) = mean(spike_rates(select));
        act_found = find( select==1 );
        % count spikes per timebin on every same condition trials
        for repeat=1:length(act_found) 
            for n=1:(x_length)
                temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                time_step=time_step+timebin;
            end
            time_step=1;                    
        end
        count_y_trial{k,i}(:,:) = temp_count;  % each trial's PSTH 
        count_y{i,k} = mean(temp_count);
        temp_max_count(i,k) = max(count_y{i,k});
    end      
end
max_count = max(max(temp_max_count));

% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count/(timebin/1000)];  % in Hz
% define figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', '1D Azimuth Tuning');
orient landscape;
title(FILE);
axis off;

xoffset=0;
yoffset=0;

% now plot
for k=1: length(unique_condition_num) 
    
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.42;
        xoffset = 0;
    end
    % output some text 
    axes('position',[0 0 1 0.9]); 
    xlim([-50,50]);
    ylim([-50,50]);
    text(-30+xoffset*100,52+yoffset*110, num2str(unique_condition_num(k)) );
%    text(-47,-40, 'Azim: 0       45       90        135        180        225        270       315');  
    text(25,-40, 'Translation');
    axis off;
    hold on;
    
    location = [6 3 2 1 4 7 8 9];
    plotAzimuths = [0:45:355];
    
    
    for j=1:length(plotAzimuths)        
        
        i = find(unique_azimuth == plotAzimuths(j));
        
        % axes('position',[0.05*i+0.01+xoffset 0.8+yoffset 0.045 0.045]);                             
        subplot_tight(3,3,location(j),[0.05 0.05]);  % HH 20130909
        
        set(bar( x_time,count_y{i,k}(1,:)/(timebin/1000),1),'Edgecolor','k','facecolor','k');   
        
        % xlabel(num2str(unique_azimuth(i)));
        hold on;
        plot( x_start, y_marker, 'r-');
        plot( x_stop,  y_marker, 'r-');
        set( gca, 'xticklabel', ' ' );
        % set the same scale for all plot
        xlim([x_length/13,x_length/1.8]);
        ylim([0,max_count/(timebin/1000)]);
     end 

    xoffset=xoffset+0.6;
%        xoffset=xoffset+0.46;
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return;

