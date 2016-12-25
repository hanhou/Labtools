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
timebin=50;
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
        % get rid off -90 and 90 cases
        if (sum(select) > 0)
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
            % get the average of the total same conditions if repetion is > 1
       %     if (length(act_found) > 1);
            dim=size(temp_count);
            if dim(1) > 1;
               count_y{i,k} = mean(temp_count);
            else
               count_y{i,k}= temp_count;     % for only one repetition cases
            end               
         else
            resp{k}(i) = 0; 
            count_y{i,k}=count_y{1,k};
         end   
         % normalize count_y
         if max(count_y{i,k})~=0;
            count_y_norm{i,k}=count_y{i,k} / max(count_y{i,k});
         else
            count_y_norm{i,k}=0;
         end
    end    
end

% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
% define figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', '3D Direction Tuning');
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
    text(-47,-40, 'Azim: 0       45       90        135        180        225        270       315');  
    text(25,-40, 'Translation');
    axis off;
    hold on;
    for i=1:length(unique_azimuth)                  
        for j=1:length(unique_elevation)
            axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]);                             
            bar( x_time,count_y{i,j,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
            hold on;
            plot( x_start, y_marker, 'r-');
            plot( x_stop,  y_marker, 'r-');
            set( gca, 'xticklabel', ' ' );
            % set the same scale for all plot
            xlim([0,x_length]);
            ylim([0,max_count]);
        end    
    end 

    xoffset=xoffset+0.46;
    
end
% %---------------------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s'];
% for i = 1 : x_length * 3
%      sprint_txt = [sprint_txt, ' %1.2f'];    
% end
% %buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_45(1,:)/8, count_y_45(2,:)/8,count_y_90(1,:)/8, count_y_90(2,:)/8,count_y_135(1,:)/8, count_y_135(2,:)/8, count_y_180(1,:), count_y_180(2,:));  
% buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_max{3} );  
% %buff = sprintf(sprint_txt, FILE, count_trial_beg,count_trial_end ); 
% 
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\DirectionTuning3D_PSTH.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);

return;

