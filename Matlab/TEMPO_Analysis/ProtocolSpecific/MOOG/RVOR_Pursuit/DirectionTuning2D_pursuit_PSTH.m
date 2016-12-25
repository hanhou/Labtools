%----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Yong, 6/27/03  Modified by Katsu, 2/1/06
%-----------------------------------------------------------------------------------------------------------------------

function Rotation_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type

temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG); 

temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);

temp_spike_data = data.spike_data(SpikeChan,:);
temp_spike_rates = data.spike_rates(SpikeChan, :); 


%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_elevation == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_elevation);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 

elevation = temp_elevation(~null_trials & select_trials);
fp_rotate = temp_fp_rotate(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);

% stim_duration = length(temp_spike_data)/length(temp_azimuth);
% spike_data = data.spike_data(1, ((BegTrial-1)*stim_duration+1):EndTrial*stim_duration);
spike_rates= temp_spike_rates(~null_trials & select_trials);
% notice that this bad_trials is the number without spon trials 

% unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_fp_rotate = munique(fp_rotate');
unique_stim_type = munique(stim_type');

h_title{1,1}='Head-fixed';
h_title{2,1}='World-fixed';
h_title{3,1}='Pursuit only';
h_title{1,2}='visual with fix';
h_title{2,2}='Back Move pursuit';
h_title{3,2}='Back stable pursuit';

% add parameters here
% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% find spontaneous trials which azimuth,elevation,fp_rotate=-9999
spon_found = find(null_trials==1);     

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_elevation);
Discard_trials = find(null_trials==1 | temp_spike_rates > 3000 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 99;
end
spike_data = temp_spike_data( temp_spike_data~=99 );

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;

for o=1: length(unique_stim_type)
for k=1: length(unique_fp_rotate)
    for j=1:length(unique_elevation)
%         for i=1: length(unique_azimuth)
            select = logical( (elevation==unique_elevation(j)) & (fp_rotate==unique_fp_rotate(k))& (stim_type==unique_stim_type(o)) );            
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                resp{k,o}(j) = mean(spike_rates(select));
                act_found = find( select==1 );
                % count spikes per timebin on every same condition trials
                for repeat=1:length(act_found) 
                    for n=1:(x_length)
                        temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                        time_step=time_step+timebin;
                    end
                    time_step=1;                    
                end
%  tameshini%               count_y_trial{k,j}(:,:) = temp_count;  % each trial's PSTH 
                % get the average of the total same conditions if repetion is > 1
           %     if (length(act_found) > 1);
                dim=size(temp_count);
                if dim(1) > 1;
                   count_y{j,k,o} = mean(temp_count);
                else
                   count_y{j,k,o}= temp_count;     % for only one repetition cases
                end
               
             else
                resp{k,o}(j) = 0; 
                count_y{j,k,o}=count_y{1,k,o};
             end   
             % normalize count_y
             if max(count_y{j,k,o})~=0;
                count_y_norm{j,k,o}=count_y{j,k,o} / max(count_y{j,k,o});
             else
                count_y_norm{j,k,o}=0;
             end
%         end
    end  
    % now find the peak
    [row_max] = find( resp{k,o}(:)==max(max(resp{k,o}(:))) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k,o}=row_max(1);
    
    if max(count_y{row_max(1), k,o})~=0;
       count_y_max{k,o} = count_y{row_max(1), k,o} / max(count_y{row_max(1), k,o});
    else
       count_y_max{k,o} =0;
    end
    % find the largest y to set scale later
    if max(count_y{row_max(1), k,o}) > max_count
        max_count = max(count_y{row_max(1), k,o});
    end
end
end


% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];

% Aihua's figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', 'RVOR_pursuit');
orient landscape;
% title(FILE);
axis off;

% now plot
% output some text 
axes('position',[0 0 1 1]); 
xlim([-50,50]);
ylim([-50,50]);

for o=1:length(unique_stim_type)
for k=1:length(unique_fp_rotate)   
    text(-49,60-13*k-50*(o-1), h_title{k,o} );
end
end
%text(-47,-35, 'Azim:             270                  225                  180                  135                  90                  45                  0                  315                  270');
text(-47,-35, 'Azim:                 0                    45                   90                  135                   180                   225                  270                 315              0');
text(-10,-40, 'RVOR / pursuit'); 
text(-40,-40, FILE); 
axis off;
hold on;

for o=1:length(unique_stim_type)
for k=1:length(unique_fp_rotate)

      for j=1:length(unique_elevation)+1
       axes('position',[(0.09+0.09*j) (1.0-0.1*k-0.4*(o-1))  0.085 0.085])
                                        % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,  
         if (j <= 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
            bar( x_time,count_y{j,k,o}(1,:) );    % j=9=1 == 0 degrees
        else
            bar( x_time,count_y{1,k,o}(1,:) ); 
        end                                
%             bar( x_time,count_y{1,j,k}(1,:) );    %no additional 9th graph% which is forward motion and the lateral edges correspond to 0 deg which is backward motion
        hold on;
        plot( x_start, y_marker, 'r-');
        plot( x_stop,  y_marker, 'r-');
        set( gca, 'xticklabel', ' ' );        
        % set the same scale for all plot
       xlim([0,x_length]);
       ylim([0,max_count]);
    end    
end 
end
% %--------------------------------------------------------------------------
% %
return;
% 
