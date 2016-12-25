% Modify %%%%%%%%DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	Katsu
%-----------------------------------------------------------------------------------------------------------------------
function Horizontal_only(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
%
temp_spike_data = data.spike_data(SpikeChan,:);
%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_tri ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
MaxAmp=max(amplitude)
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
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



% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found))
% added by Katsu 111606
spon_std = std(temp_spike_rates(spon_found))

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(null_trials==1 | temp_spike_rates > 3000 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );

% count spikes from raster data (spike_data)

zero_elevation=[0];

max_count = 1;
time_step=1;
for k=1: length(unique_condition_num)
%     for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==zero_elevation) & (condition_num==unique_condition_num(k)) );            
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
%         end
    end  
    % now find the peak
    [row_max, col_max] = find( resp{k}(:)==max(resp{k}(:)) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k}=row_max(1);
    col_m{k}=col_max(1);
    if max(count_y{col_max(1), k})~=0;
       count_y_max{k} = count_y{col_max(1), k} / max(count_y{col_max(1), k});
    else
       count_y_max{k} =0;
    end
    % find the largest y to set scale later
    if max(count_y{col_max(1), k}) > max_count
        max_count = max(count_y{col_max(1), k});
    end
end





%-------------------------------------------------------------------------
%ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-14) * length(unique_condition_num) + 1;
repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);

% first parse raw data into repetitions, including null trials
for q = 1:repetitions
   azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   condition_num_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
   spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
end

% zero_elevation=[0]


%% ADD CODE HERE FOR PLOTTING
resp_mat = [];
for i=1:length(unique_azimuth)
%     for j=1:length(unique_elevation)
        for k=1: length(unique_condition_num)
%                                                                              unique_elevation(3)=0 deg
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==zero_elevation ) & (condition_num==unique_condition_num(k)) );
            if (sum(select) > 0)                
                resp_mat(k, i) = mean(spike_rates(select));
%                 
                for q = 1 : repetitions               
                    spike_temp = spike_rates(select);                 %For one-way-ANOVA
                    resp_mat_trial{k}(q, i) = spike_temp( q );     % 
                end
                
                resp_mat_std(k, i) = std(spike_rates(select));     % calculate std between trials for later DSI usage
                resp_mat_ste(k, i) = resp_mat_std(k, i)/ sqrt(length(find( (azimuth==unique_azimuth(i)) & (elevation==zero_elevation)&(condition_num==unique_condition_num(k)) )) );
            else
%                resp_mat_trial{k}(t, j, i) = 0;
                resp_mat(k, i) = resp_mat(k,1);
%                 resp_mat_vector(k,1) =0; % for vector sum only
                resp_mat_std(k, i) = 0;
                resp_mat_ste(k, i) = 0;
            end
        end        
%     end
end


% vectorsum and calculate preferred direction
unique_elevation_s(1:length(unique_azimuth)) = 0;
for k = 1: length(unique_condition_num)
    [az(k), el(k), amp(k)] = vectorsumAngle( resp_mat(k, :), unique_azimuth, unique_elevation_s );
end

%% creat a real 3-D based plot where the center correspond to forward and both lateral edges correspond to backward
%%%% Usually, axis azimuth from left is 270-225-180-135-90--0---90 %%%% 
resp_mat_tran(:,1) = resp_mat(:,7);
resp_mat_tran(:,2) = resp_mat(:,6);
resp_mat_tran(:,3) = resp_mat(:,5);
resp_mat_tran(:,4) = resp_mat(:,4);
resp_mat_tran(:,5) = resp_mat(:,3);
resp_mat_tran(:,6) = resp_mat(:,2);
resp_mat_tran(:,7) = resp_mat(:,1);
resp_mat_tran(:,8) = resp_mat(:,8);
resp_mat_tran(:,9) = resp_mat_tran(:,1);

resp_mat_ste_tran(:,1) = resp_mat_ste(:,7);
resp_mat_ste_tran(:,2) = resp_mat_ste(:,6);
resp_mat_ste_tran(:,3) = resp_mat_ste(:,5);
resp_mat_ste_tran(:,4) = resp_mat_ste(:,4);
resp_mat_ste_tran(:,5) = resp_mat_ste(:,3);
resp_mat_ste_tran(:,6) = resp_mat_ste(:,2);
resp_mat_ste_tran(:,7) = resp_mat_ste(:,1);
resp_mat_ste_tran(:,8) = resp_mat_ste(:,8);
resp_mat_ste_tran(:,9) = resp_mat_ste_tran(:,1);
%
% %%%%%%%

%Translation direction, so Rotation Azimuth is 180 135 90 45 0 315 270 225 180
xoffset=0;
yoffset=0;
figure(4);
set(4,'Position', [5,15 980,650], 'Name', 'Horizontal only');
orient landscape;
%set(0, 'DefaultAxesXTickMode', 'manual', 'DefaultAxesYTickMode', 'manual', 'DefaultAxesZTickMode', 'manual');

% unique_elevation_270=[270;225;180;15;90;45;0;315;270];%by KT

% spon_elevation = min(unique_elevation) : 1 : max(unique_elevation);never plot by dots I want to change to line

%re-arrangement of x-axis, start at 270;225;180;15;90;45;0;315;270
%As for rotation Azi. 180;135;90;45;0;315;270;225;180


%re-a%re-arrangement of az (preferreed direction)
for k=1:length(unique_condition_num)
    if az(k)<270
        az_270(k)=abs(az(k)-270);
        az_270(k)=((az_270(k).*8)./360)+1;% x-axis convert to 1----9
    else
        az_270(k)=360+270-az(k);
        az_270(k)=((az_270(k).*8)./360)+1;% x-axis convert to 1----9
    end
end
    
%NOW! plot errorbar and plot
for k=1: length(unique_condition_num)     
    if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
        yoffset = yoffset-0.4;
        xoffset = 0;
    end
    axes('position',[0.11+xoffset 0.44+yoffset 0.32 0.24]);
    x=[1,2,3,4,5,6,7,8,9];%temporaly x axis by number, later, change to 270-----0---270
    errorbar(x, resp_mat_tran(k,:), resp_mat_ste_tran(k,:), 'bo-' );
    
    hold on;
    plot([az_270(k),az_270(k)],[min(resp_mat(k,:)),max(resp_mat(k,:))],'r-');
    hold on;
    %plot(spon_elevation, spon_resp, 'r-');
    plot([1,9],[spon_resp,spon_resp],'r-');
    hold on;
    ylabel('spikes/s');
    xlabel('Translation Azimuth');
    xlim( [1, length(unique_azimuth)+1] ); 
    set(gca,'XTickMode','manual');
    set(gca,'xtick',[1,2,3,4,5,6,7,8,9]);
    set(gca,'xticklabel','270|225|180|135|90|45|0|315|270');
    title(h_title{k});
    %set(gca, 'xtick',[unique_elevation270]);
    hold on;
    xoffset=xoffset+0.48;    
end

    
    % check significance by one-way-anova
    for k=1:length(unique_condition_num)
    
    p_1D(k) = anova1(resp_mat_trial{k},'','off');
    
    end
    
% show file name and some values in text
axes('position',[0.05,0.75, 0.9,0.1] );
xlim( [0,100] );
ylim( [0,length(unique_condition_num)] );
text(0, length(unique_condition_num), FILE);
text(13,length(unique_condition_num),'spont                  min             max                        prefer Azimuth                  p-ANOVA');

for k=1:length(unique_condition_num) 
% Now write
    text(0,length(unique_condition_num)-k,h_title{k});
    text(15,length(unique_condition_num)-k, num2str( spon_resp));
    text(27,length(unique_condition_num)-k, num2str( min(resp_mat(k,:)) ) );
    text(37,length(unique_condition_num)-k, num2str( max(resp_mat(k,:)) ) );
    text(52,length(unique_condition_num)-k, num2str(az(k)) );
    text(65,length(unique_condition_num)-k, num2str(p_1D(k)) );
end
axis off;
% 
% % %%%%%%%%%%%%%%%%%%%%%
% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
% define figure
% figure(2);
% set(2,'Position', [5,5 1000,700], 'Name', '3D Direction Tuning');
% orient landscape;
% title(FILE);
% axis off;

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
%     text(-30+xoffset*100,52+yoffset*110, h_title{k} );
    text(-47,-40, 'Azim: 270       225       180        135        90        45        0        315        270');  
    text(25,-40, 'Translation');
    axis off;
    hold on;
    for i=1:length(unique_azimuth)+1                    % aizmuth 270 are plotted two times in order to make circular data
%         for j=1:length(unique_elevation)
            axes('position',[0.05*i+0.01+xoffset (0.32-0.07)+yoffset 0.045 0.045]); 
            if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
                bar( x_time,count_y{8-i,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
            elseif(i==8)
                bar( x_time,count_y{i,k}(1,:) ); 
            else
                bar( x_time,count_y{7,k}(1,:) ); 
            end
            hold on;
            plot( x_start, y_marker, 'r-');
            plot( x_stop,  y_marker, 'r-');
            set( gca, 'xticklabel', ' ' );
            % set the same scale for all plot
            xlim([0,x_length]);
            ylim([0,max_count]);
%         end    
    end 

    xoffset=xoffset+0.46;
    
end

%---------------------------------------------------------------------------------------

return;