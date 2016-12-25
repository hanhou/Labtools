function  TiltTrans_PSTH_Old(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TiltTrans_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
%-ACH 09-29-2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TEMPO_Defs;
Path_Defs;
ProtocolDefs; 

%get the column of values for azimuth, TT_MODE
temp_azimuth=data.moog_params(ROT_AZIMUTH,:,MOOG);
temp_TT_MODE=data.moog_params(TT_MODE,:,MOOG);
temp_spike_data = data.spike_data(SpikeChan,:);
temp_spike_rates = data.spike_rates(SpikeChan,:);   %Timebin: The duration of one trial 

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)));

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
azimuth = temp_azimuth(~null_trials & select_trials);
TT_MODE=temp_TT_MODE(~null_trials & select_trials);
spike_rates= temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_TT_MODE=munique(TT_MODE');
h_title{1}='Tilt+Trans';
h_title{2}='Tilt-Trans';
h_title{3}='Tilt only';
h_title{4}='Trans only ';

% add parameters here
timebin=50;% timebin for plot PSTH
frequency=length(temp_spike_data)/length(select_trials);  % sample frequency depends on test duration
x_length = frequency/timebin;% length of x-axis
x_time=1:(frequency/timebin);% x-axis for plot PSTH

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

for j=1:length(unique_TT_MODE)
    for i=1:length(unique_azimuth)
        select=logical((azimuth==unique_azimuth(i))&(TT_MODE==unique_TT_MODE(j)));
        if(sum(select)>0)
            resp(j,i)=mean(spike_rates(select));
            act_found=find(select==1);            
            for repeat=1:length(act_found)
                for n=1:(x_length)
                    temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                    time_step=time_step+timebin;
                end
                time_step=1;   
            end
            count_y_trial{i,j}(:,:) = temp_count;  % each trial's PSTH 
            dim=size(temp_count);
            if dim(1) > 1;
                count_y{i,j} = mean(temp_count);
            else
                count_y{i,j}= temp_count;     % for only one repetition cases
            end 
        else
            resp(j,i) = 0; 
            count_y{i,j}=count_y{1,j};
        end 
        if max(count_y{i,j})~=0;
            count_y_norm{i,j}=count_y{i,j} / max(count_y{i,j});
        else
            count_y_norm{i,j}=0;
        end
    end
end

% now find the peak
[row_max, col_max] = find( resp(:,:)==max(max(resp(:,:))) );
% it is likely there are two peaks with same magnitude, choose the first one arbitraly
row_m=row_max(1);
col_m=col_max(1);
if max(count_y{col_max(1), row_max(1)})~=0;
    count_y_max = count_y{col_max(1), row_max(1)} / max(count_y{col_max(1), row_max(1)});
else
    count_y_max =0;
end
% find the largest y to set scale later
if max(count_y{col_max(1), row_max(1)}) > max_count
    max_count = max(count_y{col_max(1), row_max(1)});
end


% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
% define figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', 'Tilt_Translation');
orient landscape;
title(FILE);
axis off;

% now plot
% output some text 
axes('position',[0 0 1 0.9]); 
xlim([-50,50]);
ylim([-50,50]);
for i=1:length(unique_TT_MODE)   
    text(-48,65-22.5*i, h_title{i} );
end
text(-47,-35, 'Azim:             270                  225                  180                  135                  90                  45                  0                  315                  270');
text(-10,-40, 'Tilt / Translation'); 
axis off;
hold on;

for i=1:length(unique_azimuth)+1      % aizmuth 270 are plotted two times in order to make circular data
    for j=1:length(unique_TT_MODE)
       axes('position',[0.01+0.09*i  (1.0-0.2*j)  0.085 0.085])
        if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
            bar( x_time,count_y{8-i,j}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
        elseif(i==8)
            bar( x_time,count_y{i,j}(1,:) ); 
        else
            bar( x_time,count_y{7,j}(1,:) ); 
        end
        hold on;
        plot( x_start, y_marker, 'r-');
        plot( x_stop,  y_marker, 'r-');
        set( gca, 'xticklabel', ' ' );        
        % set the same scale for all plot
       xlim([0,x_length]);
       ylim([0,max_count]);
    end    
end 


% for i=1:length(unique_azimuth)+1  
%     for j=1:length(unique_TT_MODE)
%         if (i < 8 )              
%             figure(2);subplot(4,9, (i-1)*4+j), bar( x_time,count_y{8-i,j}(1,:) );  
%         elseif(i==8)
%             figure(2);subplot(4,9, (i-1)*4+j), bar( x_time,count_y{i,j}(1,:) ); 
%         else
%             figure(2);subplot(4,9, (i-1)*4+j), bar( x_time,count_y{7,j}(1,:) ); 
%         end
%         hold on;
%         plot( x_start, y_marker, 'r-');
%         plot( x_stop,  y_marker, 'r-');
%         set( gca, 'xticklabel', ' ' );       
%         % set the same scale for all plot
%        xlim([0,x_length]);
%        ylim([0,max_count]);
%     end    
% end 

return;

