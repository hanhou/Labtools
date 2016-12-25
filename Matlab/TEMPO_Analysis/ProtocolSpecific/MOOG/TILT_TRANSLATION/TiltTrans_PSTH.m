%----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Yong, 6/27/03 Modified by Katsu 2/14/06
%-----------------------------------------------------------------------------------------------------------------------

function Rotation_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(TT_MODE,:,MOOG); 
temp_spike_data = data.spike_data(1,:);
temp_spike_rates = data.spike_rates(SpikeChan, :); 


%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_elevation);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);

% stim_duration = length(temp_spike_data)/length(temp_azimuth);
% spike_data = data.spike_data(1, ((BegTrial-1)*stim_duration+1):EndTrial*stim_duration);
spike_rates= temp_spike_rates(~null_trials & select_trials);
% notice that this bad_trials is the number without spon trials 

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');

condition_num = stim_type;
h_title{1}='Tilt+Trans';
h_title{2}='Tilt-Trans';
h_title{3}='Tilt only';
h_title{4}='Trans only ';
% h_title{1}='Head-fixed';
% h_title{2}='World-fixed';
% h_title{3}='Pursuit only';
unique_condition_num = munique(condition_num');

% add parameters here
% timebin for plot PSTH
timebin=250;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% find spontaneous trials which azimuth,elevation,stim_type=-9999
spon_found = find(null_trials==1);     

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_elevation);
Discard_trials = find(null_trials==1 | temp_spike_rates > 3000 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
for k=1: length(unique_condition_num)
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );            
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                resp{k}(j,i) = mean(spike_rates(select));
                act_found = find( select==1 );
                % count spikes per timebin on every same condition trials
                for repeat=1:length(act_found) 
                    for n=1:(x_length)
                        temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                        time_step=time_step+timebin;
                    end
                    time_step=1;                    
                end
                count_y_trial{k,i,j}(:,:) = temp_count;  % each trial's PSTH 
                % get the average of the total same conditions if repetion is > 1
           %     if (length(act_found) > 1);
                dim=size(temp_count);
                if dim(1) > 1;
                   count_y{i,j,k} = mean(temp_count);
                else
                   count_y{i,j,k}= temp_count;     % for only one repetition cases
                end
               
             else
                resp{k}(j,i) = 0; 
                count_y{i,j,k}=count_y{1,j,k};
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
% %--------------------------------------------------------------------------
% % compare whether there is delay over recording sesseion, cut into two parts
% repetition = floor( length(spike_rates)/78); % take minimum repetition
% if repetition >=5 % only include data more than 5 repetitions
% 	max_trial = max(max( mean(count_y_trial{2,col_m{2},row_m{2}}(1:3,:)),mean(count_y_trial{2,col_m{2},row_m{2}}(4:end,:)) ));
% 	count_trial_beg = mean(count_y_trial{2,col_m{2},row_m{2}}(1:3,:)) / max_trial;
% 	count_trial_end = mean(count_y_trial{2,col_m{2},row_m{2}}(4:end,:)) / max_trial;
%     repetition
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculate PSTH for each direction according to the eccentricity from the
% % maximum direction, 45,90,135 and 180, all 5 variables
% % techniquelly, given any maximum angle, go through all the directions to
% % pick up any angles with expected difference. us norm, dot to calculate
% % the angle between two vectors in 3D. 
% count_y_45(3,100)=0;
% count_y_90(3,100)=0;
% count_y_135(3,100)=0;
% for k=1: length(unique_condition_num)
%     for j=1:length(unique_elevation)
%         for i=1: length(unique_azimuth)
%             select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );            
%             if (sum(select) > 0)
%                 [x,y,z]=sph2cart( unique_azimuth(i)*3.14159/180,unique_elevation(j)*3.14159/180,1 );
%                 direction_c=[x,y,z];
%                 [xm,ym,zm]=sph2cart( unique_azimuth(col_m{k})*3.14159/180,unique_elevation(row_m{k})*3.14159/180,1 );
%                 direction_m=[xm,ym,zm];
%                 diff_angle = 180*acos( dot(direction_c,direction_m) )/3.14159;
%                 if diff_angle > -1 & diff_angle < 1
%                     count_y_0(k,:) = count_y_norm{i,j,k};   % actually this is the same to count_y_max
%                 elseif diff_angle > 44 & diff_angle < 46
%                     count_y_45(k,:) = count_y_45(k,:) + count_y_norm{i,j,k};
%                 elseif diff_angle > 89 & diff_angle < 91
%                     count_y_90(k,:) = count_y_90(k,:) + count_y_norm{i,j,k};
%                 elseif diff_angle > 134 & diff_angle < 136
%                     count_y_135(k,:) = count_y_135(k,:) + count_y_norm{i,j,k};
%                 elseif diff_angle > 179
%                     count_y_180(k,:) = count_y_norm{i,j,k};
%                 end                
%             end
%         end
%     end
% end

% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
%---------------------------------------------------------------------------------------------------------------------------------------
% Rotation Azimuth Scale figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', 'Tilt_Translation (Rotation scale)');
orient landscape;
title(FILE);
axis off;

% output some text 
axes('position',[0 0 1 0.9]); 
xlim([-50,50]);
ylim([-50,50]);
for i=1:length(unique_condition_num)   
    text(-48,65-22.5*i, h_title{i} );
end
text(-50,-35, 'Rotation Azimuth:        270                  225                  180                  135                   90                     45                  0                  315                270');
% text(-50,-35, 'Rotation Azimuth:        0                  45                  90                 135                   180                   225                  270                 315');
text(-50,-37, 'Tilt/Trans Direction:       0                  315                    270                   225                 180                     135                 90                  45                0');
% text(-50,-37, 'Tilt/Trans Direction:   270                315                   0                   45                    90                   135                  180                 225');
text(-10,-40, 'Tilt / Translation (Rotation Azimuth Scale)');
axis off;
hold on;

% for i=1:length(unique_azimuth)     
%     for j=1:length(unique_condition_num)
%        axes('position',[0.01+0.09*i  (1.0-0.2*j)  0.085 0.085])
%                                         % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
%             bar( x_time,count_y{i,1,j}(1,:) );    % which is forward motion and the lateral edges correspond to 0 deg which is backward motion
%         hold on;
%         plot( x_start, y_marker, 'r-');
%         plot( x_stop,  y_marker, 'r-');
%         set( gca, 'xticklabel', ' ' );        
%         % set the same scale for all plot
%        xlim([0,x_length]);
%        ylim([0,max_count]);
%     end    
% end 

% for i=1:length(unique_azimuth)+1      % aizmuth 270 are plotted two times in order to make circular data
%     for j=1:length(unique_condition_num)
%        axes('position',[0.01+0.09*i  (1.0-0.2*j)  0.085 0.085])
%         if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
%             bar( x_time,count_y{8-i,j}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
%         elseif(i==8)
%             bar( x_time,count_y{i,j}(1,:) ); 
%         else
%             bar( x_time,count_y{7,j}(1,:) ); 
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
for k=1:length(unique_condition_num)      % aizmuth 270 are plotted two times in order to make circular data
    for i=1:length(unique_azimuth)+1
       axes('position',[0.01+0.09*i  (1.0-0.2*k)  0.085 0.085])
        if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
            bar( x_time,count_y{8-i,1,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
        elseif(i==8)
            bar( x_time,count_y{i,1,k}(1,:) ); 
        else
            bar( x_time,count_y{7,1,k}(1,:) ); 
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
%--------------------------------------------------------------------------------------------------------------------------------------
% Translation Scale figure
figure(3);
set(3,'Position', [5,5 1000,700], 'Name', 'Tilt_Translation (Translation Scale)');
orient landscape;
title(FILE);
axis off;

% output some text 
axes('position',[0 0 1 0.9]); 
xlim([-50,50]);
ylim([-50,50]);
for i=1:length(unique_condition_num)   
    text(-48,65-22.5*i, h_title{i} );
end
text(-50,-37, 'Rotation Azimuth:       180                  135                    90                    45                     0                     315                 270                225               180');
% text(-50,-35, 'Rotation Azimuth:        0                  45                  90                 135                   180                   225                  270                 315');
text(-50,-35, 'Tilt/Trans Direction:     270                  225                   180                   135                   90                     45                   0                 315              270');
% text(-50,-37, 'Tilt/Trans Direction:   270                315                   0                   45                    90                   135                  180                 225');
text(-10,-40, 'Tilt / Translation Direction Scale');
axis off;
hold on;

% for i=1:length(unique_azimuth)     
%     for j=1:length(unique_condition_num)
%        axes('position',[0.01+0.09*i  (1.0-0.2*j)  0.085 0.085])
%                                         % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
%             bar( x_time,count_y{i,1,j}(1,:) );    % which is forward motion and the lateral edges correspond to 0 deg which is backward motion
%         hold on;
%         plot( x_start, y_marker, 'r-');
%         plot( x_stop,  y_marker, 'r-');
%         set( gca, 'xticklabel', ' ' );        
%         % set the same scale for all plot
%        xlim([0,x_length]);
%        ylim([0,max_count]);
%     end    
% end 

for k=1:length(unique_condition_num)      % aizmuth 270 are plotted two times in order to make circular data
    for i=1:length(unique_azimuth)+1
       axes('position',[0.01+0.09*i  (1.0-0.2*k)  0.085 0.085])
        if (i < 6)                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
            bar( x_time,count_y{6-i,1,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
%         elseif(i>=2)
%             bar( x_time,count_y{10-i,j}(1,:) ); 
        else
            bar( x_time,count_y{14-i,1,k}(1,:) ); 
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
% Yong's figure


% define figure
% figure(2);
% set(2,'Position', [5,5 1000,700], 'Name', '3D Rotation Tuning');
% orient landscape;
% title(FILE);
% axis off;
% 
% xoffset=0;
% yoffset=0;
% 
% % now plot
% for k=1: length(unique_condition_num) 
%     
%     if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
%         yoffset = yoffset-0.42;
%         xoffset = 0;
%     end
%     % output some text 
%     axes('position',[0 0 1 0.9]); 
%     xlim([-50,50]);
%     ylim([-50,50]);
%     text(-30+xoffset*100,52+yoffset*110, h_title{k} );
%     text(-47,-40, 'Azim: 270        225       180       135         90         45         0         315        270');
%     text(25,-40, 'Rotation');
%     axis off;
%     hold on;
% %     for i=1:length(unique_azimuth)+1                    % aizmuth 270 are plotted two times in order to make circular data
%         for j=1:length(unique_elevation)
%             axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
% %             if (i <  )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
% %                 bar( x_time,count_y{8-i,j,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
% %             elseif(i==8)
%                 %bar( x_time,count_y{i,j,k}(1,:) ); 
%                 bar( x_time,count_y{1,j,k}(1,:) );
% %             else
% %                 bar( x_time,count_y{7,j,k}(1,:) ); 
% %             end
%             hold on;
%             plot( x_start, y_marker, 'r-');
%             plot( x_stop,  y_marker, 'r-');
%             set( gca, 'xticklabel', ' ' );
%             % set the same scale for all plot
%             xlim([0,x_length]);
%             ylim([0,max_count]);
%         end    
% %     end 
% 
%     xoffset=xoffset+0.46;
%     
% end
% %---------------------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file
% %sprint_txt = ['%s'];
% %for i = 1 : x_length * 3
% %     sprint_txt = [sprint_txt, ' %1.2f'];    
% %end
% %buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_45(1,:)/8, count_y_45(2,:)/8,count_y_90(1,:)/8, count_y_90(2,:)/8,count_y_135(1,:)/8, count_y_135(2,:)/8, count_y_180(1,:), count_y_180(2,:));  
% %buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_max{3} );  
% %buff = sprintf(sprint_txt, FILE, count_trial_beg,count_trial_end ); 
% 
% %outfile = [BASE_PATH 'ProtocolSpecific\MOOG\DirectionTuning3D_PSTH.dat'];
% %printflag = 0;
% %if (exist(outfile, 'file') == 0)    %file does not yet exist
% %    printflag = 1;
% %end
% %fid = fopen(outfile, 'a');
% %if (printflag)
% %    fprintf(fid, 'FILE\t');
% %    fprintf(fid, '\r\n');
% %end
% %fprintf(fid, '%s', buff);
% %fprintf(fid, '\r\n');
% %fclose(fid);
% 
% return;





% function  DirectionTuning2D_pursuit_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %KT 02-02-2006
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% TEMPO_Defs;
% Path_Defs;
% ProtocolDefs; 
% 
% %get the column of values for elevation, fp_rotate
% temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
% temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG);
% temp_spike_data = data.spike_data(SpikeChan,:);
% temp_spike_rates = data.spike_rates(SpikeChan,:);   %Timebin: The duration of one trial 
% 
% sti_duration=length(temp_spike_data)/length(temp_spike_rates);
% for i=1:length(temp_spike_rates)   
%     spike_data(i,:)=temp_spike_data((i-1)*sti_duration+1:i*sti_duration);
% end
% 
% real_trials=find(temp_elevation~=data.one_time_params(NULL_VALUE));
% elevation=temp_elevation(real_trials)
% %delete;
% unique_elevation = munique(elevation');
% fp_rotate = temp_fp_rotate(real_trials);
% unique_fp_rotate = munique(fp_rotate');
% 
% h_title{1}='Head-fixed';
% h_title{2}='World-fixed';
% h_title{3}='Pursuit only';
% 
% % add parameters here
% timebin=50;% timebin for plot PSTH
% x_time=1:sti_duration/timebin;%x-axis for plot PSTH
% 
% %count spike data from raster data (spike _data)
% for j=1:length(unique_fp_rotate)
%     for i=1:length(unique_elevation)
%         select=find(elevation==unique_elevation(i) & (fp_rotate==unique_fp_rotate(j)));
%         trial_mean=mean(spike_data(select,:));   
%         for n=1:length(x_time)
%             temp_count(1,n)=sum(trial_mean(1,(n-1)*timebin+1: n*timebin));            
%         end
%         count_y{i,j}=temp_count; 
%        count_y_max(i,j)=max(temp_count);
%         clear temp_count
%         if max(count_y{i,j})~=0;
%             count_y_norm{i,j}=count_y{i,j} / max(count_y{i,j});
%         else
%             count_y_norm{i,j}=0;
%         end
%     end
% end
% max_count=max(max(count_y_max));
% 
% 
% % plot PSTH now
% % get the largest count_y so that make the scale in each figures equal    
% % plot two lines as stimulus start and stop marker
% x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
% x_stop = [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
% y_marker=[0,max_count];
% % define figure
% figure(2);
% set(2,'Position', [5,5 1000,700], 'Name', 'RVOR_pursuit');
% orient landscape;
% title(FILE);
% axis off;
% 
% % now plot
% % output some text 
% axes('position',[0 0 1 0.9]); 
% xlim([-50,50]);
% ylim([-50,50]);
% for i=1:length(unique_fp_rotate)   
%     text(-48,65-22.5*i, h_title{i} );
% end
% %text(-47,-35, 'Azim:             270                  225                  180                  135                  90                  45                  0                  315                  270');
% text(-47,-35, 'Azim:                0                    45                   90                  135                    180                    225                  270                  315                  0');
% text(-10,-40, 'RVOR / pursuit'); 
% axis off;
% hold on;
% 
% for i=1:length(unique_elevation)+1      % aizmuth 0 are plotted two times in order to make circular data
%     for j=1:length(unique_fp_rotate)
%        axes('position',[0.01+0.09*i  (1.0-0.2*j)  0.085 0.085])
%         if (i < 9 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
%             bar( x_time,count_y{i,j}(1,:) );    % which is forward motion and the lateral edges correspond to 0 deg which is backward motion
%         else
%             bar( x_time,count_y{1,j}(1,:) ); 
%         end
%         hold on;
%         plot( x_start, y_marker, 'r-');
%         plot( x_stop,  y_marker, 'r-');
%         set( gca, 'xticklabel', ' ' );        
%         % set the same scale for all plot
%        xlim([0,sti_duration/timebin]);
%        ylim([0,max_count]);
%     end    
% end 
% 
% % for i=1:length(unique_elevation)+1      % aizmuth 270 are plotted two times in order to make circular data
% %     for j=1:length(unique_fp_rotate)
% %        axes('position',[0.01+0.09*i  (1.0-0.2*j)  0.085 0.085])
% %         if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
% %             bar( x_time,count_y{8-i,j}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
% %         elseif(i==8)
% %             bar( x_time,count_y{i,j}(1,:) ); 
% %         else
% %             bar( x_time,count_y{7,j}(1,:) ); 
% %         end
% %         hold on;
% %         plot( x_start, y_marker, 'r-');
% %         plot( x_stop,  y_marker, 'r-');
% %         set( gca, 'xticklabel', ' ' );        
% %         % set the same scale for all plot
% %        xlim([0,sti_duration/timebin]);
% %        ylim([0,max_count]);
% %     end    
% % end 
% 
% 
% % for i=1:length(unique_azimuth)+1  
% %     for j=1:length(unique_TT_MODE)
% %         if (i < 8 )              
% %             figure(2);subplot(4,9, (i-1)*4+j), bar( x_time,count_y{8-i,j}(1,:) );  
% %         elseif(i==8)
% %             figure(2);subplot(4,9, (i-1)*4+j), bar( x_time,count_y{i,j}(1,:) ); 
% %         else
% %             figure(2);subplot(4,9, (i-1)*4+j), bar( x_time,count_y{7,j}(1,:) ); 
% %         end
% %         hold on;
% %         plot( x_start, y_marker, 'r-');
% %         plot( x_stop,  y_marker, 'r-');
% %         set( gca, 'xticklabel', ' ' );       
% %         % set the same scale for all plot
% %        xlim([0,x_length]);
% %        ylim([0,max_count]);
% %     end    
% % end 
% 
return;