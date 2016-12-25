%----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Yong, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------

function MOOG_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

disp('hi there');
figure(10);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
TEMPO_Defs;

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_spike_data = squeeze(data.spike_data(SpikeChan,:,:));
temp_event_data = squeeze(data.event_data);
temp_spike_rates = data.spike_rates(SpikeChan, :);    
% 2008-03-23 MLM
temp_azimuth_moog = data.moog_params(HEADING,:,MOOG);
temp_azimuth_cam = data.moog_params(HEADING,:,CAMERAS);


%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth_moog == data.one_time_params(NULL_VALUE)) & (temp_azimuth_cam == data.one_time_params(NULL_VALUE)));
% For now, undo.
null_trials=logical(zeros(size(null_trials)));

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_trials ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
% stim_duration = length(temp_spike_data)/length(temp_azimuth);
% spike_data = data.spike_data(1, ((BegTrial-1)*stim_duration+1):EndTrial*stim_duration);
spike_rates= temp_spike_rates(~null_trials & select_trials);
% notice that this bad_trials is the number without spon trials 

% 2008-03-23 MLM
azimuth_moog = temp_azimuth_moog(~null_trials & select_trials);
azimuth_cam = temp_azimuth_cam(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);
unique_azimuth_moog = munique(azimuth_moog');
unique_azimuth_cam = munique(azimuth_cam');

unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Bimodal';
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
% stim_duration = length(temp_spike_data)/length(temp_azimuth);
% Discard_trials = find(null_trials==1 | trials < BegTrial | trials >EndTrial);
% for i = 1 : length(Discard_trials)
%     temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
% end
% spike_data = temp_spike_data( temp_spike_data~=9999 );
% spike_data( find(spike_data>100) ) = 1; % something is absolutely wrong 
spike_data=temp_spike_data;

% For each trial, compute a PSTH.
% Each trial is 2 seconds. Each bin is 50 ms.
nbins=round(5/(timebin/1000));
count_tmp=zeros(nbins,size(spike_data,2));
for i=1:size(spike_data,2)
    BegInd=min(find(temp_event_data(:,i) == VSTIM_ON_CD));
    EndInd=min(find(temp_event_data(:,i) == VSTIM_OFF_CD));
    if (EndInd - BegInd) > 1999 % If more than 2000 ms, take the middle 2000 ms.
        BegInd = BegInd + floor(((EndInd - BegInd) - 1999)/2);
        EndInd = BegInd + 1999;
    end
    BegInd=1;EndInd=5000;
    temp_count=spike_data( BegInd:EndInd , i );
    count_tmp(:,i)=sum(reshape( temp_count, [ nbins length(temp_count)/nbins ]),2);
end

for j=1:length(unique_azimuth_moog)
    for i=1:length(unique_azimuth_cam)

        select = logical( (azimuth_moog==unique_azimuth_moog(j)) & (azimuth_cam==unique_azimuth_cam(i)) );

        count_y{i,j} = mean(count_tmp(:, select),2);
        
    end
end

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
for j=1:length(unique_azimuth_moog)
    for i=1:length(unique_azimuth_cam)

        select = logical( (azimuth_moog==unique_azimuth_moog(j)) & (azimuth_cam==unique_azimuth_cam(i)) );

        if (sum(select) > 0)
            resp(j,i) = mean(spike_rates(select));
            act_found = find( select==1 );
            % count spikes per timebin on every same condition trials
            for repeat=1:length(act_found)
                for n=1:(x_length)
                    temp_count(repeat,n)=sum(spike_data((time_step):(n*timebin),act_found(repeat)));
                    time_step=time_step+timebin;
                end
                time_step=1;
            end
            count_y_trial{i,j}(:,:) = temp_count;  % each trial's PSTH
            % get the average of the total same conditions if repetion is > 1
            %     if (length(act_found) > 1);
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
        % normalize count_y
        if max(count_y{i,j})~=0;
            count_y_norm{i,j}=count_y{i,j} / max(count_y{i,j});
        else
            count_y_norm{i,j}=0;
        end



        % now find the peak
        [row_max, col_max] = find( resp(:,:)==max(max(resp(:,:))) )
        % it is likely there are two peaks with same magnitude, choose the first one arbitraly
        row_m=row_max(1);
        col_m=col_max(1);
        if max(count_y{col_max(1), row_max(1)})~=0;
            %  count_y_max{k} = count_y{col_max(1), row_max(1), k} / max(count_y{col_max(1), row_max(1), k}); % normalized
            count_y_max = count_y{col_max(1), row_max(1)};
            %count_y_max{k} = count_y{4, 3, k};
        else
            count_y_max =0;
        end
        % find the largest y to set scale later
        if max(count_y{col_max(1), row_max(1)}) > max_count
            max_count = max(count_y{col_max(1), row_max(1)});
        end
        
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
% 
% count_y_45(3,frequency/timebin)=0;
% count_y_90(3,frequency/timebin)=0;
% count_y_135(3,frequency/timebin)=0;
% % count_y_45(3,100)=0;
% % count_y_90(3,100)=0;
% % count_y_135(3,100)=0;
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
    text(-30+xoffset*100,52+yoffset*110, h_title{k} );
    %text(-47,-40, 'Azim: 270       225       180        135        90        45        0        315        270');  
    temp=fliplr(unique_azimuth');unique_azimuth_plot=[temp(2:end) temp(1:2)];clear temp
    text(-47,-40, ['Azim:' num2str(unique_azimuth_plot)]);  
    text(25,-40, 'Translation');
    axis off;
    hold on;
    for i=1:length(unique_azimuth_moog)
        for j=1:length(unique_azimuth_cam)
            axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
            bar( x_time,count_y{i,j}(1,:) );
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Special PSTH with velocity and acceleration information-- AC 02/17/08
figure(3)
t = 0:.05:1.85;
ampl = 0.13;
%sigma approximately equal to .17-.18
num_sigs = 6;
pos = ampl*0.5*(erf(2*num_sigs/3*(t-1)) + 1);
veloc = diff(pos)/0.05;
norm_veloc= veloc./max(veloc);
accel = diff(veloc)/0.05;
norm_accel= accel./max(accel); %normal acceleration
% norm_accel= -accel./max(accel); %flipped acceleration
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
xoffset=0;yoffset=0;

for k=1: length(unique_condition_num) % K = condition 1:vestibular, 2:visual
    for i=1:8+1
        for j=1:5
            figure(k+2);axes('position',[0.1*(i-1)+0.05+xoffset (0.92-0.1*j)+yoffset 0.09 0.09]);
            if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
                %             bar( x_time,count_y{8-i,j,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
                bar( x_time(round(x_start(1,1)):round(x_stop(1,1))),count_y{8-i,j,k}(1,round(x_start(1,1)):round(x_stop(1,1))) );
            elseif(i==8)
                bar( x_time(round(x_start(1,1)):round(x_stop(1,1))),count_y{i,j,k}(1,round(x_start(1,1)):round(x_stop(1,1))) ); 
            else
                bar( x_time(round(x_start(1,1)):round(x_stop(1,1))),count_y{7,j,k}(1,round(x_start(1,1)):round(x_stop(1,1))) );
            end
            %         plot( x_start, y_marker, 'r-','LineWidth',2.0);
            %         plot( x_stop,  y_marker, 'r-','LineWidth',2.0);
            hold on;
%             max_count=3;%Syed
            if (i==5 & j==5)
                plot([ x_start(1,1):(x_stop(1,1)-x_start(1,1))/(length(norm_veloc)-1):x_stop(1,1)],0.5*max_count*(norm_veloc),'r.','LineWidth',2.0);           
                %             plot([ x_start(1,1):(x_stop(1,1)-x_start(1,1))/(length(norm_accel)-1):x_stop(1,1)],0.5*max_count*(norm_accel),'g.','LineWidth',2.0);
                text(0,-1, h_title{k});
            else
            end
            set( gca, 'xticklabel', ' ');
            if (i>1)
                set(gca,'yticklabel',' ');
            end
            % set the same scale for all plot
            %         xlim([0,x_length]);
            xlim([round(x_start(1,1)),round(x_stop(1,1))]);            
            ylim([0,max_count]);             
            %         if(j==5 & i==5)
            %            ylim([-max_count,max_count]); 
            %         else
            %             ylim([0,max_count]);        
            %         end        
            axis off;
        end
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%---------------------------------------------------------------------------------------
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s'];
for i = 1 : x_length * 3
     sprint_txt = [sprint_txt, ' %1.2f'];    
end
%buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_45(1,:)/8, count_y_45(2,:)/8,count_y_90(1,:)/8, count_y_90(2,:)/8,count_y_135(1,:)/8, count_y_135(2,:)/8, count_y_180(1,:), count_y_180(2,:));  
buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2} );  % for 2 conditions
%buff = sprintf(sprint_txt, FILE, count_y_max{1} );   % for 1 conditions
%buff = sprintf(sprint_txt, FILE, count_trial_beg,count_trial_end ); 

outfile = [BASE_PATH 'ProtocolSpecific\MOOG\CueConflict2D\CueConflict2D_PSTH.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

return;

