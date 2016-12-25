 %----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Katsu, 06/04/07, Need Vest and Visual condition, do not run only
%vestibular condition
%-----------------------------------------------------------------------------------------------------------------------

function Rotation_PSTH_output(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Analysis);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

% % temp: save out .mat files for Tanya -- 12-2009
% save(['Z:\Users\Tanya\Nodulus_frequency_analysis\all_data\' FILE '.mat']); 
% save(['C:\Nodulus_analysis\all_data\' FILE '.mat']);


%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_spike_data = data.spike_data(SpikeChan,:);
temp_spike_rates = data.spike_rates(SpikeChan, :);    

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
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

% find spontaneous trials which azimuth,elevation,stim_type=-99
spon_found = find(null_trials==1);     

% remove null trials, bad trials (rates >3000), and trials outside Begtrial~Endtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 99;
end
spike_data = temp_spike_data( temp_spike_data~=99 );
spike_data( find(spike_data>100) ) = 1; % something is absolutely wrong 

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
%              if max(count_y{i,j,k})~=0;
%                 count_y_norm{i,j,k}=count_y{i,j,k} / max(count_y{i,j,k});
%              else
%                 count_y_norm{i,j,k}=0;
%              end
        end
    end  
    % now find the peak
    [row_max, col_max] = find( resp{k}(:,:)==max(max(resp{k}(:,:))) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k}=row_max(1);
    col_m{k}=col_max(1);
    if max(count_y{col_max(1), row_max(1), k})~=0;
%        count_y_max{k} = count_y{col_max(1), row_max(1), k} /
%        max(count_y{col_max(1), row_max(1), k});%normalized in each
%        condition
       max_count_stim{k} = max(count_y{col_max(1), row_max(1), k});
       count_y_max{k} = count_y{col_max(1), row_max(1), k};
    else
       count_y_max{k} =0;
    end
    % find the largest y to set scale later
    if max(count_y{col_max(1), row_max(1), k}) > max_count
        max_count = max(count_y{col_max(1), row_max(1), k});
    end
    % % % To do the opposite trace AB October 2007  
% find the minimum values of the responses
    [row_min, col_min] = find( resp{k}(:,:)==min(min(resp{k}(:,:))) );    
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_mi{k}=row_min(1);
    col_mi{k}=col_min(1);
    if min(count_y{col_min(1), row_min(1), k})~=0;
%        count_y_min{k} = count_y{col_min(1), row_min(1), k} /
%        min(count_y{col_min(1), row_min(1), k});%normalized in each
%        condition
       min_count_stim{k} = min(count_y{col_min(1), row_min(1), k});
    else
       count_y_min{k} =0;
    end

end
% 
% min_count_stim{:}
% max_count


%%  normalize by maximum visual stimulus = 1, this is definition
for k=1: length(unique_condition_num)  
    for j=1:length(unique_elevation)       
        for i=1: length(unique_azimuth)
            
%                 count_y_norm{i,j,k}=count_y{i,j,k} / max_count_stim{2};%k=2, visual is normalization
                  count_y_norm{i,j,k}=count_y{i,j,k} / max_count_stim{1};%k=1, vestibular signal is normalization 20\9\2007 AB
               
%                 count_y_max{k} = count_y{col_max(1), row_max(1), k} / max(count_y{col_max(1), row_max(1), k});
%        
        end
    end
%     pickupmax{k}=count_y{col_m{k},row_m{k},k}(21:60) / max_count_stim{2}; %same as below   
    pickupmax{k}=count_y_norm{col_m{k},row_m{k},k}(21:60);  
    pickupmin{k}=count_y_norm{col_mi{k},row_mi{k},k}(21:60); 
    % later, save pickupmax{k}, for normalized plottiong

end


% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
% define figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', '3D Rotation Tuning');
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
    text(-47,-40, 'Azim: 270        225       180       135         90         45         0         315        270');
    text(25,-40, 'Rotation');
    axis off;
    hold on;
    for i=1:length(unique_azimuth)+1                    % aizmuth 270 are plotted two times in order to make circular data
        for j=1:length(unique_elevation)
            axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
            if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
                bar( x_time,count_y_norm{8-i,j,k}(1,:) );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
            elseif(i==8)
                bar( x_time,count_y_norm{i,j,k}(1,:) ); 
            else
                bar( x_time,count_y_norm{7,j,k}(1,:) ); 
            end
            hold on;
            plot( x_start, y_marker, 'r-');
            plot( x_stop,  y_marker, 'r-');
            set( gca, 'xticklabel', ' ' );
            % set the same scale for all plot
            xlim([0,x_length]);
%             ylim([0,max_count/max_count_stim{2}]);%----normalized by visual (k=2) stimulus
            ylim([0,max_count/max_count_stim{1}]);%----normalized by vestibular (k=1) stimulus 20\9\2007 AB
%             ylim([0,10]);% for m3c296r1r3, [0,6], for m3c294r1r3, [0,10]
        end    
    end 

    xoffset=xoffset+0.46;
    
end
%---------------------------------------------------------------------------------------
%Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s'];
% for i = 1 : x_length * 5
%     sprint_txt = [sprint_txt, ' %1.2f'];    
% end
% % buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_45(1,:)/8, count_y_45(2,:)/8,count_y_90(1,:)/8, count_y_90(2,:)/8,count_y_135(1,:)/8, count_y_135(2,:)/8, count_y_180(1,:), count_y_180(2,:));  
% % buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_max{3} );  
% % buff = sprintf(sprint_txt, FILE, count_trial_beg,count_trial_end ); 
% %buff = sprintf(sprint_txt, FILE, pickupmax{:}); 
% buff = sprintf(sprint_txt, FILE, count_y{3, 3, 1},count_y{7, 3, 1}, count_y{1, 3, 1},count_y{5, 3, 1});
% 
% % outfile = ['Z:\Users\Syed Chowdhury\syed2\tempo_output\Peak_PSTH_output_Syed.dat'];
% outfile = ['Z:\Users\Yong\psth_rotationrollpitch.dat'];
% % % outfile = ['Z:\Users\Ayanna\tempo_output\rotation_PSTH_output_Ayanna.dat']; %saves data to my z:\users folder ab 20 Sept. 2007.
% % 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%    printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%    fprintf(fid, 'FILE\t');
%    fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
% %output_specifications  %%Tunde 9/28/07
% 
return;
% % 
