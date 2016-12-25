%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
% PSTH
%--	01/06 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cum_PSTH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);
% SpikeChan = 1;

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME, :);
temp_spike_data = data.spike_data(SpikeChan,:);   
temp_spike_rates = data.spike_rates(SpikeChan, :); 
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manually omit trials, i.e. due to bumps, lost isolation, etc. (CRF 8-2009)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% omit_trials = [537:563];
% select_trials(omit_trials) = 0;
%
% sumSpikeRates_equalzero = sum(temp_spike_rates<0.1)
% sumSpikeRates_lessthantwo = sum(temp_spike_rates<2)
% maxSpikeRate = max(temp_spike_rates)
% edges = [0 1 2 4 8 16 100];
% histc(temp_spike_rates,edges)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
total_trials = temp_total_trials( select_trials);
spike_rates = temp_spike_rates( select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');

h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';

% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );
spike_data( find(spike_data>100) ) = 1; % something is absolutely wrong 

% monkey's choice
LEFT = 1;
RIGHT = 2;
for i= 1 : length(spike_rates) 
    temp = data.event_data(1,:,i + BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LEFT;
    else
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
end
% if FILE=='m2c384r2.htb'
%    choice(889) =2; % for cell m2c384r2 % for some reason the choice is 0 for
% end

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
time_step_left=1;
time_step_right=1;
for k=1: length(unique_stim_type)
    lefttemp =find( (heading == unique_heading(5)) & (stim_type == unique_stim_type(k)) & choice==1 ) ;
    righttemp =find( (heading == unique_heading(5)) & (stim_type == unique_stim_type(k)) & choice==2 ) ;
    for i=1:length(unique_heading)
        select = logical( (heading==unique_heading(i)) & (stim_type==unique_stim_type(k)) );  
        act_found = find( select==1 );
        % count spikes per timebin on every same condition trials
        for repeat=1:length(act_found) 
            for n=1:(x_length)
                temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                time_step=time_step+timebin;
            end
            time_step=1;                    
        end
        count_y_trial{i,k}(:,:) = temp_count;  % each trial's PSTH 
     
        % get the average of the total same conditions if repetion is > 1
        dim=size(temp_count);
        count_y{i,k} = mean(temp_count);
        max_count_y(i,k) = max(count_y{i,k});
    end 
    
%     for repeat_left=1:length(lefttemp) 
%         for n=1:(x_length)
%             temp_count_left(repeat_left,n)=sum(spike_data(1,(frequency*(lefttemp(repeat)-1)+time_step):(frequency*(lefttemp(repeat)-1)+n*timebin)));
%             time_step_left=time_step_left+timebin;
%         end
%         time_step_left=1;                    
%     end
%     count_y_left(k,:) = mean(temp_count_left);
%     
%     for repeat_right=1:length(righttemp) 
%         for n=1:(x_length)
%             temp_count_right(repeat_right,n)=sum(spike_data(1,(frequency*(righttemp(repeat)-1)+time_step):(frequency*(righttemp(repeat)-1)+n*timebin)));
%             time_step_right=time_step_right+timebin;
%         end
%         time_step_right=1;                    
%     end
%     count_y_right(k,:) = mean(temp_count_right);
    % normalize PSTH to 1 for each stimulus condition
%     for i=1:length(unique_heading) 
%         count_y{i,k} = count_y{i,k} / max(max_count_y(:,k));
%     end
end

% this part find which heading is the maximum response or minimum response
for k=1: length(unique_stim_type)
    for i = 1 : length(unique_heading)
        ss(i) = sum(count_y{i,k}(x_length*3/10:x_length*5.5/10)); % only use 2 middle second data
    end
    mm=find( ss==max(ss)); 
    nn=find( ss==min(ss));
    max_index(k) = mm(1); 
    min_index(k) = nn(1);
end

% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max(max(max_count_y))];
% define figure
figure(2);
set(2,'Position', [5,5 1000,680], 'Name', 'Tuning');
orient portrait; %changed from landscape by asb (30 july 2007)
axis off;

xoffset=0;
yoffset=0;

% now plot
for k=1: length(unique_stim_type) 
    
    axes('position',[0 0 1 1]); 
    xlim([-50,50]);
    ylim([-50,50]);
    % here starts the column identification
    text(-25,45, 'vestibular'); 
    text(0,45, 'visual'); 
    text(25,45, 'combined'); 
    text(-45, 45, [FILE ', SpChan ' num2str(SpikeChan)]);
    if k == 1 
        for j = 1:length(unique_heading)
            text(-48, 45-j*8, num2str(unique_heading(j)) );
        end
    end 
    axis off;
    hold on;
    
    for i=1:length(unique_heading)        
        axes('position',[0.31*(k-1)+0.1 (0.92-0.08*i) 0.25 0.05]); %this changes the size and location of each row of figures.
         
        plot( x_time,count_y{i,k}(1,:) );  
        hold on;
        plot( x_start, y_marker, 'r-');
        plot( x_stop,  y_marker, 'r-');
        set( gca, 'xticklabel', ' ' );
        % set the same scale for all plot
        xlim([0,x_length]);
        ylim([max(max(max_count_y))*0.75,max(max(max_count_y))]);
    end 
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s']; 
% for i = 1 : 3000 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f'];    
% end
% buff = sprintf(sprint_txt, FILE, count_y{max_index(1),1},count_y{max_index(2),2} ,count_y{max_index(3),3});
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cumPSTH5s.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)   % file does not yet exist
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