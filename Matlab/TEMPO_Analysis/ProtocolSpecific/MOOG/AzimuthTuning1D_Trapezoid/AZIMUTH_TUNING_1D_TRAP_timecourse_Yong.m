 %----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Yong, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------

function AZIMUTH_TUNING_1D_TRAP_timecourse_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_spike_data = data.spike_data(SpikeChan,:);
 

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 

azimuth = temp_azimuth(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';

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
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;
end

spike_data(1,:) = temp_spike_data( 1, find(temp_spike_data(1,:)~=99) );
spike_data(1, find(spike_data>1) ) = 1; % something is absolutely wrong 

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
for k=1: length(unique_stim_type)
    for i=1: length(unique_azimuth)
        select = logical( (azimuth==unique_azimuth(i)) & (stim_type==unique_stim_type(k)) );            
        act_found = find( select==1 );
        % count spikes per timebin on every same condition trials
        for repeat=1:length(act_found) 
            for n=1:(x_length)
                temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                time_step=time_step+timebin;
            end
            time_step=1;                    
        end
        count_y{k}(i,:) = mean(temp_count);        
    end 
    maxbin_stim(k) = max(count_y{k}(i,:));
end
maxbin=max(maxbin_stim);

% plot PSTH now
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [(StartEventBin(1,1)+115)/timebin, (StartEventBin(1,1)+115)/timebin];
x_stop =  [(StopEventBin(1,1)+115)/timebin,  (StopEventBin(1,1)+115)/timebin];
y_marker=[0,maxbin];
% define figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', '1D Direction Tuning trapezoid');
orient landscape;
title([FILE ' ' 'channel=' num2str(SpikeChan)]);
axis off;

for k=1: length(unique_stim_type) 
    for i=1:length(unique_azimuth) 
        axes('position',[0.05+0.11*(i-1) 0.65-0.3*(k-1) 0.1 0.2]);
                          
        bar( x_time,count_y{k}(i,:) );   

        hold on;
        plot( x_start, y_marker, 'r-');
        plot( x_stop,  y_marker, 'r-');
        set( gca, 'xticklabel', ' ' );
        % set the same scale for all plot
        xlim([0,x_length*2.5/5]);
        ylim([0,maxbin]); 
    end     
end

% %---------------------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s'];
% for i = 1 : x_length * 3
%      sprint_txt = [sprint_txt, ' %1.2f'];    
% end
% %buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2},count_y_45(1,:)/8, count_y_45(2,:)/8,count_y_90(1,:)/8, count_y_90(2,:)/8,count_y_135(1,:)/8, count_y_135(2,:)/8, count_y_180(1,:), count_y_180(2,:));  
% %buff = sprintf(sprint_txt, FILE, count_y_max{1},count_y_max{2} );  % for 2 conditions
% buff = sprintf(sprint_txt, FILE, count_y_max{1} );   % for 1 conditions
% %buff = sprintf(sprint_txt, FILE, count_trial_beg,count_trial_end ); 
% 
% %outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\DirectionTuning3D_PSTH_Tanya.dat'];
% outfile = ['Z:\Users\HuiM\psth.dat'];
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

