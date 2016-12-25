% analysis of eyetrace for heading discrimination
% including eye position, eye velocity
% for staircase, the heading and repetition are unsymmetric, need to be
% dealt with by different method
% -GY
function HeadingDis_cum_eyetrace_staircase(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 YG

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);
temp_spike_data = data.spike_data(1,:);   % spike rasters

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
total_trials = temp_total_trials( select_trials);
azimuth = temp_azimuth(select_trials);
elevation = temp_elevation(select_trials);
stim_type = temp_stim_type(select_trials);
amplitude = temp_amplitude(select_trials);
heading = temp_heading( select_trials );
spike_rates = temp_spike_rates(select_trials);
spike_rates = 1:length(spike_rates); % creat fake data for psychophisical to avoid 0 spikes/s instead
eye_data = data.eye_data(:,:,select_trials);

%spike_rates(1:length(spike_rates)) = rand(1, length(spike_rates)); % for staircase and psychophysical

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_heading = munique(heading');

one_repetition = length(unique_heading)*length(unique_stim_type);
repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate eye velocity based on trials to replace spike_rate and then
% calculate neuronal threshold and CP, similar to the analysis of accelerometer
if sum( eye_data(3,201:300,1) ) ~= 0    % take the left coil if it exist
    eye_chan = 3;
else
    eye_chan = 1;
end
eye_chan = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus
for i = 1 : length(unique_heading)    
    select = find( heading==unique_heading(i) );  % group left headings
    for n = 1 : length(select)  
        DC_hor = mean( eye_data(eye_chan,201:250,select(n)) );
        eye_pos_hor_trial{i}(n,:) = eye_data(eye_chan,201:1000,select(n)) - DC_hor; 
    end  
end

% monkey's choice
%determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
LEFT = 1;
RIGHT = 2;
for i= 1 : length(total_trials) 
    temp = data.event_data(1,:,i + BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = RIGHT;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LEFT;
    else
     %   choice(i) = RIGHT;
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
end
select_left = find( choice==LEFT );  % group left choices
select_right = find( choice==RIGHT );  % group right choices
for l = 1:length(select_left)
    DC_hor = mean( eye_data(eye_chan,201:250,select_left(l)) );
    eye_pos_hor_trial_leftchoice(l,:) = eye_data(eye_chan,201:1000,select_left(l)) - DC_hor; 
end
for r = 1:length(select_right)
    DC_hor = mean( eye_data(eye_chan,201:250,select_right(r)) );
    eye_pos_hor_trial_rightchoice(r,:) = eye_data(eye_chan,201:1000,select_right(r)) - DC_hor; 
end
%%% plot now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);  % plot according to headings
set(2,'Position', [5,1 1000,740], 'Name', 'Heading Discrimination');
orient landscape;
for i = 1 : length(unique_heading)    
    if unique_heading(i) <0   % left headings
        subplot(1,2,1)
        plot( eye_pos_hor_trial{i}(:,:)','b-' ); 
        hold on;
    else  % right headings
        subplot(1,2,2)
        plot( eye_pos_hor_trial{i}(:,:)','r-' ); 
        hold on;
    end
    xlim([0, 800]); % should be a reasonable rangle
%    ylim([-10,10]);   % saccade window
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
end

figure(3);  % plot according to monkey's choices
set(3,'Position', [5,1 1000,740], 'Name', 'Heading Discrimination');
orient landscape;

subplot(1,2,1)
plot( eye_pos_hor_trial_leftchoice(:,:)','b-' ); 

subplot(1,2,2)
plot( eye_pos_hor_trial_rightchoice(:,:)','r-' ); 

xlim([0, 800]); % should be a reasonable rangle
%ylim([-10,10]);   % saccade window
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Ayanna will add something here -ASB(31/08/2007)
% output some text of basic parameters in the figure
figure(4);  % plot according to monkey's choices
set(3,'Position', [5,1 1000,740], 'Name', 'Heading Discrimination');
orient landscape;

targl = input('what is the LEFT target eccentricity?');
subplot(1,2,1)
plot( eye_pos_hor_trial_leftchoice(:,:)','b-' ); 
hold on
plot([500 700],[targl targl],'k')
hold off
xlabel('Leftward saccades')
xlim([500, 700]);

targr = -(targl);
subplot(1,2,2)
plot( eye_pos_hor_trial_rightchoice(:,:)','r-' ); 
hold on
plot([500 700],[targr targr],'k')
hold off
xlabel('Rightward saccades')
xlim([500, 700]); % should be a reasonable rangle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %output to text file
% sprint_txt = ['%s'];
% for i = 1 : 100
%      sprint_txt = [sprint_txt, ' %1.3f'];    
% end
% buff = sprintf(sprint_txt, FILE, eye_pos_hor_head(:,:)', eye_pos_ver_head(:,:)', eye_pos_hor_left(:,:)',eye_pos_hor_right(:,:)' );
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\Eye_velocity_raw_staircase.dat'];
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

