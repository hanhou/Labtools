%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------
function Direction2d_cue_conflict_basic_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth_moog = data.moog_params(HEADING,:,MOOG);
temp_azimuth_cam = data.moog_params(HEADING,:,CAMERAS);
temp_preferred_azimuth = data.moog_params(PREFERRED_AZIMUTH,:,MOOG);
temp_preferred_elevation = data.moog_params(PREFERRED_ELEVATION,:,CAMERAS);
preferred_azimuth = data.one_time_params(PREFERRED_AZIMUTH);
preferred_elevation = data.one_time_params(PREFERRED_ELEVATION);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth_moog == data.one_time_params(NULL_VALUE)) & (temp_azimuth_cam == data.one_time_params(NULL_VALUE)) );

% %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth_moog);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );

azimuth_moog = temp_azimuth_moog(~null_trials & select_trials);
azimuth_cam = temp_azimuth_cam(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);

spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth_moog = munique(azimuth_moog');
unique_azimuth_cam = munique(azimuth_cam');
% unique_elevation = munique(elevation');

% Grab coherence value for this experiment.
coherence = data.moog_params(COHERENCE,1,CAMERAS);
repetition = length(spike_rates)/(length(unique_azimuth_moog)*length(unique_azimuth_cam)-length(find(null_trials==1)));
%% ADD CODE HERE FOR PLOTTING
% create basic matrix represents each response vector    
resp = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_std = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_ste = zeros( length(unique_azimuth_moog) , length(unique_azimuth_cam) );
resp_trial = zeros( ...
    sum(select_trials) / ( length(unique_azimuth_moog)*length(unique_azimuth_cam) ), ...
    length(unique_azimuth_moog) , length(unique_azimuth_cam) );
for i=1:length(unique_azimuth_moog)
    for j=1:length(unique_azimuth_cam)
        select = logical( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) );
        if (sum(select) > 0) 
            resp_trial(1:sum(select),i,j) = spike_rates(select);
            resp(i, j) = mean(spike_rates(select));
            resp_std(i,j) = std(spike_rates(select));
            resp_std_root(i,j) = std(sqrt(spike_rates(select)));  
%             resp_ste(i,j) = resp_std(i,j) / sqrt(length(find( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) )) );
            resp_ste(i,j) = resp_std(i,j) / sqrt( sum(select) );
        end
    end
end 

% the first one is the conflict dataset
resp_conflict = resp(2:length(unique_azimuth_moog), 2:length(unique_azimuth_cam) );
resp_trial_conflict = resp_trial( :, 2:length(unique_azimuth_moog), 2:length(unique_azimuth_cam) );

% the second is the moog control
resp_cam = resp( 1, 2 : length(unique_azimuth_cam) );
resp_cam_ste = resp_ste( 1, 2 : length(unique_azimuth_cam) );
resp_cam_std_root = resp_std_root( 1, 2 : length(unique_azimuth_cam) );
resp_cam_std = resp_std( 1, 2 : length(unique_azimuth_cam) );
resp_trial_cam = squeeze( resp_trial( :, 1, 2 : length(unique_azimuth_cam) ) );
p_1D(2) = anova1(resp_trial_cam,'','off');
dim= size(resp_trial_cam);
repetition = dim(1);

% the third is the camera control
resp_ves = resp( 2:length(unique_azimuth_moog) , 1 );
resp_ves_ste = resp_ste( 2:length(unique_azimuth_moog) , 1 );
resp_ves_std = resp_std( 2:length(unique_azimuth_moog) , 1 );
resp_trial_ves = squeeze( resp_trial( :, 2:length(unique_azimuth_moog) , 1 ) );
p_1D(1) = anova1(resp_trial_ves,'','off');
%%

% Calcuate a few very basic measures such as modulation depth

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

% % Extract single cue tuning by averaging across the other motion cue.
% % Vestibular
% resp_conflict_ves = mean( resp_conflict,2);
% % Visual
% resp_conflict_vis = mean(resp_conflict',2);
% 
% % Calculate center of mass.
% [j1,j2]=pol2cart( unique_azimuth_cam(2:end)*pi/180 , resp_cam' );
% com_vis=cart2pol( sum(j1), sum(j2) ) * 180/pi;
% 
% [j1,j2]=pol2cart( unique_azimuth_moog(2:end)*pi/180 , resp_ves );
% com_ves=cart2pol( sum(j1), sum(j2) ) * 180/pi;
% 
% [j1,j2]=pol2cart( unique_azimuth_cam(2:end)*pi/180 , resp_conflict_vis );
% com_conflict_vis=cart2pol( sum(j1), sum(j2) ) * 180/pi;
% 
% [j1,j2]=pol2cart( unique_azimuth_moog(2:end)*pi/180 , resp_conflict_ves );
% com_conflict_ves=cart2pol( sum(j1), sum(j2) ) * 180/pi;
% 
% % Extract the peak firing rates and their directions.
% [max_ves,j2]=max( resp_ves );
% thetap_ves = unique_azimuth_moog( j2 + 1 );
% 
% [max_vis,j2]=max( resp_cam );
% thetap_vis = unique_azimuth_cam( j2 + 1 );
% 
% [max_conflict_ves,j2]=max( resp_conflict_ves );
% thetap_conflict_ves = unique_azimuth_moog( j2 + 1 );
% 
% [max_conflict_vis,j2]=max( resp_conflict_vis );
% thetap_conflict_vis = unique_azimuth_cam( j2 + 1 );
% 
% % Extract the trough firing rates and their directions.
% [min_ves,j2]=min( resp_ves );
% thetat_ves = unique_azimuth_moog( j2 + 1 );
% 
% [min_vis,j2]=min( resp_cam );
% thetat_vis = unique_azimuth_cam( j2 + 1 );
% 
% [min_conflict_ves,j2]=min( resp_conflict_ves );
% thetat_conflict_ves = unique_azimuth_moog( j2 + 1 );
% 
% [min_conflict_vis,j2]=min( resp_conflict_vis );
% thetat_conflict_vis = unique_azimuth_cam( j2 + 1 );
% 
% % Calculate tuning depths.
% depth_ves = max_ves - min_ves;
% depth_vis = max_vis - min_vis;
% depth_conflict_ves = max_conflict_ves - min_conflict_ves;
% depth_conflict_vis = max_conflict_vis - min_conflict_vis;

eye_x_left_temp(:,:) = data.eye_data(1,:,~null_trials & select_trials);
eye_y_left_temp(:,:) = data.eye_data(2,:,~null_trials & select_trials);
eye_x_right_temp(:,:) = data.eye_data(3,:,~null_trials & select_trials);
eye_y_right_temp(:,:) = data.eye_data(4,:,~null_trials & select_trials);
dim1 = size(eye_x_left_temp);
for i=1:dim1(2)
    eyeleftmaxmin(i) = abs( max(eye_x_left_temp(1:600,i))-min(eye_x_left_temp(1:600,i)) );
    eyerightmaxmin(i) = abs( max(eye_x_right_temp(1:600,i))-min(eye_x_right_temp(1:600,i)) );
end
for i=1:200
    eye_x_left(i,:) = eye_x_left_temp(i+322,:)-mean(eye_x_left_temp(222:322,:));
    eye_y_left(i,:) = eye_y_left_temp(i+322,:)-mean(eye_y_left_temp(222:322,:));
    eye_x_right(i,:) = eye_x_right_temp(i+322,:)-mean(eye_x_right_temp(222:322,:));
    eye_y_right(i,:) = eye_y_right_temp(i+322,:)-mean(eye_y_right_temp(222:322,:));
    deviation_left(i,:) = sqrt(eye_x_left(i,:).^2+eye_y_left(i,:).^2);
    deviation_right(i,:) = sqrt(eye_x_right(i,:).^2+eye_y_right(i,:).^2);
end
if median(eyeleftmaxmin)<median(eyerightmaxmin) % left eye has no signal
    deviation(1,:) = median(deviation_right(:,:),2);
    eye_pos = deviation_right(:,:);
else
    deviation(1,:) = median(deviation_left(:,:),2);
    eye_pos = deviation_left(:,:);
end
dim = size(eye_pos);
count = 0;
for j=1:dim(2)    
    gauss2 = normpdf(1:1:200,100,1); % smooth data at a SD of 5ms (1 point)
    eye_temp = conv(eye_pos(:,j), gauss2);
    eye_pos_smooth(:,j) = eye_temp(100:end-ceil(200/2)); 
    eye_vel(:,j) = abs(diff(eye_pos_smooth(:,j)))*200; 
    saccadefind = find( eye_vel(:,j)>10 );
    if length(saccadefind)>=1 
       tempfind = find(diff(saccadefind)~=1);
       if length(tempfind)>=1
           saccade_count = length(tempfind)+1; % the number of microscaddes            
           tempfind2 = [0 tempfind' length(saccadefind)];
           for n=1:saccade_count
               peak_vel(n+count) = max( eye_vel(saccadefind(tempfind2(n)+1):saccadefind(tempfind2(n+1)),j) );
           end
       else
           saccade_count = 1;           
           peak_vel(1+count) = max(eye_vel(:,j));
       end
    else
        saccade_count = 0;        
    end   
    count = count+saccade_count; 
end
saccade_rate=count/dim(2);
saccade_vel = hist(log10(peak_vel), 1:0.05:2)/count;

% Write important values to a file.
sprint_txt = ['%s\t'];
for i = 1 : 1000
     sprint_txt = [sprint_txt, ' %1.3f\t'];    
end
% buff = sprintf(sprint_txt, FILE, coherence, spon_resp, p_1D, resp_cam,resp_cam_ste,resp_cam_std_root );
% outfile = ['Z:\Users\Yong\Paper_coherence\LowHighCoherenceVisual.dat'];
% buff = sprintf(sprint_txt, FILE,8,repetition, squeeze(resp_trial(1:repetition,2:9,1)),squeeze(resp_trial(1:repetition,1,2:9)) );
% buff = sprintf(sprint_txt,FILE,spon_resp); 
% outfile = ['Z:\Users\Yong\noisecorrelation882_spontaneous.dat'];

buff = sprintf(sprint_txt, FILE, saccade_rate, deviation, saccade_vel);
outfile = ['Z:\Users\Yong\eyetrain.dat'];

printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
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