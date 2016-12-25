%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------
function Direction2d_cue_conflict_basics(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
resp_trial_cam = squeeze( resp_trial( :, 1, 2 : length(unique_azimuth_cam) ) );

% the third is the camera control
resp_ves = resp( 2:length(unique_azimuth_moog) , 1 );
resp_ves_ste = resp_ste( 2:length(unique_azimuth_moog) , 1 );
resp_trial_ves = squeeze( resp_trial( :, 2:length(unique_azimuth_moog) , 1 ) );
%%

% Calcuate a few very basic measures such as modulation depth

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));

% Extract single cue tuning by averaging across the other motion cue.
% Vestibular
resp_conflict_ves = mean( resp_conflict,2);
% Visual
resp_conflict_vis = mean(resp_conflict',2);

% Calculate center of mass.
[j1,j2]=pol2cart( unique_azimuth_cam(2:end)*pi/180 , resp_cam' );
com_vis=cart2pol( sum(j1), sum(j2) ) * 180/pi;

[j1,j2]=pol2cart( unique_azimuth_moog(2:end)*pi/180 , resp_ves );
com_ves=cart2pol( sum(j1), sum(j2) ) * 180/pi;

[j1,j2]=pol2cart( unique_azimuth_cam(2:end)*pi/180 , resp_conflict_vis );
com_conflict_vis=cart2pol( sum(j1), sum(j2) ) * 180/pi;

[j1,j2]=pol2cart( unique_azimuth_moog(2:end)*pi/180 , resp_conflict_ves );
com_conflict_ves=cart2pol( sum(j1), sum(j2) ) * 180/pi;

% Extract the peak firing rates and their directions.
[max_ves,j2]=max( resp_ves );
thetap_ves = unique_azimuth_moog( j2 + 1 );

[max_vis,j2]=max( resp_cam );
thetap_vis = unique_azimuth_cam( j2 + 1 );

[max_conflict_ves,j2]=max( resp_conflict_ves );
thetap_conflict_ves = unique_azimuth_moog( j2 + 1 );

[max_conflict_vis,j2]=max( resp_conflict_vis );
thetap_conflict_vis = unique_azimuth_cam( j2 + 1 );

% Extract the trough firing rates and their directions.
[min_ves,j2]=min( resp_ves );
thetat_ves = unique_azimuth_moog( j2 + 1 );

[min_vis,j2]=min( resp_cam );
thetat_vis = unique_azimuth_cam( j2 + 1 );

[min_conflict_ves,j2]=min( resp_conflict_ves );
thetat_conflict_ves = unique_azimuth_moog( j2 + 1 );

[min_conflict_vis,j2]=min( resp_conflict_vis );
thetat_conflict_vis = unique_azimuth_cam( j2 + 1 );

% Calculate tuning depths.
depth_ves = max_ves - min_ves;
depth_vis = max_vis - min_vis;
depth_conflict_ves = max_conflict_ves - min_conflict_ves;
depth_conflict_vis = max_conflict_vis - min_conflict_vis;

%%

save([BASE_PATH 'ProtocolSpecific\MOOG\Cueconflict2D\mat\' FILE(1:end-4) '.mat'], ...
    'azimuth_cam','azimuth_moog',...
    'com_conflict_ves','com_conflict_vis','com_ves','com_vis',...
    'depth_conflict_ves','depth_conflict_vis','depth_ves','depth_vis',...
    'max_conflict_ves','max_conflict_vis','max_ves','max_vis',...
    'min_conflict_ves','min_conflict_vis','min_ves','min_vis',...
    'resp',...
    'resp_cam','resp_cam_ste',...
    'resp_conflict',...
    'resp_conflict_ves','resp_conflict_vis',...
    'resp_std','resp_ste',...
    'resp_trial','resp_trial_cam','resp_trial_conflict','resp_trial_ves',...
    'resp_ves','resp_ves_ste',...
    'spon_resp',...
    'thetap_conflict_ves','thetap_conflict_vis','thetap_ves','thetap_vis',...
    'thetat_conflict_ves','thetat_conflict_vis','thetat_ves','thetat_vis',...
    'unique_azimuth_cam','unique_azimuth_moog');

%%

% Write important values to a file.
sprint_txt = ['%s\t'];
for i = 1 : 500 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %g\t'];   
end

outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Cueconflict2D\Basic.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)   % change headings here if diff conditions varied
    fprintf(fid, 'FILE\t coherence\t ');
    fprintf(fid, 'com_ves\t com_vis\t com_conflict_ves\t com_conflict_vis\t ');
    fprintf(fid, 'depth_ves\t depth_vis\t depth_conflict_ves\t depth_conflict_vis\t ');
    fprintf(fid, 'max_ves\t max_vis\t max_conflict_ves\t max_conflict_vis\t ');
    fprintf(fid, 'min_ves\t min_vis\t min_conflict_ves\t min_conflict_vis\t ');
    fprintf(fid, 'spon_resp\t ');
    fprintf(fid, 'thetap_ves\t thetap_vis\t thetap_conflict_ves\t thetap_conflict_vis\t ');
    fprintf(fid, 'thetat_ves\t thetat_vis\t thetat_conflict_ves\t thetat_conflict_vis\t ');
    fprintf(fid, '\r\n');
end

buff = sprintf( sprint_txt, ...
    FILE, coherence, ...
    com_ves, com_vis, com_conflict_ves, com_conflict_vis, ...
    depth_ves, depth_vis, depth_conflict_ves, depth_conflict_vis, ...
    max_ves, max_vis, max_conflict_ves, max_conflict_vis, ...
    min_ves, min_vis, min_conflict_ves, min_conflict_vis, ...
    spon_resp, ...
    thetap_ves, thetap_vis, thetap_conflict_ves, thetap_conflict_vis, ...
    thetat_ves, thetat_vis, thetat_conflict_ves, thetat_conflict_vis ...
    );
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

%%

return;