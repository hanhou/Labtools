%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function PURSUIT_HEADING_eyemovement(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_total_trials = data.misc_params(OUTCOME, :);
temp_mask_status = data.moog_params(MASK_STATUS,:,MOOG);
temp_mask_radius = data.moog_params(MASK_RADIUS,:,MOOG);
temp_microstim = data.misc_params(MICROSTIM,:);
temp_fix_x = data.moog_params(FIX_X,:,MOOG); 
pursuit_speed = data.moog_params(PURSUIT_VELOCITY,1,MOOG)

trials = 1:length(temp_heading);		% a vector of trial indices
select_trials =  (trials >= BegTrial) & (trials <= EndTrial) ;
%select_trials =  ((trials >= BegTrial) & (trials<=540)) | ((trials >= 1081) & (trials <= EndTrial));

stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
total_trials = temp_total_trials( select_trials);
mask_status= temp_mask_status( select_trials );
mask_radius= temp_mask_radius( select_trials );
microstim = temp_microstim(select_trials );
fix_x = temp_fix_x( select_trials);

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
unique_mask_status = munique(mask_status');
unique_mask_radius = munique(mask_radius');
unique_microstim = munique(microstim');
unique_fix_x = munique(fix_x');

correct_rate = [];
yy=[];
figure;
for i=1:1080
    plot(data.eye_data(3,100:800,i),'-');
    title(num2str(i));
end
for n = 1:length(unique_fix_x)
    for k = 1:length(unique_stim_type)   
        select = find( (stim_type==unique_stim_type(k)) & (fix_x == unique_fix_x(n)) ) ;
  %      offset_x = mean( data.eye_data(1,180:200,select) ); % horizontal
  %      offset_y = mean( data.eye_data(4,180:200,select) ); % vertical
        for jj = 1:length(select)
            resp_x{k,n}(jj,:) = data.eye_data(3,101:700,select(jj));  % horizontal   
            resp_y{k,n}(jj,:) = data.eye_data(4,101:700,select(jj));  %
        end
    end
end 
figure;
for i=1:200
    plot(data.eye_data(7,101:700,i),'-');
 %   hold on;    
end
%---------------------------------------------------------------------------------------
return;