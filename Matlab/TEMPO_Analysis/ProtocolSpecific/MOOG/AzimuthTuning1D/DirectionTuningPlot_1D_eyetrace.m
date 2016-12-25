% AZIMUTH_TUNING_1D.m -- Plots response as a function of azimuth for
% heading task
%--	GY, 07/12/04
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_1D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_interocular_dist = data.moog_params(INTEROCULAR_DIST,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             
%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);
motion_coherence = temp_motion_coherence(~null_trials & select_trials);
interocular_dist = temp_interocular_dist(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_motion_coherence = munique(motion_coherence');
unique_interocular_dist = munique(interocular_dist');

% descide whether loop is stim_type or disparity (when there is vary of
% disparity, it's the visual only condition
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_azimuth)
        select = find( (stim_type==unique_stim_type(k)) & (azimuth==unique_azimuth(i)) );
        select2 = find( (stim_type==unique_stim_type(1)) & (azimuth==unique_azimuth(1)) );
%         for n = 1 : length(select)  
%             DC_begin = mean( (data.eye_data(1,201:300,select(n))-data.eye_data(3,201:300,select(n))) );
%             DC_end = mean( (data.eye_data(1,501:600,select(n))-data.eye_data(3,501:600,select(n))) );
%             DC = mean([DC_begin, DC_end]);
%             verg_azi(n) = mean( (data.eye_data(1,301:500,select(n))-data.eye_data(3,301:500,select(n))) );
%             verg_azi(n) = verg_azi(n) - DC; % subtract DC offset
%             verg_trace_trial(n,:) = data.eye_data(1,201:600,select(n))-data.eye_data(3,201:600,select(n));
%             verg_trace_trial(n,:) = verg_trace_trial(n,:) - DC;            
%         end
        for nn = 1: length(select2)
            Left_x(nn,:) = data.eye_data(1,201:600,select2(nn)) - mean(data.eye_data(1,201:600,select2(nn)));
            Left_y(nn,:) = data.eye_data(2,201:600,select2(nn)) - mean(data.eye_data(2,201:600,select2(nn)));
            Right_x(nn,:) = data.eye_data(3,201:600,select2(nn)) - mean(data.eye_data(3,201:600,select2(nn)));
            Right_y(nn,:) = data.eye_data(4,201:600,select2(nn)) - mean(data.eye_data(4,201:600,select2(nn)));
        end
%         % devide data into left, right, forward and backward sections
%         verg_mean{k}(i) = median(verg_azi);
%         verg_trace{k}(i,:) = median( verg_trace_trial(:,:) ) ;
    end
%     if length(unique_azimuth)==10 % plus 67.5 and 112.5 directions
%         verg_trace_exp{k}(1,:) = median(verg_trace{k}(3:5,:));
%         verg_trace_con{k}(1,:) = verg_trace{k}(9,:);
%         verg_trace_0{k}(1,:) = verg_trace{k}(1,:);
%         verg_trace_180{k}(1,:) = verg_trace{k}(7,:);
%     else
%         verg_trace_exp{k}(1,:) = verg_trace{k}(3,:);
%         verg_trace_con{k}(1,:) = verg_trace{k}(7,:);
%         verg_trace_0{k}(1,:) = verg_trace{k}(1,:);
%         verg_trace_180{k}(1,:) = verg_trace{k}(5,:);
%     end
end
figure(2);
set(2,'Position', [5,1 1000,740], 'Name', 'Heading Discrimination-Vestibular');
axes('position',[0.1,0.3, 0.45,0.5] );
orient portrait;
subplot(4,1,1);
plot(Left_x','b-');
ylim([-1,1]);
ylabel('Left X');
title(FILE);
subplot(4,1,2);
plot(Left_y','b-');
ylim([-1,1]);
ylabel('Left Y');
subplot(4,1,3);
plot(Right_x','b-');
ylim([-1,1]);
ylabel('R X');
subplot(4,1,4);
plot(Right_y','b-');
ylim([-1,1]);
ylabel('R Y');

% %---------------------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s'];
% for i = 1 : 1300
%      sprint_txt = [sprint_txt, ' %1.2f'];    
% end
% buff = sprintf(sprint_txt, FILE, verg_trace_exp{1}, verg_trace_con{1}, verg_trace_0{1}, verg_trace_180{1}   );
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\DirectionTuning1D_eyetrace.dat'];
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
% %---------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
return;