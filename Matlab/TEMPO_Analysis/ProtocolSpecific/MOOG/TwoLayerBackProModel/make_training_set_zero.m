%-----------------------------------------------------------------------------------------------------------------------
%-- Train a two layer network with basis functions specified by
%--   tuning model and parameters using the delta rule:
%--     W(t+1) = W(t) + learning_rate * FR_real * (FR_ideal - FR_real)
%-- Created - pwatkins, 6/04
%-----------------------------------------------------------------------------------------------------------------------
function [training_inputs] = make_training_set(fun,fun_eye,Pt,P_eye,...
    range_az,range_el,range_gaze_az,range_gaze_el)

Curvefit_defines;

% extract number of parameters of our model and number of neurons in the
% basis layer from the basis layer parameters.
[num_neurons num_params] = size(Pt);
[num_eye_neurons num_eye_params] = size(P_eye);

create_sample_points;

% create a training set for all three stimulus types
% CURVEFIT_OCCULAR_STIM = 2;
% CURVEFIT_VESTIBULAR_STIM = 1;
% CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM = 3;
num_stim_types = 3;

% create the training set
num_training_trials = num_gaze_angles;
training_inputs = zeros(2*num_neurons + 2*num_eye_neurons, num_training_trials);
    for j=1:num_gaze_angles
        
        % this is idiotic, not sure why i did this...
        %fix_az = pi/2 + unique_gaze_angles_xr(j);   % because horizontal gaze angles vary from straight ahead
        fix_az = unique_gaze_angles_xr(j);  % doing?
        fix_el = unique_gaze_angles_yr(j);
                
        % compute eye centered inputs and ideal head centered response
        for k=1:num_neurons
            % set visual inputs to zero
            training_inputs(k,j) = 0;

            % set vestibular inputs to zero
            training_inputs(num_neurons+2*num_eye_neurons+k,j) = 0;
        end            
        
        % compute horizontal gaze angle inputs
        for k=1:num_eye_neurons
            training_inputs(num_neurons+k,j) = ...
                feval( fun_eye, P_eye(k,1:num_eye_params), unique_gaze_angles_xr(j) );
        end
        
        % compute vertical gaze angle inputs
        for k=1:num_eye_neurons
            training_inputs(num_neurons+num_eye_neurons+k,j) = ...
                feval( fun_eye, P_eye(k,1:num_eye_params), unique_gaze_angles_yr(j) );
        end
        
        %         % for testing training
        %         for i = beg_i:end_i
        %             ec_inputs = training_inputs(1:num_neurons,i);
        %             hc_outputs = training_outputs(1:num_neurons,i);
        %             basis_model_contour_plots(unique_azimuth, unique_elevation, ...
        %                 unique_point_azimuth, unique_point_elevation, ec_inputs, hc_outputs);
        %             %         sprintf('%e FOE_az %.2f FOE_el %.2f gaze_az %.2f gaze_el %.2f', SSE(t), FOE_az/pi*180, ...
        %             %             FOE_el/pi*180, gaze_az/pi*180, gaze_el/pi*180)
        %             pause;        
        %         end
            
    end % for each sampled gaze angle
end % for each stim type

% % for testing training
% for i=1:num_unique_points
%     %j = repmat( [repmat([1:26],1,num_gaze_angles/2) zeros(1,num_unique_points*num_gaze_angles/2)], 1, num_stim_types );
%     j = repmat( [zeros(1,num_unique_points*num_gaze_angles/2) repmat([1:26],1,num_gaze_angles/2)], 1, num_stim_types );
%     %gaze_inputs = training_inputs(num_neurons+1:num_neurons+num_eye_neurons,j==i);
%     gaze_inputs = training_inputs(num_neurons+num_eye_neurons+1:num_neurons+2*num_eye_neurons,j==i);
%     for m=1:num_stim_types
%         figure(888);
%         plot(gas/pi*180, gaze_inputs(1:num_gaze_angles/2,(m-1)*num_gaze_angles/2+1:m*num_gaze_angles/2));
%     end
%     pause;
% end

