%-----------------------------------------------------------------------------------------------------------------------
%-- Train a two layer network with basis functions specified by
%--   tuning model and parameters using the delta rule:
%--     W(t+1) = W(t) + learning_rate * FR_real * (FR_ideal - FR_real)
%-- Created - pwatkins, 6/04
%-----------------------------------------------------------------------------------------------------------------------
function [training_inputs, training_outputs] = make_training_set(fun,fun_eye,Ptvis,Ptves,P_eye,...
    range_az,range_el,range_gaze_az,range_gaze_el)

Curvefit_defines;

% extract number of parameters of our model and number of neurons in the
% basis layer from the basis layer parameters.
[num_neurons num_params] = size(Ptvis);
[num_eye_neurons num_eye_params] = size(P_eye);

create_sample_points;

% define a minimum response:
% for the visual inputs to in vestibular only case 
min_rsp_vis = -1;
% for the vestibular inputs to in the visual only case.
min_rsp_ves = -1;

% create a training set for all three stimulus types
% CURVEFIT_OCCULAR_STIM = 2;
% CURVEFIT_VESTIBULAR_STIM = 1;
% CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM = 3;
num_stim_types = 3;

% create the training set
num_training_trials = num_stim_types*num_unique_points*num_gaze_angles;
training_inputs = zeros(2*num_neurons + 2*num_eye_neurons, num_training_trials);
training_outputs = zeros(num_neurons, num_training_trials);
for m=1:num_stim_types
            
    for j=1:num_gaze_angles
        
        % this is idiotic, not sure why i did this...
        %fix_az = pi/2 + unique_gaze_angles_xr(j);   % because horizontal gaze angles vary from straight ahead
        fix_az = unique_gaze_angles_xr(j);  % doing?
        fix_el = unique_gaze_angles_yr(j);
        
        % if you want eye-centered output, run following
        % gaze angles vary hc coords
%         [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(unique_point_azimuth_r,...
%             unique_point_elevation_r,-fix_az,-fix_el);
%         az_hc = azimuth_gaze; el_hc = elevation_gaze;
%         az_ec = unique_point_azimuth_r; el_ec = unique_point_elevation_r;
        % if you want head-centered output, run following
        % gaze angles vary ec coords
        [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(unique_point_azimuth_r,...
            unique_point_elevation_r,fix_az,fix_el);
        az_hc = unique_point_azimuth_r; el_hc = unique_point_elevation_r;
        az_ec = azimuth_gaze; el_ec = elevation_gaze;
        
        % start and stop indices for this trial set
        beg_i = (m-1)*num_gaze_angles*num_unique_points + (j-1)*num_unique_points + 1;
        end_i = (m-1)*num_gaze_angles*num_unique_points + j*num_unique_points;
        
        % compute eye centered inputs and ideal head centered response
        for k=1:num_neurons
            Ptmpvis = Ptvis(k,1:num_params);
            
            % input layer does not vary amplitude as a function of gaze angle
            Ptmpvis(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE)) = ...
                Ptmpvis(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_OFFSET));
            
            % input layer does not vary dc offset as a function of gaze angle
            Ptmpvis(curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET)) = ...
                Ptmpvis(curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSETB));

            Ptmpves = Ptves(k,1:num_params);
            
            % input layer does not vary amplitude as a function of gaze angle
            Ptmpves(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE)) = ...
                Ptmpves(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_OFFSET));
            
            % input layer does not vary dc offset as a function of gaze angle
            Ptmpves(curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET)) = ...
                Ptmpves(curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSETB));
                        
            % compute our fit data with the given parameters
            if m == CURVEFIT_VESTIBULAR_STIM
                %training_inputs(k,beg_i:end_i) = zeros(1,num_unique_points);
                training_inputs(k,beg_i:end_i) = min_rsp_vis*ones(1,num_unique_points);
            else
                training_inputs(k,beg_i:end_i) = feval( fun, Ptmpvis, az_ec, el_ec );
                %training_inputs(k,beg_i:end_i) = feval( fun, Ptmpvis, az_hc, el_hc );
            end
            
            if m == CURVEFIT_VESTIBULAR_STIM
                % compute our output fit data in hc coords (or ec coords)
                % if you want eye-centered output, run following
                %training_outputs(k,beg_i:end_i) = feval( fun, Ptmpvis, az_ec, el_ec );
                % if you want head-centered output, run following
                training_outputs(k,beg_i:end_i) = feval( fun, Ptmpves, az_hc, el_hc );
            else
                % compute our output fit data in hc coords (or ec coords)
                % if you want eye-centered output, run following
                %training_outputs(k,beg_i:end_i) = feval( fun, Ptmpvis, az_ec, el_ec );
                % if you want head-centered output, run following
                training_outputs(k,beg_i:end_i) = feval( fun, Ptmpvis, az_hc, el_hc );
            end
            
            % add vestibular inputs that are head centered
            if m == CURVEFIT_OCCULAR_STIM
                %training_inputs(num_neurons+2*num_eye_neurons+k,beg_i:end_i) = ...
                %    zeros(1,num_unique_points);
                training_inputs(num_neurons+2*num_eye_neurons+k,beg_i:end_i) = ...
                    min_rsp_ves*ones(1,num_unique_points);
            else
                training_inputs(num_neurons+2*num_eye_neurons+k,beg_i:end_i) = ...
                    feval( fun, Ptmpves, az_hc, el_hc );
                %training_inputs(num_neurons+2*num_eye_neurons+k,beg_i:end_i) = ...
                %    feval( fun, Ptmpves, az_ec, el_ec );
            end        
        end            
        
        % compute horizontal gaze angle inputs
        for k=1:num_eye_neurons
            %%if m == CURVEFIT_VESTIBULAR_STIM
            %if m == CURVEFIT_OCCULAR_STIM
            if 0
            %if 1
                training_inputs(num_neurons+k,beg_i:end_i) = zeros(1,num_unique_points);
                %training_inputs(num_neurons+k,beg_i:end_i) = normrnd(15,5,1,num_unique_points);
            else
                training_inputs(num_neurons+k,beg_i:end_i) = ...
                    repmat(feval( fun_eye, P_eye(k,1:num_eye_params), unique_gaze_angles_xr(j) ), ...
                    1, num_unique_points);
            end
        end
        
        % compute vertical gaze angle inputs
        for k=1:num_eye_neurons
            %%if m == CURVEFIT_VESTIBULAR_STIM
            %if m == CURVEFIT_OCCULAR_STIM
            if 0
            %if 1
                training_inputs(num_neurons+num_eye_neurons+k,beg_i:end_i) = zeros(1,num_unique_points);
                %training_inputs(num_neurons+num_eye_neurons+k,beg_i:end_i) = normrnd(15,5,1,num_unique_points);
            else
                training_inputs(num_neurons+num_eye_neurons+k,beg_i:end_i) = ...
                    repmat(feval( fun_eye, P_eye(k,1:num_eye_params), unique_gaze_angles_yr(j) ), ...
                    1, num_unique_points);
            end
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

