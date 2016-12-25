%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_eye_centered.m -- This function computes the error over three
%--   gaze angles for a eye centered model using the specified function
%--   as the model to fit each gaze angle individually.
%-- Created - pwatkins, 4/04
%-----------------------------------------------------------------------------------------------------------------------
function F = Curvefit_eye_centered(xp,fun,azimuth,elevation,trial_data, ...
    num_complete_repititions,gaze_angle,fixation_type)

Curvefit_defines;
% remap the parameters from the default so that non-linearities are
% constant over all gaze angles.
curvefit_parameter_mapping(CURVEFIT_PARAM_AZIMUTH) = 1;
curvefit_parameter_mapping(CURVEFIT_PARAM_ELEVATION) = 2;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN) = 3;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN_AZ) = 3;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN180) = 4;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN_EL) = 4;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE) = 5;
curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET) = 6;
curvefit_parameter_mapping(CURVEFIT_PARAM_WEIGHT180) = 7;

% The eye centered model varies the rotation angles for each of
% the gaze angles, corresponding to the actual gaze angle and direction.
% Assume first parameters in the parameter vector are the fixed ones.
num_params = length(xp);
num_gazes = length(gaze_angle);
if mod(num_params-CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG,num_gazes) ~= 0
    error('bad number of parameters in Curvefit_eye_centered');
else
    num_dyn_params = (num_params-CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG) / num_gazes;
end

P = zeros(1,CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG+num_dyn_params);  % the current parameter set
F = 0;   % total sum squared error
for f=1:num_gazes
    
    % this is confusing, but all it is doing is creating the parameters for
    % the current gaze angle by taking the fixed parameters in the
    % beginning of the whole parameter vector, and then taking the dynamic
    % parameters after that, depending on the current gaze angle being fit.
    P(1:CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG) = xp(1:CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG);
    P(CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG+1:num_dyn_params+CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG) = ...
        xp(num_dyn_params*(f-1)+CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG+1: ...
           num_dyn_params*f+CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG); 

    % modify the appropriate rotation angle, depending on the gaze angle
    % and direction.  assume that first parameter is azimuth rotation angle
    % and second parameter is elevation rotation angle.
    if fixation_type == CURVEFIT_VARY_FIXATION_X
        % this probably should work with the correct rotation
        %[P(1) P(2)] = Curvefit_rotate_coords(P(1),P(2),gaze_angle(f),0);
        % these should be equivalent for azimuth gaze shift?
        %P(1) = P(1) + gaze_angle(f);
        [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(azimuth,elevation,gaze_angle(f),0);
    elseif fixation_type == CURVEFIT_VARY_FIXATION_Y
        % not sure this will work, since amount of horizontal rotation
        % allowed for a vertical fixation change will be different
        % depending on the elevation.
        %[P(1) P(2)] = Curvefit_rotate_coords(P(1),P(2),0,gaze_angle(f));
        % clearly this does not allow for shift in azimuth.
        %P(2) = P(2) + gaze_angle(f);
        % this allows for shift in azimuth, depending on the elevation
        [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(azimuth,elevation,0,gaze_angle(f));

        %% test other gaze angles
        %hooch = [pi/4 0 -pi/4];
        %[azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(azimuth,elevation,0,hooch(f));
    else
        error('invalid fixation type in Curvefit_eye_centered');
    end
        
    % compute our fit data with the given parameters
    fit_data = feval( fun, P, azimuth_gaze, elevation_gaze );
    %fit_data = feval( fun, P, azimuth, elevation );
    curvefit_gaze_fits(f,:) = fit_data;  % kludge to save fit for each gaze
    
    % compute the error over all repititions for this gaze.
    % See DirectionTuningPlot_Curvefit for a description of the trial_data
    % structure.
    fit_data = repmat( fit_data, 1, num_complete_repititions );
    gaze_data(:,:) = trial_data(f,:,2:num_complete_repititions+1);  % reduce dims
    curvefit_gaze_sse(f) = ... % kludge to save sse for each gaze
        sum((Curvefit_sqrt_abs(fit_data) - Curvefit_sqrt_abs(gaze_data(:)')).^2);
    F = F + curvefit_gaze_sse(f);

    mean_data = trial_data(f,:,1); % kludge to save residual for each gaze
    curvefit_gaze_residuals(f,:) = curvefit_gaze_fits(f,:) - mean_data;
    curvefit_gaze_mean_sse(f) = sum((Curvefit_sqrt_abs(curvefit_gaze_fits(f,:)) - ...
        Curvefit_sqrt_abs(mean_data)).^2);  % kludge to save residuals versus mean
end

% make sure we return a real and finite value
if ~isreal(F) | ~all(isfinite(F))
    error(sprintf( 'Curvefit_eye_centered returning nonreal or nonfinite value' ));
end
