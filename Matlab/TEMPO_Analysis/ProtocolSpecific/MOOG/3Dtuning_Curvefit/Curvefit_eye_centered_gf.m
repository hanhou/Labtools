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
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN180) = 4;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE) = 5;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_SLOPE) = 5;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_OFFSET) = 6;
curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET) = 7;
curvefit_parameter_mapping(CURVEFIT_PARAM_WEIGHT180) = 8;

% The eye centered model varies the rotation angles for each of
% the gaze angles, corresponding to the actual gaze angle and direction.
% Assume first parameters in the parameter vector are the fixed ones.
num_params = length(xp);
num_gazes = length(gaze_angle);
if mod(num_params-CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF,num_gazes) ~= 0
    error('bad number of parameters in Curvefit_eye_centered');
else
    num_dyn_params = (num_params-CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF) / num_gazes;
end

P = zeros(1,CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF+num_dyn_params);  % the current parameter set
F = 0;   % total sum squared error
for f=1:num_gazes
    
    % this is confusing, but all it is doing is creating the parameters for
    % the current gaze angle by taking the fixed parameters in the
    % beginning of the whole parameter vector, and then taking the dynamic
    % parameters after that, depending on the current gaze angle being fit.
    P(1:CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF) = xp(1:CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF);
    P(CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF+1:num_dyn_params+CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF) = ...
        xp(num_dyn_params*(f-1)+CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF+1: ...
           num_dyn_params*f+CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF); 

    % modify the appropriate rotation angle, depending on the gaze angle
    % and direction.
    if fixation_type == CURVEFIT_VARY_FIXATION_X
        [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(azimuth,elevation,gaze_angle(f),0);
    elseif fixation_type == CURVEFIT_VARY_FIXATION_Y
        [azimuth_gaze elevation_gaze] = Curvefit_rotate_coords(azimuth,elevation,0,gaze_angle(f));
    else
        error('invalid fixation type in Curvefit_eye_centered');
    end

    % modify the amplitude, depending on the gaze angle to emulate gain fields between gaze angles.    
    P(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE)) = ...
        P(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_SLOPE))*gaze_angle(f) + ...
        P(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_OFFSET));
    
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
