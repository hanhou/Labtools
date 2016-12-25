%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_head_centered.m -- This function computes the error over three
%--   gaze angles for a head centered model using the specified function
%--   as the model to fit each gaze angle individually.
%-- Created - pwatkins, 4/04
%-- Use only one non-linearity for both terms - Z. Briggs 6/09
%-----------------------------------------------------------------------------------------------------------------------
function F = Curvefit_head_centered(xp,fun,azimuth,elevation,trial_data, ...
    num_complete_repititions,gaze_angle,fixation_type)

Curvefit_defines;
% remap the parameters from the default so that non-linearities are
% constant over all gaze angles.
curvefit_parameter_mapping(CURVEFIT_PARAM_AZIMUTH) = 1;
curvefit_parameter_mapping(CURVEFIT_PARAM_ELEVATION) = 2;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN) = 3;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE) = 4;
curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET) = 5;
curvefit_parameter_mapping(CURVEFIT_PARAM_WEIGHT180) = 6;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN2) = 7;

% The head centered model keeps the rotation angles constant for each of
% the gaze angles.  
% Assume first parameters in the parameter vector are the fixed ones.
num_params = length(xp);
num_gazes = length(gaze_angle);
if mod(num_params-CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG,num_gazes) ~= 0
    error('bad number of parameters in Curvefit_head_centered');
else
    num_dyn_params = (num_params-CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG) / num_gazes;
end

P = zeros(1,CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG+num_dyn_params);  % the current parameter set
F = 0;   % total sum squared error
for f=1:num_gazes
    
    % this is confusing, but all it is doing is creating the parameters for
    % the current gaze angle by taking the fixed parameters in the
    % beginning of the whole parameter vector, and then taking the dynamic
    % parameters after that, depending on the current gaze angle being fit.
    P(1:CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG) = xp(1:CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG);
    P(CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG+1:num_dyn_params+CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG) = ...
        xp(num_dyn_params*(f-1)+CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG+1: ...
           num_dyn_params*f+CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG);
    
    % compute our fit data with the given parameters
    fit_data = feval( fun, P, azimuth, elevation );
    [primary,secondary]=Curvefit_cos_tuning_7p_asp_doublehalf_seperator(P,azimuth, elevation);
    %select=logical(fit_data>(max(trial_data(:,:,1))));
    %fit_data(select)=max(trial_data(:,:,1));
    %select=logical(fit_data<(min(trial_data(:,:,1))));
    %fit_data(select)=min(trial_data(:,:,1));
    curvefit_gaze_fits(f,:) = fit_data;  % kludge to save fit for each gaze
    curvefit_primary_fits(f,:) = primary;
    curvefit_secondary_fits(f,:) = secondary;
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
    error(sprintf( 'Curvefit_head_centered returning nonreal or nonfinite value' ));
end
