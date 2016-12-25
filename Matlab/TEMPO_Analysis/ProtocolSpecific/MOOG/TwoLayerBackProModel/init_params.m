%-----------------------------------------------------------------------------------------------------------------------
%-- Initialize parameters for the input and output neurons, based
%--   on some tuning models.
%-- Created - pwatkins, 6/04
%-----------------------------------------------------------------------------------------------------------------------
function [Pvis, Pves, num_neurons, P_eye, num_eye_neurons] = ...
    init_params(fun,range_az,range_el,range_gaze_az,range_gaze_el)

% % create a distribution of eye position - Curvefit_gaus1_tuning_3p
% num_eye_neurons = 5;
% fixed_params = repmat([15/180*pi 15], num_eye_neurons, 1);
% P_eye = [range_gaze_az(1):(range_gaze_az(2)-range_gaze_az(1))...
%         /(num_eye_neurons-1):range_gaze_az(2)];
% %gaze_el = rand(1,1).*(range_gaze_el(2)-range_gaze_el(1)) + range_gaze_el(1);
% P_eye = [P_eye' fixed_params];

% create a distribution of eye position - Curvefit_linear_tuning_2p
num_eye_neurons = 6;
%P_eye = [40/pi 15; -40/pi 15; 20/pi 15; -20/pi 15; 10/pi 15];
%P_eye = [90/pi 15; -90/pi 15; 45/pi 15; -45/pi 15; 27/pi 15];
%P_eye = [315/4/pi 17.5; -315/4/pi 17.5; 315/8/pi 17.5; -315/8/pi 17.5; 315/16/pi 17.5];
%P_eye = [9/2/pi 0; -9/2/pi 0; 9/4/pi 0; -9/4/pi 0; 9/8/pi 0];

%P_eye = [9/2/pi 0; -9/2/pi 0; 27/4/pi 0; -27/4/pi 0; 9/4/pi 0; -9/4/pi 0; 9/8/pi 0; -9/8/pi 0];
P_eye = [9/2/pi 0; -9/2/pi 0; 27/8/pi 0; -27/8/pi 0; 9/4/pi 0; -9/4/pi 0];
%P_eye = [9/2/pi 0; -9/2/pi 0; 27/4/pi 0; -27/4/pi 0];
%P_eye = [9/2/pi 0; -9/2/pi 0];

% create parameters using Curvefit_cos_tuning_5p model
Curvefit_defines;
curvefit_parameter_mapping(CURVEFIT_PARAM_AZIMUTH) = 1;
curvefit_parameter_mapping(CURVEFIT_PARAM_ELEVATION) = 2;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN) = 3;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE) = 4;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_SLOPE) = 4;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE_OFFSET) = 5;
curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET) = 6;
curvefit_parameter_mapping(CURVEFIT_PARAM_DC_SLOPE) = 6;
curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSETB) = 7;
create_sample_points;
num_neurons = num_unique_points;
% fix nonlin, amplitude slope and offset, dc slope and offset
%fixed_params = repmat([1 0 25 0 10], num_neurons, 1);
%fixed_params = repmat([2 0 1 0 0], num_neurons, 1);
fixed_params = repmat([1 0 1.46211715726001 0 -0.46211715726001], num_neurons, 1); % [-1 1]
%fixed_params = repmat([1 0 1.27935251260251 0 -0.279352512602509], num_neurons, 1); % [-0.75 1]
%fixed_params = repmat([2 0 1.76159415595577 0 -0.761594155955765], num_neurons, 1); % [-1 1]
Pvis = [unique_point_azimuth_r' unique_point_elevation_r' fixed_params];

fixed_params = repmat([1 0 1.46211715726001 0 -0.46211715726001], num_neurons, 1); % [-1 1]
%fixed_params = repmat([1 0 1.27935251260251 0 -0.529352512602509], num_neurons, 1); % [-1 0.75]
%fixed_params = repmat([1 0 1.09658786794501 0 -0.596587867945007], num_neurons, 1); % [-1 0.5]
%fixed_params = repmat([2 0 1.76159415595577 0 -0.761594155955765], num_neurons, 1); % [-1 1]
Pves = [unique_point_azimuth_r' unique_point_elevation_r' fixed_params];
