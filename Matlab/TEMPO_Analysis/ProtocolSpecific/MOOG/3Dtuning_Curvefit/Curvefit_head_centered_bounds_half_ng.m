%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_head_centered_bounds.m -- Create bounds for head centered model.
%-- Created - pwatkins, 4/04
%-----------------------------------------------------------------------------------------------------------------------
function R = Curvefit_head_centered_bounds(trial_data, num_gazes)

Curvefit_defines;

dyn_terms = 0;  % no define - dependent on this module anyways,
                % and computed elsewhere as: 
                % (total_terms - FIXED_PARAMS) / num_gazes
terms = CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG + dyn_terms*num_gazes;

% create bounds for the fixed parameters for this model.
% the fixed parameters are azimuth and elevation rotation angles.
lb = [-pi/2  -pi/2   0   0   0   -1 0];
ub = [ 3*pi/2   pi/2  50   5e3 5e3 1 50];

% Compute a guess at the amplitude by taking the amplitude of the average
% data over all repititions.
% See DirectionTuningPlot_Curvefit for a description of the trial_data
% structure.
if length(trial_data) == 0
    amplitude = 100*rand(1,num_gazes);
else
    amplitude = (max(trial_data(:,:,1)) - min(trial_data(:,:,1)))/2;
end

% create a completely random but tightly bounded initial x0 guess.

% create initial guess for the fixed parameters for this model.
% the fixed parameters are azimuth and elevation rotation angles.
x0 = rand(1,CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG).*[2*pi pi 2 0 0 2 2] - [pi pi/2 0 -amplitude -amplitude 1 0];



R{CURVEFIT_BOUNDS_X0} = x0;
R{CURVEFIT_BOUNDS_LB} = lb;
R{CURVEFIT_BOUNDS_UB} = ub;
