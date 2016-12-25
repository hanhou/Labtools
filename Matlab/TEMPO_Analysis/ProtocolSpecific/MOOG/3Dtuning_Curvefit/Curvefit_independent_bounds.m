%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_independent_bounds.m -- Create bounds for independent parameters model.
%-- Created - pwatkins, 4/04
%-----------------------------------------------------------------------------------------------------------------------
function R = Curvefit_independent_bounds(trial_data, num_gazes)

Curvefit_defines;

dyn_terms = 7;  % no define - dependent on this module anyways,
                % and computed elsewhere as: 
                % (total_terms - FIXED_PARAMS) / num_gazes
terms = CURVEFIT_NUM_FIXED_PARAMETERS_IND + dyn_terms*num_gazes;

% create bounds for the fixed parameters for this model.
% the fixed parameters are azimuth and elevation rotation angles.
%lb = [];
%ub = [];

% create bounds for the dynamic parameters for this model.
for j=1:num_gazes
    i = (j-1)*dyn_terms + CURVEFIT_NUM_FIXED_PARAMETERS_IND+1;
    lb(i:i+dyn_terms-1) = [-2*pi  -2*pi 0     0     0 -50 -50];
    ub(i:i+dyn_terms-1) = [ 2*pi   2*pi 5e3   5e3   1  50  50];
end

% Compute a guess at the amplitude by taking the amplitude of the average
% data over all repititions.
% See DirectionTuningPlot_Curvefit for a description of the trial_data
% structure.
for i=1:num_gazes
    amplitude(i) = (max(trial_data(i,:,1)) - min(trial_data(i,:,1)))/2;
end

% create a completely random but tightly bounded initial x0 guess.

% create initial guess for the fixed parameters for this model.
% the fixed parameters are azimuth and elevation rotation angles.
%x0 = rand(1,CURVEFIT_NUM_FIXED_PARAMETERS_IND).*[] - [];

% create initial guess for the dynamic parameters for this model.
for j=1:num_gazes
    i = (j-1)*dyn_terms + CURVEFIT_NUM_FIXED_PARAMETERS_IND+1;
    x0(i:i+dyn_terms-1) = rand(1,dyn_terms).*[2*pi pi amplitude(j) ...
            amplitude(j) 1 2 2] - [pi pi/2 0 0 0 1 1];
end

R{CURVEFIT_BOUNDS_X0} = x0;
R{CURVEFIT_BOUNDS_LB} = lb;
R{CURVEFIT_BOUNDS_UB} = ub;
