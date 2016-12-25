%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_eye_centered_bounds.m -- Create bounds for eye centered model.
%-- Created - pwatkins, 4/04
%-----------------------------------------------------------------------------------------------------------------------
function R = Curvefit_eye_centered_bounds(trial_data, num_gazes)

Curvefit_defines;

dyn_terms = 0;  % no define - dependent on this module anyways,
                % and computed elsewhere as: 
                % (total_terms - FIXED_PARAMS) / num_gazes
terms = CURVEFIT_NUM_FIXED_PARAMETERS_EC_5GF + dyn_terms*num_gazes;

% create bounds for the fixed parameters for this model.
% the fixed parameters are azimuth and elevation rotation angles.
lb = [-2*pi  -2*pi   0  -50  0    -50  0   ];
ub = [ 2*pi   2*pi  50   50  5e3   50  5e3 ];

% % create bounds for the dynamic parameters for this model.
% for j=1:num_gazes
%     i = (j-1)*dyn_terms + CURVEFIT_NUM_FIXED_PARAMETERS_EC_5GF+1;
%     lb(i:i+dyn_terms-1) = [0     0    ];
%     ub(i:i+dyn_terms-1) = [5e3   5e3  ];
% end

% Compute a guess at the amplitude by taking the amplitude of the average
% data over all repititions.
% See DirectionTuningPlot_Curvefit for a description of the trial_data
% structure.
if length(trial_data) == 0
    amplitude = 100*rand(1,num_gazes);
else
    for i=1:num_gazes
        amplitude(i) = (max(trial_data(i,:,1)) - min(trial_data(i,:,1)))/2;
    end
end

% create a completely random but tightly bounded initial x0 guess.

% create initial guess for the fixed parameters for this model.
% the fixed parameters are azimuth and elevation rotation angles.
x0 = rand(1,CURVEFIT_NUM_FIXED_PARAMETERS_EC_5GF).*[2*pi pi 2 100 0 100 0] ...
    - [pi pi/2 0 50 -amplitude(2) 50 -amplitude(2)];

% % create initial guess for the dynamic parameters for this model.
% for j=1:num_gazes
%     i = (j-1)*dyn_terms + CURVEFIT_NUM_FIXED_PARAMETERS_EC_5GF+1;
%     x0(i:i+dyn_terms-1) = rand(1,dyn_terms).*[amplitude(j) amplitude(j)] - [0 0];
% end

R{CURVEFIT_BOUNDS_X0} = x0;
R{CURVEFIT_BOUNDS_LB} = lb;
R{CURVEFIT_BOUNDS_UB} = ub;
