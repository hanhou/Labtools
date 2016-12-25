%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_defines.m -- Defines specific for curve fitting (see DirectionTuningPlot_Curvefit.m)
%-- Created - pwatkins, 4/04
%-----------------------------------------------------------------------------------------------------------------------

CURVEFIT_MAX_REPITITIONS = 50;
CURVEFIT_NUM_FIXED_PARAMETERS_HC = 2;
CURVEFIT_NUM_FIXED_PARAMETERS_EC = 2;
CURVEFIT_NUM_FIXED_PARAMETERS_HC_NG = 4;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG = 4;
CURVEFIT_NUM_FIXED_PARAMETERS_HC_5NG = 3;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_5NG = 3;
CURVEFIT_NUM_FIXED_PARAMETERS_HC_6P_NG = 6;
CURVEFIT_NUM_FIXED_PARAMETERS_HC_7P_ASP_NG = 7;
CURVEFIT_NUM_FIXED_PARAMETERS_HC_8P_ASP_NG = 8;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_NG_VF = 5;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_GF = 6;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_GFDC = 8;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_GFB = 10;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_GFADC = 7;
CURVEFIT_NUM_FIXED_PARAMETERS_EC_5GF = 7;
CURVEFIT_NUM_FIXED_PARAMETERS_IND = 0;
CURVEFIT_BOUNDS_X0 = 1;
CURVEFIT_BOUNDS_LB = 2;
CURVEFIT_BOUNDS_UB = 3;
CURVEFIT_MAX_PARAMETERS = 50;
CURVEFIT_ROTATION_FIT = 1;
CURVEFIT_ROTATION_GAZE = 2;
CURVEFIT_MAX_OPT_FILE_SIZE = 65535;
CURVEFIT_INVALID_GAZE_ANGLE = 361;

% this is the default parameter mapping
CURVEFIT_PARAM_AZIMUTH = 1;
CURVEFIT_PARAM_ELEVATION = 2;
CURVEFIT_PARAM_AMPLITUDE = 3;
CURVEFIT_PARAM_DC_OFFSET = 4;
CURVEFIT_PARAM_WEIGHT180 = 5;
CURVEFIT_PARAM_NONLIN = 6;
CURVEFIT_PARAM_NONLIN180 = 7;
CURVEFIT_PARAM_WEIGHTFIX = 8;
CURVEFIT_PARAM_AMPLITUDE_SLOPE = 9;
CURVEFIT_PARAM_AMPLITUDE_OFFSET = 10;
CURVEFIT_PARAM_DC_SLOPE = 11;
CURVEFIT_PARAM_DC_OFFSETB = 12;
CURVEFIT_PARAM_WEIGHT180_SLOPE = 13;
CURVEFIT_PARAM_WEIGHT180_OFFSET = 14;
CURVEFIT_PARAM_NONLIN2 = 15;
CURVEFIT_PARAM_NONLIN_AZ = 16;
CURVEFIT_PARAM_NONLIN_EL = 17;
global curvefit_parameter_mapping;

% these are global variables, all kludges to pass data.
global curvefit_gaze_fits;
global curvefit_primary_fits;
global curvefit_secondary_fits;
global curvefit_gaze_residuals;
global curvefit_gaze_sse;
global curvefit_gaze_mean_sse;
global curvefit_R2_export;
global curvefit_Chi2_export;
global curvefit_fitting_stims_export;
global curvefit_fitting_param_export;
global curvefit_repititions_export;
global curvefit_gaze_angle_export;
global curvefit_fixation_type_export;
global curvefit_R2_distrib_sig_export;
global gendata_noise_type;
global curvefit_bootstrap_fits;

% enumeration of all possibilites of r^2 distrib comparisons
CURVEFIT_R2_DISTRIB_EQUAL_MEANS = 0;
CURVEFIT_R2_DISTRIB_1_NSIG = -1;
CURVEFIT_R2_DISTRIB_2_NSIG = -2;
CURVEFIT_R2_DISTRIB_1_SIG = 1;
CURVEFIT_R2_DISTRIB_2_SIG = 2;

% xxx - check if these are defined somewhere already
CURVEFIT_VARY_FIXATION_X = 0;
CURVEFIT_VARY_FIXATION_Y = 1;
CURVEFIT_OCCULAR_STIM = 2;
CURVEFIT_VESTIBULAR_STIM = 1;
CURVEFIT_OCCULAR_AND_VESTIBULAR_STIM = 3;
global curvefit_stim_string;
curvefit_stim_string = {'vestibular only' 'occular only' 'vestibular + occular'};
