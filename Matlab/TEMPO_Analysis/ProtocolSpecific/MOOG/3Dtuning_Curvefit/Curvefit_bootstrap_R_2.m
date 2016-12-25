%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_bootstrap_R_2.m -- Callback function for bootstrapping of Vary_Fixation protocol data.
%     Computes R_2 correlation coefficient of re-sampled data, passed in by boostrap.
%--	pwatkins, 5/04
%-----------------------------------------------------------------------------------------------------------------------
function R_2 = Curvefit_bootstrap_R_2(bootstrap_data)

Curvefit_defines;

% use mean across repititions of the re-sampled data to compute R_2
mean_data = mean(bootstrap_data, 1);

% compute ss_yy using the mean data at each point, instead of trial data
mean_value = mean(Curvefit_sqrt_abs(mean_data));
SS_YY = sum(sum((Curvefit_sqrt_abs(mean_data) - mean_value).^2));

% compute SSE using mean data at each point, over all gaze angles
SSE = sum((Curvefit_sqrt_abs(curvefit_bootstrap_fits) - Curvefit_sqrt_abs(mean_data)).^2);

% compute the r^2 correlation coefficient.
R_2 = 1 - SSE / SS_YY;
