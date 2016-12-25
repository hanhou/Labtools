%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_cos_halfrect_tuning_7p.m -- This is the 7 parameter model of cos tuning as 
%--   a function of the position on the sphere (azimuth, elevation).  Model consists of 
%--   two components, 180 degrees opposed with same nonlinearities applied and both half 
%--   wave rectified at 0, before adding dc_offset.
%-- Created - pwatkins, 4/04
%-- Modified for single nonlinearity - Z. Briggs 
%-- Modified to clean up half wave rectify - pwatkins, 8/04
%--    xxx - need to verify this is what we want.
%--      previous code from Zack only rectified azimuth component of 180 degree term.
%--      rectification of F was performed after dc offset is added in, which does nothing for dc_offset > 1.
%-----------------------------------------------------------------------------------------------------------------------
function F = Curvefit_cos_tuning_7p_halfrectmodel(xp,azimuth,elevation)

Curvefit_defines;

% Just to be cryptic, and make the formulas below more concise, 
% use these variables:
t = azimuth;        % theta, angle in azimuth
p = elevation;      % phi, angle in elevation
t0 = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_AZIMUTH));  % theta naught, rotation angle in azimuth
p0 = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_ELEVATION)); % phi naught, rotation angle in elevation
num_points = length(azimuth);  % should be same as length(elevation)

amplitude = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE));
dc_offset = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET));
weight180 = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_WEIGHT180));
nonlin_az = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN_AZ));
nonlin_el = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN_EL));

% The rotation angles are used to compute the tuning function in a rotated 
% coordinate frame.  This is different than adding a phase angle to both 
% azimuth and elevation, because it allows for an axis between the peaks or 
% troughs of the fit that is not parallel to latitude or longitude planes 
% through the circle.
[rott rotp] = Curvefit_rotate_coords(t,p,t0,p0);

% compute the function value in the rotated coordinate frame.
%g = cos(rott) .* cos(rotp);

% compute the same function 180 degrees out of phase
%h = cos(rott - pi) .* cos(rotp);

% use max of function to normalize after nonlinearity is applied.
max_nonlin_az = 1;
max_nonlin_el = 1;

% apply a variable non-linearity to the function to amplify peaks or troughs.
% apply seperately to azimuth and elevation terms (aspect ratio).

e = 1e-8;  % make sure max_nonlin is not zero

if abs(nonlin_az) > e
    ga = ( exp(nonlin_az*cos(rott)) - 1 ) / (nonlin_az);
    ha = ( exp(nonlin_az*cos(rott-pi)) - 1 ) / (nonlin_az);
    max_nonlin_az = ( exp(nonlin_az*max_nonlin_az) - 1 ) / nonlin_az;
else
    ga = cos(rott);
    ha = cos(rott-pi); %sel = (ha < 0); ha(sel) = 0; % zack's rectifier
end

if abs(nonlin_el) > e
    gb = ( exp(nonlin_el*cos(rotp)) - 1 ) / (nonlin_el);
    hb = ( exp(nonlin_el*cos(rotp)) - 1 ) / (nonlin_el);
    max_nonlin_el = ( exp(nonlin_el*max_nonlin_el) - 1 ) / nonlin_el;
else
    gb = cos(rotp);
    hb = cos(rotp);
end

% multiply to get cosine tuning in azimuth and elevation.
% half wave rectify the main term and 180 degree term separately.
g = ga.*gb; sel = (g < 0); g(sel) = 0;
h = ha.*hb; sel = (h < 0); h(sel) = 0;
max_nonlin = max_nonlin_az*max_nonlin_el;

% kludge for visualizing two components of this model.
curvefit_primary_fits = g;
curvefit_secondary_fits = h; 

% use a weighted sum of the cosine terms.
% normalize the sum of these terms to [-1,1]
% so that the amplitude term still reflects the overall amplitude.
F = (g + weight180*h) / (max_nonlin + abs(weight180)*max_nonlin);

% add in an amplitude and dc offset.
F = amplitude*F + dc_offset;

% NEED floor here GCD

% make sure we return a real and finite value
if ~isreal(F) | ~all(isfinite(F))
    error(sprintf( 'Curvefit_cos_halfrect_tuning_7p returning nonreal or nonfinite value' ));
end
