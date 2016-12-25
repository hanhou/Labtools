%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_cos_tuning_7p.m -- This is the 6 parameter model of cos tuning as 
%--   a function of the position on the sphere (azimuth, elevation).
%-- Created - pwatkins, 4/04
%-- Modified for single nonlinearity - Z. Briggs 
%-----------------------------------------------------------------------------------------------------------------------
function F = Curvefit_cos_tuning_6p(xp,azimuth,elevation)

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
nonlin = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN));

% The rotation angles are used to compute the tuning function in a rotated 
% coordinate frame.  This is different than adding a phase angle to both 
% azimuth and elevation, because it allows for an axis between the peaks or 
% troughs of the fit that is not parallel to latitude or longitude planes 
% through the circle.
[rott rotp] = Curvefit_rotate_coords(t,p,t0,p0);

% compute the function value in the rotated coordinate frame.
g = cos(rott) .* cos(rotp);

% compute the same function 180 degrees out of phase
h = cos(rott - pi) .* cos(rotp);

% apply a variable non-linearity to the function to amplify peaks or troughs
if nonlin ~= 0
    g = ( exp(nonlin*g) - 1 ) / nonlin;
    h = ( exp(nonlin*h) - 1 ) / nonlin;
end


% use a weighted sum of the amplified harmonic terms.
% normalize the sum of these terms to [-1,1]
% so that the amplitude term still reflects the overall amplitude.
F = g + weight180*h;
% once in a blue moon the gods align to destroy my life and make F identically 0.
if ~all(F == 0)
    F = F ./ max(abs(F));
end

% add in an amplitude and dc offset.
F = amplitude*F + dc_offset;

% make sure we return a real and finite value
if ~isreal(F) | ~all(isfinite(F))
    error(sprintf( 'Curvefit_cos_tuning_6p returning nonreal or nonfinite value' ));
end
