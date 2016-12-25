%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_cos_tuning_7p.m -- This is the 6 parameter model of cos tuning as 
%--   a function of the position on the sphere (azimuth, elevation).
%-- Created - pwatkins, 4/04
%-- Modified for single nonlinearity - Z. Briggs 
% NOTE: this version needs to take vector inputs for azimuth and elevation,
% GCD 9/8/04
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
nonlin = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN));
nonlin2 = xp(curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN2));
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

% apply a variable non-linearity and aspect ratio to the function to amplify peaks or troughs
if (nonlin ~= 0)
    ga=( exp(nonlin*cos(rott)) - 1 ) / (nonlin);
    ha=( exp(nonlin*cos(rott-pi)) - 1 ) / (nonlin);
    
else
    ga = cos(rott);
    ha = cos(rott-pi);
end

if (nonlin2 ~= 0)
    gb = ( exp(nonlin2*cos(rotp)) - 1 ) / (nonlin2);
    hb = ( exp(nonlin2*cos(rotp)) - 1 ) / (nonlin2);
    
else
    gb = cos(rotp);
    hb = cos(rotp);
end
g=ga.*gb;
selecta=logical(ha<0);
selectb=logical(hb<0);
ha(selecta)=0;
%hb(selectb)=0;
h=ha.*hb;

% use a weighted sum of the amplified harmonic terms.
% normalize the sum of these terms to [-1,1]
% so that the amplitude term still reflects the overall amplitude.

F = g + weight180*h;

% commented out for use in d' function - commented lines take place in
% var_mean_fit.m

% once in a blue moon the gods align to destroy my life and make F identically 0.
 if ~all(F == 0)
     F = F ./ max(abs(F));
 end
 
% % add in an amplitude and dc offset.
F = amplitude*F + dc_offset;
select=logical(F<0);
F(select)=0;

% make sure we return a real and finite value
if ~isreal(F) | ~all(isfinite(F))
    error(sprintf( 'Curvefit_cos_tuning_6p returning nonreal or nonfinite value' ));
end
