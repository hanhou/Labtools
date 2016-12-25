%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_gaus1_tuning_2p.m -- This is the 3 parameter model of gaussian
%--   tuning in one dimension.
%-- Created - pwatkins, 4/04
%-----------------------------------------------------------------------------------------------------------------------
function F = Curvefit_gaus1_tuning_3p(xp,position)

Curvefit_defines;

mean_p = xp(1);
var_p = xp(2);
amplitude = xp(3);

F = amplitude * exp( -(position - mean_p).^2 / (2*var_p.^2) );

% make sure we return a real and finite value
if ~isreal(F) | ~all(isfinite(F))
    error(sprintf( 'Curvefit_gaus1_tuning_2p returning nonreal or nonfinite value' ));
end
