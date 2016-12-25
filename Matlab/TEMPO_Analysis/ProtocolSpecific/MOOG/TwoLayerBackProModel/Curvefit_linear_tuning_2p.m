%-----------------------------------------------------------------------------------------------------------------------
%-- Curvefit_gaus1_tuning_2p.m -- This is the 2 parameter model of gaussian
%--   tuning in one dimension.
%-- Created - pwatkins, 4/04
%-----------------------------------------------------------------------------------------------------------------------
function F = Curvefit_linear_tuning_2p(xp,position)

Curvefit_defines;

slope = xp(1);
intercept = xp(2);

F = slope*position + intercept;

% make sure we return a real and finite value
if ~isreal(F) | ~all(isfinite(F))
    error(sprintf( 'Curvefit_gaus1_tuning_2p returning nonreal or nonfinite value' ));
end
