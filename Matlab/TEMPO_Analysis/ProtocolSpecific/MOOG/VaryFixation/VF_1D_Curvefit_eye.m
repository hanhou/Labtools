% VF_1D_Curvefit_eye.m -- function for fitting 1D heading tuning data in vary_fixation experiment:
% fitting all three gaze angles simultaneously, while constraining the mu parameter (peak) to be eye-centered

function F = VF_1D_Curvefit_eye(X,xdata)

% Function #4: Charlie Special (variant wrapped gaussian that can also be sinusoidal)
% 6 params: x = [A mu sigma K K-sig DC]
% but now 16 params: all are independent (i.e., x3) except mu

F_minus = X(1,1) * ( exp(-2*(1-cos(xdata-(X(2,2)+pi/8)))/(X(1,5)*X(1,3))^2) + X(1,4)*exp(-2*(1-cos(xdata-(X(2,2)+pi/8)-pi))/X(1,3)^2) ) + X(1,6);
F_zero = X(2,1) * ( exp(-2*(1-cos(xdata-X(2,2)))/(X(2,5)*X(2,3))^2) + X(2,4)*exp(-2*(1-cos(xdata-X(2,2)-pi))/X(2,3)^2) ) + X(2,6);
F_plus = X(3,1) * ( exp(-2*(1-cos(xdata-(X(2,2)-pi/8)))/(X(3,5)*X(3,3))^2) + X(3,4)*exp(-2*(1-cos(xdata-(X(2,2)-pi/8)-pi))/X(3,3)^2) ) + X(3,6);

F = [F_minus; F_zero; F_plus];

% threshold the fitted values (don't allow less than zero)
F(F < 0) = 0;

return;