% VF_1D_Curvefit_err_eye.m -- computes the error between the square root of the
% data and the square root of the values from VF_1D_Curvefit

function err = VF_1D_Curvefit_err_eye(X)

global xdata ydata_merged stimtype

yfit_eye = VF_1D_Curvefit_eye(X,xdata);

% err = norm( sqrt(yfit_eye)-sqrt(ydata_merged{stimtype}).^2 );  does not converge, and does not perfectly follow R^2
err = sum(sum((sqrt(yfit_eye) - sqrt(ydata_merged{stimtype})).^2));

return;