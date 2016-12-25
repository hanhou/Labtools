% VF_1D_Curvefit_err_head.m -- computes the error between the square root of the
% data and the square root of the values from VF_1D_Curvefit_head

function err = VF_1D_Curvefit_err_head(X)

global xdata ydata_merged stimtype

yfit_head = VF_1D_Curvefit_head(X,xdata);

% err = norm( sqrt(yfit_head)-sqrt(ydata_merged{stimtype}).^2 );   does not converge, and does not perfectly follow R^2
err = sum(sum((sqrt(yfit_head) - sqrt(ydata_merged{stimtype})).^2));

return;