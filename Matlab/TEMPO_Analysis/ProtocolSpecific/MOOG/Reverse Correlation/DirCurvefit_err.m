function err = DirCurvefit_err(x)

global xData yData
yfit = DirCurvefit(x,xData);

%threshold the fitted values (don't allow less than zero)
% yfit(yfit < 0) = 0;

err = norm(sqrt(yfit)-sqrt(yData))^2;

return;