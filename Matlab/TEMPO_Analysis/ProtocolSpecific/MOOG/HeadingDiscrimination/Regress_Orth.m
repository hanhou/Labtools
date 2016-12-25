% calculate the min err
function [err] = Regress_Orth(estimate_f)
global a_regress;

y = estimate_f(1,1)*a_regress(:,1) + estimate_f(1,2);
y_diff = abs( y - a_regress(:,2) );   % the vertical diffrence
yy_temp = sqrt( (y(end)-y(1))^2 + (a_regress(end,1)-a_regress(1,1))^2 );
err_r = y_diff*abs(a_regress(end,1)-a_regress(1,1)) / yy_temp;
%err = sum(err_r.^2);
err = sum(y_diff.^2);

