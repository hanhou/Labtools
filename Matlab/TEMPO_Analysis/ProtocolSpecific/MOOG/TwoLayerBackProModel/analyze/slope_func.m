% calculate the min err
function err = slope_func(estimate)
global maxx;
x_real = [-40,-20,0,20,40];
y = estimate(1,1)*x_real + estimate(1,2);
yy_temp = sqrt( (y(4)-y(2))^2 + 40^2 );
err_r = abs(maxx(:) - y(:))*40 / yy_temp;
err = sum(err_r.^2);