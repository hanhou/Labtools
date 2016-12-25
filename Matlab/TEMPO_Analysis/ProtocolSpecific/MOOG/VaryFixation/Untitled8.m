xdata = 0:0.01:2*pi;

% 6 params: x = [A mu sigma K K-sig DC]
x = [10 pi pi/3 0 1 0];

F = x(1) * ( exp(-2*(1-cos(xdata-x(2)))/(x(5)*x(3))^2) + x(4)*exp(-2*(1-cos(xdata-x(2)-pi))/x(3)^2) ) + x(6);

figure; plot(xdata,F); title(num2str(x))