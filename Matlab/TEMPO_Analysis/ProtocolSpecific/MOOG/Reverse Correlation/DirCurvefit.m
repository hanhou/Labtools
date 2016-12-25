function F = DirCurvefit(x,xData)

% % Function #3: wrapped Gaussian, or two (Yang and Maunsell 2004)
% 
% % re-write exponential numerator according to wrapped Gaussian conditions (Yang and Maunsell 2004):
% % if <= pi, exp_num =  ( [0 + xData] - x(2)) ^2;
% % if >= pi, exp_num = ([-2*pi + xData] - x(2)) ^2;
% 
% % xdata_mod = xData;
% % xdata_mod(abs(xData - x(2)) > pi) = -2*pi + xData(abs(xData - x(2)) > pi);
% % exp_num = (xdata_mod - x(2));
% 
% %exp_num = [];
% xdata_size = size(xData);
% exp_num = zeros(xdata_size);
% %for t = 1: length(xData(1,:))
% for t = 1: prod(xdata_size)
%     if abs(xData(t)-x(2)) <= pi
%         exp_num(t) = xData(t) - x(2);
%     elseif abs(xData(t)-x(2)) > pi
%         exp_num(t) = 2*pi - xData(t) + x(2);
%     else
%         error('something is wrong with exp_num');
%     end
% end
% 
% % 4 params: x = [A mu sigma DC]
% F = x(1) * exp(-(exp_num.^2)/(x(3)^2)) + x(4);
% F(F < 0) = 0;

F = x(1) * exp(-2*(1-cos(xData-x(2)))/x(3)^2)  + x(4);
%**********************************************************************************

% threshold the fitted values (don't allow less than zero)
F(F < 0) = 0;

return;