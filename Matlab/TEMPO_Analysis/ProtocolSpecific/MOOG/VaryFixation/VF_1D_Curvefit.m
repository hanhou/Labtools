% VF_1D_Curvefit.m -- function for fitting 1D heading tuning data in vary_fixation experiment -- CRF

function F = VF_1D_Curvefit(x,xdata)

%**********************************************************************************
%
% % Function #1: modified sinusoid (Nguyenkim and DeAngelis 2003), 
% % with second term 180 deg out of phase and half-wave rectified:
% 
% % F = A1 * [ G(sin(f*theta + phi)) + {K * G(sin(f*theta + phi + 180))}_+ ] + DC
% % where G(x) = (exp(n*x)-1) / n, f is constrained to 1, and K is bounded between 0 and 1.
% 
% % Thus 5 params: x = [A n phi K DC]
% 
% F = x(1) * ( (exp(x(2)*(sin(xdata + x(3))))-1) / x(2) + ...
%     x(4) * halfrect((exp(x(2)*(sin(xdata + x(3) + pi)))-1) / x(2)) ) + x(5); 
% 
% 
%**********************************************************************************
% 
% % Function #2: two gaussians with independent widths, where the second 
% % has a mean 180 deg from the first, and a K value that can range from -1 to 1.
% 
% % Thus 6 params: x = [A mu sigma1 K sigma2 DC]
% 
% % (to add/subtract 180 appropriately for 2nd peak)
% if x(2) > pi
%     temp_pi = -pi;
% else
%     temp_pi = pi;
% end
% 
% F = x(1) * ( exp((-(xdata - x(2)).^2) / (2 * x(3)^2)) + ... 
%     x(4) * exp((-(xdata - (x(2) + temp_pi)).^2) / (2 * x(5)^2)) ) + x(6);
% 
% 
%**********************************************************************************
% Function #3: wrapped Gaussian, or two (Yang and Maunsell 2004)

% % re-write exponential numerator according to wrapped Gaussian conditions (Yang and Maunsell 2004):
% % if <= pi, exp_num =  ( [0 + xdata] - x(2)) ^2;
% % if >= pi, exp_num = ([-2*pi + xdata] - x(2)) ^2;
% 
% % xdata_mod = xdata;
% % xdata_mod(abs(xdata - x(2)) > pi) = -2*pi + xdata(abs(xdata - x(2)) > pi);
% % exp_num = (xdata_mod - x(2));
% 
% %exp_num = [];
% xdata_size = size(xdata);
% exp_num = zeros(xdata_size);
% %for t = 1: length(xdata(1,:))
% for t = 1: prod(xdata_size)
%     if abs(xdata(t)-x(2)) <= pi
%         exp_num(t) = xdata(t) - x(2);
%     elseif abs(xdata(t)-x(2)) > pi
%         exp_num(t) = 2*pi - xdata(t) + x(2);
%     else
%         error('something is wrong with exp_num');
%     end
% end
% 
% 
% % 4 params: x = [A mu sigma DC]
% % F = x(1) * exp(-(exp_num.^2)/(x(3)^2)) + x(4);
% 
% % 6 params: x = [A mu sigma1 K sigma2 DC]
% if x(2) > pi   % (to add/subtract 180 appropriately, for 2nd peak)
%     temp_pi = -pi;
% else
%     temp_pi = pi;
% end
% F = x(1) * ( exp(-(exp_num.^2)/(x(3)^2)) + x(4)*exp(-((exp_num+temp_pi).^2)/(x(5)^2))) + x(6);

%**********************************************************************************
% Function #4: Charlie Special (modified wrapped gaussian)

% 6 params: x = [A mu sigma K K-sig DC]
F = x(1) * [ exp(-2*(1-cos(xdata-x(2)))/(x(5)*x(3))^2) + x(4)*exp(-2*(1-cos(xdata-x(2)-pi))/x(3)^2) ] + x(6);


%**********************************************************************************

% threshold the fitted values (don't allow less than zero)
F(F < 0) = 0;

return;