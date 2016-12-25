function v=vonMises(r,p)
% function v=vonMises(r,p)
% von Mises distribution
% v(r,p)=exp(p(2)*cos(r-p(1)))/(2*pi*I0(p(2)));
% I0(p(2)):modified bessel function of order 0
% r: angle in rad
% p(1): mean angle in rad
% P(2): parameter of concentration

v=exp(p(2)*cos(r-p(1)))/(2*pi*besseli(0,p(2)));
