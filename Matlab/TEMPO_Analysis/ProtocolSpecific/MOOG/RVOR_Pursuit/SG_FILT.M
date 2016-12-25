function [c] = sg_filt(N,M,ld)

% Savitzky-Golay filter (see Press et al, 1992)
%
% poly =a0 + a1*i1 + a2*i2 + a3*i3 + ... aM*iM     polynomial in i
% to be fitted to the values f-nL, ..., fnR
%
% design matrix Aij = i**j  where  i= -nL, ... nR  (evaluated at data points)
%                                  j = 0, ..., M   (basis functions)
%
% Normal equations: (AT * A) * a = AT * f
%
% cn = Sum(0 to M) {inv(AT * A)}0m n**m
%
% M                  order of polynomial in i
% ld                 ld =0: smoothing
%
% march 1996

 nR = N;
 nL = nR;                   % symetrical moving window: 2*nR+1 points

 indx = (-nL:nR)';
 A=[indx.^0];
 for j = 1:M
     A = [A indx.^j];       % design matrix
 end

% smoothing
 c = zeros(1,nL+nR+1);
 B = inv(A' * A);
 for n = 0:nL
     for m = 0:M
         c(nL+1-n) = c(nL+1-n) + B(ld+1,m+1)*(n)^m;
     end
 end

 for n = 1:nR
     for m = 0:M
         c(nL+1+n) = c(nL+1+n) + B(ld+1,m+1)*(-n)^m;
     end
 end
