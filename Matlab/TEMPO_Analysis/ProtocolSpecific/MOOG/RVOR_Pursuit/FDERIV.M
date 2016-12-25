function [y] = fderiv(x,M,step)
% first derivative using Savitzky-Golay smoothing filters
% x can be in either deg or rads (must ve a VECTOR not MATRIX)
% M is the number of points to smooth
% compute SG-coefficients:
% c = [cM, cM-1, ..., c0, c-1, ..., c-M]'   in descending order
% x = [x1,...,xN], N > 2M                  vector of function values
% derivative: polynomial multiplication (convolution)
%             left shift by M points

% step = sampl_freq = 833.33;
  if(2*M < length(x))
    c = sg_filt(M,2,1);                % 1st derivative of 2nd order polynomial
    y=conv(c,x).* step;                % length(y) = length(x)+length(c)-1
                                       %           = length(x)+ 2M
    y=y(1+M:length(x)+M);              % left-shift by M points; adjust size
  else
%    y = x;
%    fprintf('0-order derivative\n')
    y = NaN * ones(size(x));
    fprintf('size too small, replace vector by NaNs\n')
  end
