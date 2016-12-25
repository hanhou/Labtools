function ci = nlparci(beta,resid,J)
%NLPARCI Confidence intervals on parameters of nonlinear models.
%   CI = NLPARCI(BETA,RESID,J) returns the 95% confidence interval CI
%   on the nonlinear least squares parameter estimate BETA, given the 
%   residuals RESID, and the Jacobian matrix J at the solution.
%
%   The confidence interval calculation is valid for systems where 
%   the number of rows of J exceeds the length of BETA. 
%
%   NLPARCI uses the outputs of NLINFIT for its inputs.
%   Example:
%      [beta,resid,J]=nlinfit(input,output,'f',betainit);
%      ci = nlparci(beta,resid,J);
%
%   See also NLINFIT.
%

%   Bradley Jones 1-28-94
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.12 $  $Date: 2002/01/17 21:31:28 $

%initialization
if nargin < 3
   error('Requires three inputs.');
end;

resid = resid(:);
[m,n] = size(J);
if m <= n
   error('The number of observations must exceed the number of parameters.');
end;

if length(beta) ~= n
   error('The length of beta must equal the number of columns in J.')
end

% approximation when a column is zero vector
temp = find(max(abs(J)) == 0);
if ~isempty(temp)
   J(temp,:) = J(temp,:) + sqrt(eps);
end;

%calculate covariance
[Q R] = qr(J,0);
Rinv = R\eye(size(R));
diag_info = sum((Rinv.*Rinv)')';
% covariancematrix = (Rinv.*Rinv)' %% UNCOMMENT THIS LINE TO SEE THE COVARIANCE MATRIX
v = m-n;
rmse = sqrt(sum(resid.*resid)/v);

% calculate confidence interval
delta = sqrt(diag_info) .* rmse*tinv(0.975,v);
ci = [(beta(:) - delta) (beta(:) + delta)];

%--end of nlparci.m---

