function [a,aerr,chisq,yfit] = fitnonlin(x,y,sig,fitfun, a0)

% FITNONLIN Fit a nonlinear function to data using the gradient search 
%           method discussed in Bevington and Robinson in Section 8.4.
%    [a,aerr,chisq,yfit] = fitnonlin(x,y,sig,fitfun,a0) 
%
%    Inputs:  x -- the x data to fit
%             y -- the y data to fit
%             sig -- the uncertainties on the data points
%             fitfun -- the name of the function to fit to
%             a0 -- the initial guess at the parameters 
%
%    Outputs: a -- the best fit parameters
%             aerr -- the errors on these parameters
%             chisq -- the final value of chi-squared
%             yfit -- the value of the fitted function
%                     at the points in x
%  
%    Note: "fitfun" should be in a .m file similar to the
%    following example.
%
%          The following lines are saved in a file called
%          "sinfit.m", and the routine is invoked with
%          the fitfun parameter equal to 'sinfit' (including
%          the quotes)
%
%          function f = sinfit(x,a)
%          f = a(1)*sin(a(2)*x+a(3));
%

%*****************************************
%*** Parameters you may need to modify ***
%*****************************************
stepsize = [.1 .01 .01];
%abs(a0)*.01; % the amount parameters will be varied by in each iteration
chicut = .0000000001;  % Maximum differential allowed between successive chi^2 values
silent = 0;
%*****************************************
% These parameters can be varied if you have reason to believe your fit is
% converging to quickly or that you are in a local minima of the chi
% square.

stepdown = .5;
a = a0;

chi2 = calcchi2(x,y,sig,fitfun,a);
chi1 = chi2+chicut*2;

% keep looking while the value of chi^2 is changing
i=0;
while (abs(chi2-chi1))>chicut
    i=i+1;
    % Unless silent=1, the following is printed for each iteration:
    %   the current best fit parameters "a"
    %   the current chi-square "ChiSq"
    %   the change in the chi-square "diff"
  if silent~=1
    fprintf(1,'a = ');
    fprintf(1,'%f ',a);
    fprintf(1,'\t ChiSq = %f', chi2);
    fprintf(1,'\t diff = %f\n', abs(chi2-chi1));
    fprintf(1,'\n');
  end
 
  [anew,stepsum] = gradstep(x,y,sig,fitfun,a,stepsize,stepdown);
  a = anew;
  stepdown = stepsum;
  chi1 = chi2;
  chi2 = calcchi2(x,y,sig,fitfun,a);
end
i
%Unless silent=1, prints the last values after the minimum chi-sq has been
%found
if silent~=1
    fprintf(1,'a = ');
    fprintf(1,'%f ',a);
    fprintf(1,'\t ChiSq = %f', chi2);
    fprintf(1,'\t diff = %f\n', abs(chi2-chi1));
end


% calculate the returned values
aerr = sigparab(x,y,sig,fitfun,a,stepsize);
chisq = calcchi2(x,y,sig,fitfun,a);
yfit = feval(fitfun,x,a);

%------------------- 
% the following function calculates the (negative) chi^2 gradient at
% the current point in parameter space, and moves in that direction
% until a minimum is found returns the new value of the parameters and 
% the total length travelled

function [anew,stepsum] = gradstep(x,y,sig,fitfun,a,stepsize, stepdown)

chi2 = calcchi2(x,y,sig,fitfun,a);

grad = calcgrad(x,y,sig,fitfun,a,stepsize);

chi3 = chi2*1.1;
chi1 = chi3;

% cut down the step size in parameter space until a single step in 
% the direction of the negative gradient yields a decrease in chi^2
stepdown = stepdown*2;
while chi3>chi2
  stepdown = stepdown/2;
  anew = a+stepdown*grad;
  chi3 = calcchi2(x,y,sig,fitfun,anew);
end
stepsum = 0;

% keep going in this direction until a minimum is passed
while chi3<chi2
  stepsum = stepsum+stepdown;
  chi1 = chi2;
  chi2 = chi3;
  anew = anew+stepdown*grad;
  chi3 = calcchi2(x,y,sig,fitfun,anew);
end
% approximate the minimum as a parabola (See Bevington p. 147).
step1 = stepdown*((chi3-chi2)/(chi1-2*chi2+chi3)+.5);
anew = anew - step1*grad;

%--------------------
% this function just calculates the value of chi^2
function chi2 = calcchi2(x,y,sig,fitfun,a)

chi2 = sum( ((y-feval(fitfun,x,a)) ./sig).^2);


%--------------------
% this function calculates the (negative) gradient at a point in 
% parameter space (See Bevington p. 154).
function grad = calcgrad(x,y,sig,fitfun,a, stepsize)

f = 0.01;
[dum, nparm] = size(a);

grad = a;
chisq2 = calcchi2(x,y,sig,fitfun,a);  

for i=1:nparm
  a2 = a;
  da = f*stepsize(i);
  a2(i) = a2(i)+da;
  chisq1 = calcchi2(x,y,sig,fitfun,a2);
  grad(i) = (chisq2-chisq1);
end

t = sum(grad.^2);
grad = stepsize.*grad/sqrt(t);

%-----------------
% this function calculates the errors on the final fitted 
% parameters by approximating the minimum as parabolic
% in each parameter (See Bevington, p. 147).
function err=sigparab(x,y,sig,fitfun,a,stepsize)

[dum, nparm] = size(a);

err = a;
chisq2 = calcchi2(x,y,sig,fitfun,a);  

for i=1:nparm

  a2 = a;
  da = stepsize(i);
  a2(i) = a2(i)+da;
  chisq3 = calcchi2(x,y,sig,fitfun,a2);
  a2(i) = a2(i)-2*da;
  chisq1 = calcchi2(x,y,sig,fitfun,a2);
  err(i)=da*sqrt(2/(chisq1-2*chisq2+chisq3));

end





