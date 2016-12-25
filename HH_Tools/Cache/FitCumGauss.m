function para = FitCumGauss(x,y,w)
% [bias,sigma] = FitCumGauss(x,y,w)
% This function fits the psychometric curve using cumulative gauss function
% by HH 2012/10/20

if nargin < 3
    w = ones(size(x));
end

options = optimset('MaxIter',5000,'MaxFunEvals',5000);
errorCumGauss = inline('norm((normcdf(x,para(1),para(2))-y).*w)','para','x','y','w');  % Sum square error

% Initial guess to avoid divergence
biasGuess = [-10 -1 0 1 10];
sigmaGuess = [0.5 1 2 4 8 16 32 100];
errorGuess = zeros(length(biasGuess),length(sigmaGuess));
for i = 1:length(biasGuess)
    for j = 1:length(sigmaGuess)
        errorGuess(i,j) = errorCumGauss([biasGuess(i),sigmaGuess(j)],x,y,w);
    end
end
[minId1,minId2] = find(errorGuess == min(errorGuess(:)));
bias0 = biasGuess(minId1);
sigma0 = sigmaGuess(minId2);

% Do the fitting
[para,~,flag] = fminsearch(@(para)errorCumGauss(para,x,y,w),[bias0,sigma0],options);

% Deal with exceptions
if flag ~= 1 || para(2) > 100000
    para = [NaN NaN];
end
end
