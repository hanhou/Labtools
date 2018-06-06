% to evaluate how well the models fit neural responses
% BIC (Bayesian Information Criterion)
% RSS: Residual sum of squares
% n: number of data points
% p: number of function parameters
% the lower the BIC value, the model is better
% from ( Laurens et al., 2017, elife ; Konishi S and Kitagawa G. 2008.)
% LBY 20170328

function BIC =  BIC_fit(n,RSS,p)

BIC = n*log(RSS/n) + p*log(n);

end