function ci = getCI ( inference, cut, p )
% ci = getCI ( inference, cut )
%
% Get the confidence interval from the inference object
%
% inference should be a struct as generated either by BayesInference or BootstrapInference
%
% cut should be the index of the desired cut
%
% p should be the desired coverage of the confidence interval (p<1)
%
%
% This file is part of psignifit3 for matlab (c) by Ingo Fründ

notin = 1-p;
probs = [0.5*notin,1-0.5*notin];

if strcmp ( inference.call, 'bootstrap' )
    bias = inference.bias_thres(cut);
    acc  = inference.acc_thres(cut);
    probs = normcdf( bias + ( norminv(probs) + bias ) ./ (1-acc*(norminv(probs) + bias )) );
end;

ci = prctile ( inference.mcthres(:,cut), 100*probs );
