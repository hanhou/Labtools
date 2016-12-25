% bootweib takes in the bootstrap samples, and refits the sample to a
% weibull_bs_fit.  The weibull threshold will be returned.
function [thresh alpha beta gamma] = bootweib(bstrp_trials);

global Data;
TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

coherence = Data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence(bstrp_trials)');

for i = 1:length(unique_coherence)
    ok_values{i} = bstrp_trials(find (coherence(bstrp_trials) == unique_coherence(i)));
    c = sum(Data.misc_params(OUTCOME, ok_values{i}));
    t = length(ok_values{i});
    fit_data(i,1) = unique_coherence(i);
    fit_data(i,2) = c/t;
    fit_data(i,3) = t;
end

[alpha beta gamma] = weibull_bs_fit(fit_data);
thresh = weibull_bs_threshold([alpha beta gamma]);

return;

