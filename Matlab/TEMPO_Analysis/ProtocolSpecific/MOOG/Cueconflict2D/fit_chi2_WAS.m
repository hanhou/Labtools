function  y = fit_chi2( x , p )
% This function is a helper for doing chi2 for conflict fitting.
y = p.func( x , p.weights );