% This file performs some fake unit tests of the matlab interface.
% Essentially this only tests whether all functions can be called

setuptestenvironment

% 2 afc
disp ( 'testing 2afc setting' );
disp ( 'MapEstimate' );
mapest = MapEstimate ( data2afc, priors );
disp ( 'BootstrapInference' );
boots  = BootstrapInference ( data2afc, priors, 'samples', 200 );
disp ( 'BayesInference: explicite stepwidths' );
mcmc   = BayesInference ( data2afc, priors, 'samples', 200, 'stepwidths', std(boots.mcestimates,1) );
disp ( 'BayesInference: pilot sample' );
mcmcp  = BayesInference ( data2afc, priors, 'samples', 200, 'pilot', boots.mcestimates );
disp ( 'BayesInference: generic' );
mcmcg  = BayesInference ( data2afc, priors, 'samples', 200, 'pilot', boots.mcestimates, 'generic' );

% GoodnessOfFit calls Diagnostics internally
disp ( 'Goodness of fit:' );
disp ( 'Bootstrap');
GoodnessOfFit ( boots );
disp ( 'mcmc' );
GoodnessOfFit ( mcmc );
disp ( 'mcmc(pilot)');
GoodnessOfFit ( mcmcp );
disp ( 'mcmc(generic)' );
GoodnessOfFit ( mcmcg );
disp ( 'clean up');
close all;

% 1 afc
disp ( '=================================================================================')
disp ( 'testing yes no setting' );
disp ( 'MapEstimate' );
mapest = MapEstimate ( data1afc, priors, 'nafc', 1 );
disp ( 'BootstrapInference' );
boots  = BootstrapInference ( data1afc, priors, 'nafc', 1, 'samples', 200 );
disp ( 'BayesInference: explicite stepwidths' );
mcmc   = BayesInference ( data1afc, priors, 'nafc', 1, 'samples', 200, 'stepwidths', std(boots.mcestimates,1) );
disp ( 'BayesInference: pilot sample' );
mcmcp  = BayesInference ( data1afc, priors, 'nafc', 1, 'samples', 200, 'pilot', boots.mcestimates );
disp ( 'BayesInference: generic' );
mcmcg  = BayesInference ( data1afc, priors, 'nafc', 1, 'samples', 200, 'pilot', boots.mcestimates, 'generic' );

% GoodnessOfFit calls Diagnostics internally
disp ( 'Goodness of fit:' );
disp ( 'Bootstrap');
GoodnessOfFit ( boots );
disp ( 'mcmc' );
GoodnessOfFit ( mcmc );
disp ( 'mcmc(pilot)');
GoodnessOfFit ( mcmcp );
disp ( 'mcmc(generic)' );
GoodnessOfFit ( mcmcg );
disp ( 'clean up');
close all;

% 1 afc with gamma==lambda
disp ( '=================================================================================')
disp ( 'testing yes no with gamma==lambda' );
disp ( 'MapEstimate' );
mapest = MapEstimate ( data1afc, priors, 'nafc', 1, 'gammaislambda' );
disp ( 'BootstrapInference' );
boots  = BootstrapInference ( data1afc, priors, 'nafc', 1, 'gammaislambda', 'samples', 200 );
disp ( 'BayesInference: explicite stepwidths' );
mcmc   = BayesInference ( data1afc, priors, 'nafc', 1, 'gammaislambda', 'samples', 200, 'stepwidths', std(boots.mcestimates,1) );
disp ( 'BayesInference: pilot sample' );
mcmcp  = BayesInference ( data1afc, priors, 'nafc', 1, 'gammaislambda', 'samples', 200, 'pilot', boots.mcestimates );
disp ( 'BayesInference: generic' );
mcmcg  = BayesInference ( data1afc, priors, 'nafc', 1, 'gammaislambda', 'samples', 200, 'pilot', boots.mcestimates, 'generic' );

% GoodnessOfFit calls Diagnostics internally
disp ( 'Goodness of fit:' );
disp ( 'Bootstrap');
GoodnessOfFit ( boots );
disp ( 'mcmc' );
GoodnessOfFit ( mcmc );
disp ( 'mcmc(pilot)');
GoodnessOfFit ( mcmcp );
disp ( 'mcmc(generic)' );
GoodnessOfFit ( mcmcg );
disp ( 'clean up');
close all;
