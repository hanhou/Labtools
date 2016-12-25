% bootstrap_2afc.m
%
% This file illustrates the analysis of 2afc data using constrained
% maximum likelihood and bootstrap analyzes in MATLAB.
%
% The analysis is explained in more detail in the "Quick start to psignifit"
% that can be found at http://psignifit.sourceforge.net/TUTORIAL.html


% Set the priors
priors.m_or_a = 'None';
priors.w_or_b = 'None';
priors.lambda = 'Uniform(0,0.1)';

% Type of data
nafc = 2;

% Define the data
stimulus_intensities = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
number_of_correct = [34, 32, 40, 48, 50, 48];
number_of_trials  = [50, 50, 50, 50, 50, 50];
data = [stimulus_intensities; number_of_correct; number_of_trials]';

% Run the bootstrap inference
results = BootstrapInference(data, priors);

% Print thresholds and slops for all cuts
for i = 1:length(results.cuts)
    th = sprintf('Threshold(%.2f)  = \t %f', results.cuts(i), getThres(results, i));
    disp(th)
    sl = sprintf('Slope(%.2f) \t = \t %f', results.cuts(i), getSlope(results, i));
    disp(sl)
end

% Plot the results
plotPMF(results)
