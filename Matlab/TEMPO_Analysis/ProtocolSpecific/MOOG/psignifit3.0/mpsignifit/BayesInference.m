function results = BayesInference ( data, priors, varargin )
% results = BayesInference ( data, priors, ... )
%
% Bayesian inference for psychometric functions using Makov Chain Monte Carlo.
%
% Data should be an array with three columns and one row for each block of trials.
%    The first column should contain the stimulus intensity in the respective block,
%    the second column should contain the number of correctly identified trials in
%    the respective block, and the third column should contain the total number of
%    trials presented in the respective block.
%
%    A valid data variable would be:
%
%    >> data = [0, 1, 2, 3; 3, 5, 9; 10, 10, 10]''
%
%
% Priors should be a struct with the fields m_or_a, w_or_b, lambda, and gamma. 'None'
%    can be used to specify "No prior" i.e. an improper, flat prior. A valid prior would
%    be:
%
%    >> priors.m_or_a = 'None';
%    >> priors.w_or_b = 'None';
%    >> priors.lambda = 'Uniform(0,.1)';
%    >> priors.gamma  = 'Uniform(0,.1)';
%
%    For more information on the specification of priors for psychometric functions, see
%
%    http://psignifit.sourceforge.net/BAYESINTRO.html#specification-of-prior-distributions
%
% In the default configuration, the function draws samples from the posterior distribution
% of model parameters given the data using the Metropolis-Hastings algorithm. The following
% parameters can be used to modify the behavior of the function. To specify one of these
% parameters, call the function as
%
% >> BayesInference ( data, priors, 'PARAM', VALUE )
%
% Some parameters don't require a value (for instance the 'gammaislambda' setting). In this
% case call the function as
%
% >> BayesInference ( data, priors, 'PARAM' )
%
% In the folloing list, all valid parameters are given. If a parameter requires a value argument,
% the type of this value argument is given after the name of the parameter and separated by a comma.
%
% Parameters
% ----------
%
% 'nafc', integer               Default: 2
%       number of alternatives in the analyzed task. If nafc is 1 this means that a yes-no task was used.
%       In this case, the lower asymptote is considered a free parameter and is fitted as a fourth parameter.
%       See gammaislambda below too.
%
% 'sigmoid', char               Default: 'logistic'
%       Name of the sigmoid object to be used for fitting. psignifit includes a large number of sigmoids.
%       These include (but are not restricted to) 'logistic', 'gauss', or 'gumbel'. See
%
%       http://psignifit.sourceforge.net/MODELSPECIFICATION.html#specifiing-the-shape-of-the-psychometric-function
%
%       for further information.
%
% 'core', char                  Default: 'mw0.1'
%       Name of the core object to be used for fitting. psignifit includes a large number of cores. These
%       include (but are not restricted to) 'ab', 'mw0.1'. See
%
%       http://psignifit.sourceforge.net/MODELSPECIFICATION.html#specifiing-the-shape-of-the-psychometric-function
%
%       for further information.
%
% 'gammaislambda'               Default: false
%       In yes-no tasks, it is sometimes of interest to constrain the model to have the same upper and lower
%       asymptote. This can be achieved in psignifit by constraining gamma to be lambda, effectively reducing the
%       number of free parameters by one.
%
% 'verbose'                     Default: false
%       print out more about what is going on.
%
% 'cuts', double (vector)       Default: [0.25, 0.5, 0.75]
%       Cuts are the performances at which thresholds should be determined. Thus, setting cuts to 0.5 will define
%       the threshold to be the signal level that corresponds to a performance half way between the upper and the
%       lower asymptote.
%
% 'samples', int                Default: 2000
%       Number of MCMC samples to be generated
%
% 'pilot', double array         Default: false
%       An array of pilot samples to tune the markov chain monte carlo procedure. A good idea is often to use
%       a small number of bootstrap samples to tune the sampler, such as
%
%       >> pilot = BootstrapInference ( data, priors, 'samples', 200 );
%       >> mcmc  = BayesInference ( data, priors, 'pilot', pilot.mcestimates );
%
%       By default, NO automatic tuning is performed!
%
% 'stepwidths', double vector   Default: [0.1, 0.1, 0.01]
%       standard deviations of the different components of the proposal distribution. Also see 'pilot'
%
% 'generic'                     Default: false
%       Use the generic Metropolis Hastings algorithm proposed by Raftery and Lewis. This is often a good choice if you
%       use a pilot sample.
%
%
% This file is part of psignifit3 for matlab. (c) 2010 by Ingo Fründ

% Check data format
if size ( data, 2 ) ~= 3
    error ( 'data should have three columns' );
end

% default values
nafc = 2;
sigmoid = 'logistic';
core    = 'mw0.1';
gil = '';
gammaislambda = false;
verbosity = '';
cuts = [0.25,0.5,0.75];
samples = 2000;
verbose = false;
stepwidths = [0.1,0.1,0.01];
pilot = false;
generic = '';

% If priors for lambda and gamma are set as 'None', display warning
if exist ( 'priors.lambda' ) 
    if priors.lambda == 'None'
        disp('No prior set on lambda; this may lead to non-plausible value for the lapse rate')
    end
end

if exist ( 'priors.gamma')
    if priors.gamma == 'None'
        disp('No prior set on gamma; this may lead to non-plausible value for the guess rate')
    end
end

% Set a default prior if none is present
if exist ( 'priors' ) ~= 1;
    priors.m_or_a = 'None';
    priors.w_or_b = 'None';
    priors.lambda = 'Uniform(0,.1)';
    priors.gamma  = 'Uniform(0,.1)';
end


% Check input
while size(varargin,2) > 0
    [opt,varargin] = popoption ( varargin );
    switch opt
    case 'nafc'
        [nafc,varargin] = popoption(varargin);
    case 'sigmoid'
        [sigmoid,varargin] = popoption(varargin);
    case 'core'
        [core,varargin] = popoption(varargin);
    case 'gammaislambda'
        gil = '-e';
        gammaislambda = true;
    case 'verbose'
        % verbosity = '-v';   % Matlab system call merges stdout and stderr --- so this does not well
        verbose = true;
    case 'cuts'
        [cuts,varargin] = popoption(varargin);
    case 'samples'
        [samples,varargin] = popoption(varargin);
    case 'stepwidths'
        [stepwidths,varargin] = popoption(varargin);
    case 'pilot'
        pilot = true;
        pilotf = tempname;
        [pilotsample,varargin] = popoption(varargin);
        f = fopen ( pilotf, 'w' );
        fprintf ( f, '\n# mcestimates\n' );
        for k = 1:size(pilotsample,1)
            fprintf ( f, '%s\n', num2str(pilotsample(k,:)) );
        end
        fclose(f);
    case 'generic'
        generic = '-generic';
    otherwise
        warning ( sprintf ( 'unknown option: %s !\n' , opt ) );
    end
end

if pilot
    stepwidths_or_pilot = pilotf;
else
    stepwidths_or_pilot = sprintf('"%s',num2str(stepwidths(:)','%f,'));
    stepwidths_or_pilot(end) = '"';
end

% Store the data
dataf = tempname;
save ( '-ascii', dataf, 'data' );

% Fiddle around with the fourth prior. Do we need it?
if nafc > 1
    prior4 = '';
elseif gammaislambda
    prior4 = '';
else
    prior4 = sprintf ( '-prior4 "%s"', getfield ( priors, 'gamma' ) );
end

% Determine cuts
scuts = sprintf ( '"%s', num2str ( cuts, '%f,') );
scuts(length(scuts)) = '"';

% Write the command
cmd = sprintf ( 'psignifit-mcmc %s %s --matlab -prior1 "%s" -prior2 "%s" -prior3 "%s" %s -nsamples %d -nafc %d -s %s -c %s %s -cuts %s -proposal %s %s', ...
    verbosity, dataf, ...
    getfield(priors,'m_or_a'), getfield(priors,'w_or_b'), getfield(priors,'lambda'), prior4, ...
    samples, nafc, sigmoid, core, gil, scuts, stepwidths_or_pilot, generic );

if verbose
    cmd
end

% Do the real work
[status,output] = system ( cmd );
eval ( output );

% Store paradigm
results.call = 'bayes';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;
results.priors = priors;
results.params_estimate = mean ( results.mcestimates(samples/2:end,:) );
results.burnin = samples/2;
results.nsamples = samples;

% Clean up
delete ( dataf );
if exist ( stepwidths_or_pilot, 'file' )
    delete ( stepwidths_or_pilot );
end
