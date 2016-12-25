function results = BootstrapInference ( data, priors, varargin )
% results = BootstrapInference ( data, constraints, ... )
%
% Determine bias corrected and accelerated bootstrap confidence intervals for a psychometric function
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
% Constraints should be a struct with the fields m_or_a, w_or_b, lambda, and gamma. 'None'
%    can be used to specify "No constraint". Valid constraints would
%    be:
%
%    >> priors.m_or_a = 'None';
%    >> priors.w_or_b = 'None';
%    >> priors.lambda = 'Uniform(0,.1)';
%    >> priors.gamma  = 'Uniform(0,.1)';
%
%    Internally, the constraints are used like bayesian priors. This means that constraints
%    could in principle be any probability distribution. Psignifit implements a number of
%    probability distributions to be used as constraints / priors. For more information on
%    the specification of constraints / priors for psychometric functions, see
%
%    http://psignifit.sourceforge.net/BAYESINTRO.html#specification-of-prior-distributions
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
% 'nonparametric'               Default: false
%     perform nonparametric bootstrap instead of the default parametric bootstrap
%
% Examples of usage:
% ------------------
%
% results = BootstrapInference ( data, priors, 'nafc', 1, 'gammaislambda' )
%    fits a yes-no task, constraining gamma to equal lambda
%
% results = BootstrapInference ( data, priors, 'nonparametric' )
%    performs nonparametric bootstrap
%
% This function is part of psignifit3 for matlab (c) 2010 by Ingo Fründ

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
npr = '';
verbosity = '';
cuts = [0.25,0.5,0.75];
samples = 2000;
verbose = false;

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
    case 'nonparametric'
        npr = '-nonparametric';
    case 'verbose'
        % verbosity = '-v';   % Matlab system call merges stdout and stderr --- so this does not well
        verbose = true;
    case 'cuts'
        [cuts,varargin] = popoption(varargin);
    case 'samples'
        [samples,varargin] = popoption(varargin);
    otherwise
        warning ( sprintf ( 'unknown option: %s !\n' , opt ) );
    end
end

% Get the point estimate
if gammaislambda
    mapest = MapEstimate ( data, priors, 'nafc', nafc, 'sigmoid', sigmoid, 'core', core, 'cuts', cuts, 'gammaislambda' );
else
    mapest = MapEstimate ( data, priors, 'nafc', nafc, 'sigmoid', sigmoid, 'core', core, 'cuts', cuts  );
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
cmd = sprintf ( 'psignifit-bootstrap %s %s --matlab -prior1 "%s" -prior2 "%s" -prior3 "%s" %s -nsamples %d -nafc %d -s %s -c %s %s %s -cuts %s', ...
    verbosity, dataf, ...
    getfield(priors,'m_or_a'), getfield(priors,'w_or_b'), getfield(priors,'lambda'), prior4, ...
    samples, nafc, sigmoid, core, gil, npr, scuts );

if verbose
    cmd
end

% Do the real work
[status,output] = system ( cmd );
eval ( output );

% Store paradigm
results.call = 'bootstrap';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;
results.priors = priors;
results.params_estimate = mapest.params_estimate;
results.mapest = mapest;
results.burnin = 1;
results.nsamples = samples;

% Clean up
delete ( dataf );
