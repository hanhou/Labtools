function results = MapEstimate ( data, priors, varargin )
% results = MapEstimate ( data, priors, ... )
%
% Determine a point estimate of the psychometric function.
%
% The point estimate is taken to be the Maximum A Posteriori in a bayesian setting.
% If all priors are 'None', this is the same as the maximum likelihood estimator.
% For proper priors, the MAP estimator takes into account the constraints imposed
% by the priors.
%
% data should be an array with three columns and one row for each block of trials.
%    The first column should contain the stimulus intensity in the respective block,
%    the second column should contain the number of correctly identified trials in
%    the respective block, and the third column should contain the total number of
%    trials presented in the respective block.
%
% priors should be a struct with the fields m_or_a, w_or_b, lambda, and gamma. 'None'
%    can be used to specify "No prior" i.e. an improper, flat prior. A valid prior would
%    be
%
%    >> priors.m_or_a = 'None';
%    >> priors.w_or_b = 'None';
%    >> priors.lambda = 'Uniform(0,.1)';
%    >> prior.gammas  = 'Uniform(0,.1)';
%
%    For more information on the specification of priors for psychometric functions, see
%
%    http://psignifit.sourceforge.net/BAYESINTRO.html#specification-of-prior-distributions
%
% The behavior of the function can be modified by setting a number of parameters for the fitting
% procedure. To specify one of these parameters, call the function as
%
% >> MapEstimate ( data, priors, 'PARAM', VALUE )
%
% Some parameters don't require a value (for instance the 'gammaislambda' setting). In this
% case call the function as
%
% >> MapEstimate ( data, priors, 'PARAM' )
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
verbosity = '';
cuts = [0.25,0.5,0.75];
verbose = false;

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
        gil = '-gammaislambda';
        gammaislambda = true;
    case 'verbose'
        % verbosity = '-v';   % Matlab system call merges stdout and stderr --- so this does not well
        verbose = true;
    case 'cuts'
        [cuts,varargin] = popoption(varargin);
    otherwise
        warning ( sprintf ( 'unknown option: %s !\n' , opt ) );
    end
end

% Store the data
dataf = tempname;
save ( '-ascii', dataf, 'data' );

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
cmd = sprintf ( 'psignifit-mapestimate %s --matlab -prior1 "%s" -prior2 "%s" -prior3 "%s" %s -nafc %d -s %s -c %s -cuts %s %s', ...
    dataf, getfield(priors,'m_or_a'), getfield(priors,'w_or_b'), getfield(priors,'lambda'), prior4, ...
    nafc, sigmoid, core, scuts, verbosity);

if verbose
    cmd
end

% Do the real work
[status,output] = system ( cmd );
eval ( output );

results.call = 'mapestimate';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;
results.priors = priors;
results.burnin = 1;
results.nsamples = 0;

% clean up
delete ( dataf );
