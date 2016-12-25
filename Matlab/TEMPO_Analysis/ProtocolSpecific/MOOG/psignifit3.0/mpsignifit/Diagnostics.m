function results = Diagnostics ( data, parameters, varargin )
% results = Diagnostics ( data, parameters, ... )
%
% Determine some parameters of fitted psychometric function parameters.
% This involves evaluation of the function for different x values or determining
% some diagnostic values for goodness of fit.
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
% Parameters should a vector of parameters at which the diagnostic values should
%    be evaluated.
%
% The behavior of this function can be modified by specifiing certain addition
% parameters. To specify one of these parameters, call the function as
%
% >> Diagnostics ( data, priors, 'PARAM', VALUE )
%
% Some parameters do not require a value (for instance the 'gammaislambda' setting). In this
% case call the function as
%
% >> Diagnostics ( data, priors, 'PARAM' )
%
% In the following list, all valid parameters are given. If a parameter requires a value argument,
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
    otherwise
        warning ( sprintf ( 'unknown option: %s !\n' , opt ) );
    end
end

% Store the data
dataf = tempname;
save ( '-ascii', dataf, 'data' );

sparams = sprintf ( '"%s', num2str ( parameters, '%f,') );
sparams(end) = '"';

% Determine cuts
scuts = sprintf ( '"%s', num2str ( cuts, '%f,') );
scuts(end) = '"';

% Write the command
cmd = sprintf ( 'psignifit-diagnostics %s --matlab -c %s -s %s -params %s -cuts %s -nafc %d %s %s', ...
    dataf, core, sigmoid,sparams,scuts,nafc,verbosity, gil );

if verbose
    cmd
end

[status,output] = system ( cmd );
eval ( output );

results.call = 'diagnostics';
results.nafc = nafc;
results.sigmoid = sigmoid;
results.core = core;
results.gammaislambda = gammaislambda;
results.cuts = cuts;
results.data = data;

% Clean up
delete ( dataf );
