function axhandle = plotPMF ( inference, varargin )
% axes = PlotPMF ( inference, ... )
%
% Plots an inference object as given by one of the commands BayesInference, BootstrapInference, or MapEstimate
%
% Additional Parameters can be set by calling
%
% axes = PlotPMF ( inference, 'PARAMETER', 'VALUE', ... )
%
% Parameters
% ----------
%
% verbose       show many status messages
% axes          axes handle where the plot should go
% color         use this color for the plot
% xlabel        label for the x axis
% ylabel        label for the y axis
% 
% This file is part of psignifit3 for matlab (c) by Ingo FrÃ¼nd

% Check data format
if size ( inference.data, 2 ) ~= 3
    error ( 'data should have three columns' );
end

if ~isstruct ( inference )
    error ( 'inference should be a struct' );
end

% Set default parameters
verbose = false;
axhandle = gca;
color = [0,0,1];
xname = 'Stimulus intensity';
yname = 'Proportion correct';

% posterior samples only make sense in a bayesian framework
if strcmp(inference.call, 'bayes')
    nposteriorsamples = 20;
elseif strcmp(inference.call, 'bootstrap')
    nposteriorsamples = 0;
elseif strcmp(inference.call, 'mapestimate')
    nposteriorsamples = 0;
else
    error ( 'I have no idea what kind of inference was performed' );
end

% Check input
while size(varargin,2) > 0
    [opt,varargin] = popoption ( varargin );
    switch opt
    case 'verbose'
        verbose = true;
    case 'axes'
        [axhandle,varargin] = popoption ( varargin );
    case 'color'
        [color,varargin] = popoption ( varargin );
    case 'xlabel'
        [xname,varargin] = popoption ( varargin );
    case 'ylabel'
        [yname,varargin] = popoption ( varargin );
    otherwise
        warning ( sprintf ( 'unknown option: %s !\n' , opt ) );
    end
end

% Range in which the deviance should be colored
burnin = inference.burnin;
drange = [min(inference.mcdeviance(burnin:end)),min(inference.mcdeviance(burnin:end))+2*std(inference.mcdeviance(burnin:end))];

% Diagnostics of the point estimate
if inference.gammaislambda
    diag = Diagnostics ( inference.data, inference.params_estimate, ...
        'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts, 'gammaislambda' );
else
    diag = Diagnostics ( inference.data, inference.params_estimate, ...
        'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts );
end


% Prepare the plot
cla ( axhandle );
hold on;

% If we are in a bayesian framework, we want a couple of examples from the posterior
for k = 1:nposteriorsamples
    whichpmf = floor(inference.burnin + (inference.nsamples-inference.burnin)*rand());

    % Evaluate the sample
    if inference.gammaislambda
        pmf = Diagnostics ( inference.data, inference.mcestimates ( whichpmf, : ), ...
            'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts, 'gammaislambda' );
    else
        pmf = Diagnostics ( inference.data, inference.mcestimates ( whichpmf, : ), ...
            'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts );
    end

    % Determine appropriate color
    if pmf.deviance>drange(2)
        deviance = drange(2);
    else
        deviance = pmf.deviance;
    end
    dpos = (deviance-drange(1))./(drange(2)-drange(1));
    if isnan(dpos)
        dpos = 0;
    elseif dpos < 0
        dpos = 0;
    elseif dpos > 1
        dpos = 1;
    end

    % Plot the example
    plot ( axhandle, pmf.pmf(:,1), pmf.pmf(:,2), '-', 'color', dpos*[1,1,1]+(1-dpos)*color, 'linewidth', 1 );
end

% In any case we need data and the pointestimate
plot ( axhandle, inference.data(:,1), inference.data(:,2)./inference.data(:,3), 'o', 'color', color );
plot ( axhandle, diag.pmf(:,1), diag.pmf(:,2), '-', 'color', color, 'linewidth', 3 );

% And we need the confidence intervals
for cut = 1:length(inference.cuts)
    ci = getCI ( inference, cut, 0.95 );
    if inference.nafc>1
        guess = 1./inference.nafc;
    else
        guess = inference.params_estimate(end);
    end
    h = guess + inference.cuts(cut)*(1-guess-inference.params_estimate(3));
    plot ( axhandle, ci, [h,h], '-', 'color', color );
end

% Set the ranges
xmin = min ( inference.data(:,1) );
xmax = max ( inference.data(:,1) );
xrange = xmax-xmin;
set(axhandle, 'xlim', [xmin-.05*xrange,xmax+.05*xrange], 'ylim', [0,1] );

% Label axes
xlabel ( xname );
ylabel ( yname );


% Add info about psychometric function as multi-line text object
textstr(1) = {['sigmoid: ', inference.sigmoid]};
textstr(2) = {['core: ', inference.core]};
textstr(3) = {['nAFC: ', num2str(inference.nafc)]};
textstr(4) = {['Deviance: ', num2str(diag.deviance)]};

text(xmin, 0.9, textstr, 'HorizontalAlignment', 'left');

% Releast the plot
hold off;

