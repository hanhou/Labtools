function createFit(x,t)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(X,T)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1


% Data from dataset "t vs. x":
%    X = x:
%    Y = t:
%    Unweighted
%
% This function was automatically generated on 10-Sep-2013 02:06:00

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[737 271 688 485]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;


% --- Plot data originally in dataset "t vs. x"
x = x(:);
t = t(:);
h_ = line(x,t,'Parent',ax_,'Color',[0.333333 0 0.666667],...
    'LineStyle','-', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(x));
xlim_(2) = max(xlim_(2),max(x));
legh_(end+1) = h_;
legt_{end+1} = 't vs. x';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
    set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-0.19412532624559101, 5.8141706237448787]);
end


% --- Create fit "fit 1"
fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Algorithm','Levenberg-Marquardt','MaxFunEvals',5000);
ok_ = isfinite(x) & isfinite(t);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [30 10 0.47999999999999998 0.60599999999999998 ];
set(fo_,'Startpoint',st_);
ft_ = fittype('a*exp(-2*(1-cos(x-pref))/sigma^2)+b',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b', 'pref', 'sigma'});

% Fit this model using new data
cf_ = fit(x(ok_),t(ok_),ft_,fo_)

% Or use coefficients from the original fit:
if 0
    cv_ = { -47.598632596144483, 52.239567166113389, -0.3547186384542293, -2.4477034573613983};
    cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';

% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
