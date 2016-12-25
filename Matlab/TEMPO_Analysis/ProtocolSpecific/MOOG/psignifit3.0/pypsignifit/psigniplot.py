#/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__docformat__ = "restructuredtext"

import pylab as p
import numpy as N
import re
from scipy import stats
import psignidata
import swignifit.interface_methods as interface

from matplotlib.ticker import MaxNLocator

__all__ = ["GoodnessOfFit","ConvergenceMCMC","ParameterPlot","ThresholdPlot","plotSensitivity","plotInfluential","plotMultiplePMFs"]
__warnred = [.7,0,0]

import warnings
spineswarning = """your axes object does not support spines.
The most probable reason for this is that you are using an old version of matplotlib.
spines allow for more beautiful plots and are a new feature in matplotlib 1.0.0.
"""

class parameterdict ( dict ):
    def __add__ ( self, other ):
        out = parameterdict ( self )
        for k,i in other.iteritems ():
            out.setdefault ( k, i )
        return out

class DefaultParameters ( object ):
    def __init__ ( self ):
        self.alltext = parameterdict()
        self.text    = parameterdict(fontsize=8)
        self.title   = parameterdict(fontsize=12)
        self.label   = parameterdict(fontsize=10)
        self.allplots= parameterdict(color='b')
        self.line    = parameterdict()
        self.highlight = parameterdict(color='r')
        self.warning = parameterdict(fontsize=16, color=[.7,0,0])

rc = DefaultParameters()

def drawaxes ( ax, xtics=None, xfmt=None, ytics=None, yfmt=None, xname=None, yname=None ):
    """Draw x and y axes that look nicer than standard matplotlib

    This function deletes the pylab default axes and draws visually more
    pleasing axes. In addition it takes care that all graphics remain within
    the area covered by the axes.

    :Parameters:
        *ax* :
            the axes object in which the graphics should stay
        *xtics* :
            an array of x-ticks
        *xfmt* :
            format string for the x-ticks
        *ytics* :
            an array of y-ticks
        *yfmt* :
            format string for the y-ticks
        *xname* :
            label for the x-axis
        *yname* :
            label for the y-axis
    """

    if getattr ( ax, 'spines', False ):
        # New Implementation using spines
        for loc, spine in ax.spines.iteritems():
            if loc in ['left','bottom']:
                spine.set_position( ('outward', 10) )   # Outward by 10 points
            elif loc in ['right','top']:
                spine.set_color('none')                 # no 'spine'
            else:
                raise ValueError ( 'unknown spine location: %s' % loc )
    else:
        warnings.warn ( spineswarning, DeprecationWarning )
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_xlabel ( xname, **(rc.label+rc.alltext) )
    ax.set_ylabel ( yname, **(rc.label+rc.alltext) )

def prepare_axes ( ax, haveon=("bottom","left" ) ):
    """Prepare an axes object to look nicer than standard matplotlib

    :Parameters:
        *ax* :
            axes object that should be prepared
        *haveon* :
            axes that should be shown

    :Return:
        the prepared axes object
    """
    if getattr ( ax, 'spines', False ):
        for loc,spine in ax.spines.iteritems():
            if loc in haveon:
                spine.set_position ( ("outward",10) )
            else:
                spine.set_color ( "none" )
    else:
        warnings.warn ( spineswarning, DeprecationWarning )

    if "bottom" in haveon:
        ax.xaxis.set_ticks_position ( "bottom" )
    elif "top" in haveon:
        ax.xaxis.set_ticks_position ( "top" )
    else:
        ax.xaxis.set_ticks_position ( "none" )
        ax.xaxis.set_ticklabels ( "" )
    if "left" in haveon:
        ax.yaxis.set_ticks_position ( "left" )
    elif "right" in haveon:
        ax.yaxis.set_ticks_position ( "right" )
    else:
        ax.yaxis.set_ticks_position ( "none" )
        ax.yaxis.set_ticklabels ( "" )

    return ax

def axes_array_h ( fig, naxes, axsize, lowerleft=(0.1,0.1), dist=0.05, showally=True, nox=False ):
    """Draw a horizontal array of axes

    :Parameters:
        *fig*, matplotlib.figure instance :
            the figure in which to plot the axes
        *naxes*, integer :
            how many axes should be generated
        *axsize*, tuple:
            size of each axes system
        *lowerleft*, tuple:
            lower left corner of the first axes system
        *dist*, float:
            horizontal separation of adjacent axes
        *showally*, bool:
            indicates whether all y-axes should be shown or not
        *nox*, bool:
            if True, no x-axes are shown

    :Return:
        a sequence of the newly generated axes objects
    """
    xsize,ysize = axsize
    xdist,ydist = lowerleft
    step = xsize+dist

    if nox:
        axs = [prepare_axes ( fig.add_axes ( [xdist,ydist,xsize,ysize] ), haveon="left" )]
    else:
        axs = [prepare_axes ( fig.add_axes ( [xdist,ydist,xsize,ysize] ) )]
    for n in xrange(1,naxes):
        xdist += step
        haveon = []
        if not nox:
            haveon.append ( "bottom" )
        if showally:
            haveon.append ( "left" )
        axs.append ( prepare_axes ( fig.add_axes ( [xdist,ydist,xsize,ysize] ), haveon=haveon ) )

    return axs

def plotRd ( InferenceObject, ax=None, regressor="p" ):
    """plot deviance residuals against a regressor

    Deviance residuals are used plotted agains either predicted performance or
    block index to check for systematic deviations of the data from the fitted
    function.

    :Parameters:
        *InferenceObject* :
            a BootstrapInference or BayesInference object containing
            the actual inference data
        *ax* :
            an axes object where the plot should go
        *regressor* :
            plot deviance residuals against model prediction (p) or
            against block index (k)
    """
    if ax==None:
        ax = prepare_axes ( p.axes() )
    else:
        ax = prepare_axes ( ax )

    # Plot the data points
    if regressor=="p":
        psi = InferenceObject.evaluate ( InferenceObject.data[:,0] )
    elif regressor=="k":
        psi = N.arange(len(InferenceObject.data[:,0]))+1
    else:
        raise ValueError,"regressor %s is unknown" % regressor
    psilims = N.array([psi.min(),psi.max()])
    devianceresiduals = InferenceObject.devianceresiduals
    ax.plot ( psi, devianceresiduals, "o", color=rc.allplots['color'] )

    # Linear regression
    A = N.ones((len(psi),2),'d')
    A[:,1] = psi
    a,b = N.linalg.lstsq(A,devianceresiduals)[0]
    ax.plot(psilims,a+b*psilims,':', color=rc.allplots['color'] )

    if regressor=="p":
        if InferenceObject.model["nafc"]==1:
            ax.set_xlim (0,1)
        else:
            ax.set_xlim(1./InferenceObject.model["nafc"],1)

        # In this case predictions larger than 1 and less than 0 are impossible
        xtics = list(ax.get_xticks())
        for k in xrange(len(xtics)):
            if xtics[k]>1. or xtics[k]<0.:
                xtics.pop(k)
        ax.set_xticks(xtics)

    # Generate the respective labels
    if regressor=="p":
        ax.text(psilims.mean(),ax.get_yticks()[-2],"Rpd=%.3f" % ( InferenceObject.Rpd, ), **rc.text )
        xname = "model prediction"
    elif regressor=="k":
        ax.text(psilims.mean(),ax.get_yticks()[-2],"Rkd=%.3f" % ( InferenceObject.Rkd, ), **rc.text )
        xname = "block index"

    ax.set_ylabel ( "deviance residuals", **(rc.label+rc.text) )
    ax.set_xlabel ( xname, **(rc.label+rc.text) )

    return ax

def plotppScatter ( simdata, observed, quantity, shortname=None, ax=None ):
    """plot a scatter diagram to compare observed and predicted data

    :Parameters:
        *simdata* :
            data simulated from the model (typically data from posterior predictives)
        *obseved* :
            observed data (transformed in the same way as the posterior predictives)
        *quantity* :
            name of the quantity that is checked
        *shortname* :
            abbreviation of the plotted quantity
        *ax* :
            pylab axes object where the plot should go.
    """
    if ax==None:
        ax = p.gca()

    ax.plot ( simdata, observed, '.', color=rc.allplots['color'] )
    xl = list(ax.get_xlim())
    yl = list(ax.get_ylim())
    axmin = N.min ( list(xl)+list(yl) )
    axmax = N.max ( list(xl)+list(yl) )
    ax.plot ( [axmin,axmax],[axmin,axmax], 'k:' )
    ax.set_xlim ( axmin, axmax )
    ax.set_ylim ( axmin, axmax )

    ax.set_xlabel ( "observed "+quantity, **(rc.label+rc.alltext) )
    ax.set_xlabel ( "predicted "+quantity, **(rc.label+rc.alltext) )

    # Write diagnostics
    pval = N.mean( (simdata-observed)>=0 )
    ax.set_title ( "Bayesian p (%s)=%.3f" % (shortname,pval), **(rc.text+rc.alltext) )

    if pval<0.975 and pval>0.025:
        return True
    else:
        return False

def plotHistogram ( simdata, observed, xname, shortname=None, ax=None, hideobserved=False, reference="bootstrap" ):
    """plot a histogram and compare observed data to it

    :Parameters:
        *simdata* :
            an array of monte-carlo samples of the parameter of interest
        *observed* :
            observed value of the parameter of interest (for MCMC samples, it is often
            reasonable to use this as the value of 'no effect' or something)
        *xname* :
            name of the paramter of interest
        *shortname* :
            short name of the parameter of interest
        *ax* :
            axes object defining the area where the plot should go
        *hideobserved* :
            if this is True, the observed value is not plotted
        *reference* :
            reference of the data. Could be either a string 'bootstrap'/'bayes' or a number
            against which the histogram is tested

    :Output:
        returns a boolean value indicating whether or not the Null-Hypothesis that
            observed was drawn from the same distribution as simdata is true
    """
    if ax is None:
        ax = p.axes()

    if reference.lower()[:5]==  "boots":
        reference = observed
    elif reference.lower()[:5]=="bayes":
        reference = 0

    # Remove nan
    simdata = N.nan_to_num ( simdata )
    simdata = simdata[simdata!=0]

    # Make sure we have a useful shortname
    if shortname is None:
        shortname = xname

    # Correlations plots should be treated differently
    if shortname[0] == "R":
        ax.hist ( simdata, bins=N.arange(-1,1,.1) )
        ax.set_xlim ( -1, 1 )
    else:
        ax.hist ( simdata, bins=20 )

    # Get the tics and ranges
    # xtics = p.getp(ax,"xticks")
    # ytics = p.getp(ax,"yticks")
    # xr = xtics.max()-xtics.min()
    # yy = [ytics.min(),ytics.max()+0.02*xr]
    yy = p.array(ax.get_ylim ())
    yy[1] += 0.02*(yy[1]-yy[0])

    # Plot percentile bars
    if not hideobserved:
        c = parameterdict ( {"color": rc.highlight['color']} )
        ax.plot ( [observed]*2, yy, **(c+rc.line+rc.allplots) )
    if shortname=="D":
        p95 = p.prctile ( simdata, 95 )
        ax.plot ( [p95]*2, yy, ':', color=rc.highlight["color"] )
    else:
        p25,p975 = p.prctile ( simdata, (2.5,97.5) )
        ax.plot ( [p25]*2, yy, 'r:', [p975]*2, yy, ':', color=rc.highlight['color'] )

    # Draw the full plot
    ax.set_ylim ( yy )

    # Write diagnostics
    if shortname=="D":
        ax.set_title ( "%s=%.3f, %s_crit=%.3f" % (shortname, observed, shortname, p95 ), **(rc.text+rc.alltext) )
        if reference < p95:
            return True
        else:
            return False
    else:
        ax.set_title ( "%s=%.3f, c(2.5%%)=%.3f, c(97.5%%)=%.3f" % (shortname,observed,p25,p975), **(rc.text+rc.alltext) )

        if reference>p25 and reference<p975:
            return True
        else:
            return False

def plotPMF ( InferenceObject, xlabel_text="Stimulus intensity", ylabel_text=None,ax=None, showaxes=True, showdesc=False, **kwargs ):
    """Show the psychometric function and data in an axes system

    This function plots the best fitting psychometric function and with the
    corresponding data points. If data points are labelled influential, they
    are plotted as red squares, if data points are labelled as outliers, they
    are plotted as red triangles.
    The function uses its internal knowledge about the task (nAFC or Yes/No)
    to put the correct labels to the y-axis.

    :Parameters:
        *ax* :
            axes object in which the plot should go (default: current)
        *xlabel_text* :
            label for the x-axis
        *ylabel_text* :
            label for the y-axis, if this is None, the functions
            determines the correct label from its internal knowledge
            about the task
        *showaxes* :
            if this is set to False, no axes will be drawn
        *showdesc* :
            if this is set to False, no convergence description is drawn

    :Return:
        returns a tuple (line,points,lims)
        *line* :
            the matplotlib.lines.Line2D object representing the fitted curve
        *points* :
            the matplotlib.collections.CircleCollection object representing
            the fitted data points
        *lims* :
            limits of the drawn x axis.

    :Example:
    You can use this function to plot multiple psychometric functions. However,
    keep in mind that the function plotMultiplePMFs does the same job more
    conveniently for you. However, plotPMF will typically allow for more control
    over the plotting process. this is demonstrated below:

    >>> d0 = [[0, 28, 50], [2, 33, 50], [4, 38, 50], [6, 45, 50], [8, 45, 50], [10, 49, 50]]
    >>> d1 = [[0, 22, 50], [2, 34, 50], [4, 31, 50], [6, 42, 50], [8, 42, 50], [10, 46, 50]]
    >>> d2 = [[0, 26, 50], [2, 31, 50], [4, 38, 50], [6, 47, 50], [8, 49, 50], [10, 49, 50]]
    >>> constraints = ("","","Uniform(0,.1)")
    >>> B0 = BootstrapInference ( d0, priors=constraints, plotting={'color': 'r'} )
    >>> B1 = BootstrapInference ( d1, priors=constraints, plotting={'color': 'b'} )
    >>> B2 = BootstrapInference ( d2, priors=constraints, plotting={'color': 'k'} )
    >>> plotPMF ( B0, showaxes=False )
    >>> plotPMF ( B1, showaxes=False )
    >>> plotPMF ( B2, showaxes=True  )

    Note that the last call to plotPMF sets showaxes to True and thus draws the axes.
    """
    if ax==None:
        ax = p.gca()

    # Plot the psychometric function
    xmin = InferenceObject.data[:,0].min()
    xmax = InferenceObject.data[:,0].max()
    x = N.mgrid[xmin:xmax:100j]
    psi = InferenceObject.evaluate ( x )
    pmfline = ax.plot(x,psi,
            color     = kwargs.setdefault ( 'color', InferenceObject.color ),
            linestyle = kwargs.setdefault ( 'linestyle', InferenceObject.linestyle ),
            linewidth = kwargs.setdefault ( 'linewidth', InferenceObject.linewidth ),
            label     = kwargs.setdefault ( 'label', InferenceObject.label )
            )

    # Plot the data
    xd = InferenceObject.data[:,0]
    pd = InferenceObject.data[:,1].astype("d")/InferenceObject.data[:,2]
    nd = InferenceObject.data[:,2]
    pmfpoints = ax.scatter ( xd, pd, s=nd, c=kwargs.setdefault ( 'color', InferenceObject.color ),
            marker=kwargs.setdefault("markertype", InferenceObject.marker) )

    # Check axes limits
    ymin,ymax = -.05,1.05
    if ylabel_text is None:
        if InferenceObject.model["nafc"]>1:
            ylabel_text = "P(correct)"
        else:
            ylabel_text = "P(Yes)"

    # Determine tics
    ax.set_ylim ( ymin, ymax )
    ytics = list(ax.get_yticks())
    # Clean up ytics
    newytics = []
    while len(ytics):
        yt = ytics.pop(0)
        if yt>=0 and yt<=1.01:
            newytics.append(yt)
    ax.set_yticks ( newytics )

    # Write some model information
    if showdesc:
        txt = InferenceObject.desc
        if not InferenceObject.deviance is None:
            txt = txt+"\nD=%g" % (InferenceObject.deviance,)
        ax.text ( 0.3*(xmin+xmax),ymin+.1,txt, **(rc.text+rc.alltext) )

    return pmfline,pmfpoints,ax

def plotThres ( InferenceObject, ax=None, color="b" ):
    """Plot thresholds and confidence intervals

    :Parameters:
        *InferenceObject* :
            either a BootstrapInference object or a BayesInference object
        *ax* :
            a pylab.axes object to be used for the plot.
        *color* :
            a pylab color marker
    """
    if ax == None:
        ax = p.gca()

    # Determine the range where the data live
    datarange = InferenceObject.data[:,0].min(),InferenceObject.data[:,0].max()
    dataw = datarange[1]-datarange[0]

    for k,cut in enumerate(InferenceObject.cuts):
        c25,c975 = InferenceObject.getCI ( cut=k, conf=(.025,.975) )
        thres = InferenceObject.getThres ( cut )
        ylev = InferenceObject.evaluate ( [thres] )
        if c25 < datarange[0]-dataw*0.2:
            bar = [ datarange[0]-dataw*0.2 ]
            markers = ["<"]
            c25out = True
            ax.text ( datarange[0]-dataw*0.2,ylev,"%g"%(c25,), horizontalalignment="center", fontsize=7 )
        else:
            bar = [ c25 ]
            c25out = False
            markers = ["|"]
        if thres > datarange[0] and thres < datarange[1]:
            bar.append(thres)
            markers.append("|")
        if c975 > datarange[1]+dataw*0.2:
            bar.append(datarange[1]+dataw*0.2)
            c975out = True
            ax.text ( datarange[1]+dataw*0.2,ylev,"%g"%(c975,), horizontalalignment="center", fontsize=7 )
            markers.append(">")
        else:
            bar.append(c975)
            c975out = False
            markers.append("|")
        ax.plot ( bar,[ylev]*len(bar), '-', color=color )
        for x,m in zip(bar,markers):
            ax.plot ( [x], [ylev], marker=m, color=color )
    return ax

def GoodnessOfFit ( InferenceObject, warn=True ):
    """Draw a diagnostic figure to help assessing goodness of fit

    This graphic is intended to help the user determine how well the fitted function describes
    the data. The plot has 6 fields:

    +-----+-----+-----+
    |  1  |  3  |  5  |
    +-----+-----+-----+
    |  2  |  4  |  6  |
    +-----+-----+-----+

    The fields provide the following information:
    1.  The data and the fitted psychometric function. "fitted" here means the parameters are
        the mean of the posterior for BayesInference objects and the (constrained)
        maximum likelihood fit for BootstrapInference objects. To get an idea of the posterior
        resp. bootstrap distribution, credibility intervals are plotted at some positions (the
        location and width of the credibility intervals is given in the constructor). To make
        the posterior distribution for BayesInference objects really "plastic", a number of
        samples from the posterior distribution over psychometric functions are also drawn in
        light blue. The saturation of blue also codes the deviance of the respective function:
        the more saturated, the better the fit. For BootstrapInference objects, outliers and
        influential observations are marked as red triangles and red squares.
    2.  A histogram to approximate the posterior resp. bootstrap distribution of deviances.
        For BootstrapInference objects this histogram provides important information. It estimates
        the distribution of deviance that would be expected if the fitted model were perfectly
        valid. If the deviance of the fitted model is far in the tails of the deviance histogram,
        this typically indicates a bad fit. In that case, a warning is displayed if warn==True.
    3.  A plot of model predictions against deviance residuals. If there is no obvious
        interrelation between model prediction and deviance residuals, this indicates that the
        model describes the data reasonably well. To get an idea of the interrelation between
        model prediction and deviance residuals, the best fitting line is plotted as a dotted line.
    4.  A histogram of samples from the posterior resp. bootstrap distribution of correlations
        between model prediction and deviance residuals. The interpretation of this histogram
        differs for BootstrapInference and for BayesInference. For BayesInference the distibution
        should include 0. If the distribution is clearly shifted away from 0, this is strong
        evidence, that something is wrong with your model or your data. For BootstrapInference,
        The distribution shown corresponds to the distribution that would be expected if your
        fitted psychometric function would perfectly describe the data. Thus, if the maximum
        likelihood estimate (the vertical bold red line) is in the extremes of the distribution,
        this is strong evidence, that something is wrong with your model or your data.
    5,6 Similar to 3 and 4 but form correlations between block index and deviance residuals.
        Correlations between block index and deviance residuals indicate nonstationary
        data as should be found during e.g. perceptual learning.

    :Parameters:
        *warn* :
            if warn is set to True, red warning messages are displayed
            whenever the fit does not seem to describe the data well.
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Goodness of fit diagnostics require monte carlo samples. Try to call the sample() method of your inference object."

    fig = p.figure ( figsize=(10,8) )

    ax_plot,ax_rpd,ax_rkd       = axes_array_h ( fig, 3, (.22,.3), (.1,.6), dist=0.1, showally=True )
    ax_deviance,ax_rpdh,ax_rkdh = axes_array_h ( fig, 3, (.22,.3), (.1,.1), dist=0.1, showally=True )

    infer = InferenceObject.__repr__().split()[1]
    if infer not in ["BayesInference","BootstrapInference", "ASIRInference"]:
        raise ValueError, "Unknown InferenceObject: %s" % (InferenceObject.__repr__().split()[1],)

    # First plot about deviance
    if infer in ["BayesInference","ASIRInference"]:
        InferenceObject.drawposteriorexamples ( ax=ax_plot )
    plotThres ( InferenceObject, ax=ax_plot )
    plotPMF   ( InferenceObject, ax=ax_plot, showdesc=True )
    if infer in ["BayesInference","ASIRInference"]:
        distname = "posterior"
        observed = -2*N.log ( InferenceObject.nullevidence )
        good = plotppScatter ( InferenceObject.ppdeviance, InferenceObject.mcdeviance, "deviance", "D", ax_deviance )
    elif infer == "BootstrapInference":
        distname = "bootstrap"
        observed = InferenceObject.deviance
        good = plotHistogram ( InferenceObject.mcdeviance, observed, "bootstrap deviance", "D", ax_deviance )
    if warn and not good:
        ax_deviance.text ( N.array(ax_deviance.get_xlim()).mean(), N.array(ax_deviance.get_ylim()).mean(),
                "The fitted model is a bad\ndescription of the data!",
                horizontalalignment="center", verticalalignment="center", rotation=45, **(rc.warning+rc.alltext) )

    # The other two plots are in a loop: Rpd, Rkd
    ax = [ax_rpd,ax_rkd]
    axh = [ax_rpdh,ax_rkdh]
    index = ["p","k"]
    warningtext = ["Simulated Rpd differs from observed!\nTry other sigmoid?",
            "Simulated Rkd differs from observed!\nData are nonstationary!"]

    for k in xrange ( 2 ):
        plotRd ( InferenceObject, ax[k], index[k] )
        name = "R%sd" % (index[k],)
        if infer in ["BayesInference","ASIRInference"]:
            good = plotppScatter ( eval("InferenceObject.pp%s" % (name,)), eval("InferenceObject.mc%s"%(name,)), name,name, axh[k] )
        else:
            good = plotHistogram ( eval("InferenceObject.mc%s" % (name,)), eval("InferenceObject.%s"%(name,)), "bootstrap "+name, name, axh[k] )
        if warn and not good:
            axh[k].text ( p.array(axh[k].get_xlim()).mean(), p.array(axh[k].get_ylim()).mean() , warningtext[k], \
                     horizontalalignment="center", verticalalignment="center", rotation=45, **(rc.warning+rc.alltext) )

def plotGeweke ( BayesInferenceObject, parameter=0, ax=None, warn=True ):
    """Geweke plot of moving average of samples

    :Parameters:
        *BayesInferenceObject* :
            a BayesInference object that contains all the
            infromation about the sampling process
        *parameter* :
            index of the model parameter of interest
        *ax* :
            the pylab.axes object where the plot should go
        *warn* :
            should a warning message be displayed if non stationarity
            of the samples is observed?
    """

    if BayesInferenceObject.mcestimates is None:
        raise ValueError, "Geweke MCMC convergence diagnostic requires monte carlo samples. Try to call the sample() method of your inference object."

    stationary,z,bad = BayesInferenceObject.geweke ( parameter )

    if ax is None:
        ax = prepare_axes (p.axes())

    x = N.arange (z[:,0].shape[0])+1
    for k in xrange ( z.shape[-1] ):
        ax.plot ( x, z[:,k], 'o-' )
    ax.plot ( ax.get_xlim(), [-2]*2, 'k:' )
    ax.plot ( ax.get_xlim(), [ 2]*2, 'k:' )

    if warn and not stationary:
        nsegments = z.shape[0]
        ax.text(0.5*nsegments,0,"chains did not converge", rotation=45,
                verticalalignment="center", horizontalalignment="center", **(rc.warning+rc.alltext) )

    ax.set_yticks ( N.array( (-3,-2,-1,0,1,2,3) ) )
    ax.set_xticks ( x )
    ax.set_ylim ( -3,3 )
    ax.set_xlim ( 0.5, z[:,k].shape[0]+.5 )

    ax.set_xlabel ( "chain segment", **(rc.label+rc.alltext) )
    ax.set_ylabel ( "z-score", **(rc.label+rc.alltext) )

    return ax

def plotChains ( BayesInferenceObject, parameter=0, ax=None, raw=False, warn=True ):
    """Simply plot all chains for a single parameter

    :Parameters:
        *parameter* :
            index of the model parameter to plot
        *raw* :
            plot raw samples instead of thinned samples after burnin
        *ax* :
            axes in which to print
        *warn* :
            if True, warnings are written into the plot
    """

    if BayesInferenceObject.mcestimates is None:
        raise ValueError, "Plotting MCMC chains requires monte carlo samples. Try to call the sample() method of your inference object."

    # Do we have an appropriate axis?
    if ax==None:
        ax = prepare_axes ( p.axes() )

    # Plot the chains
    for c in xrange(BayesInferenceObject.nchains):
        samples = BayesInferenceObject.getsamples ( c, raw=raw )
        ax.plot ( samples[:,parameter] )

    # Learn something about the axes
    xtics = N.array(ax.get_xticks())
    x0    = xtics.min()
    xr    = xtics.max()-xtics.min()
    ytics = ax.get_yticks()
    y0    = ytics.min()
    yr    = N.array(ytics.max()-ytics.min())

    if BayesInferenceObject.nchains>2:
        ax.text(x0+0.6*xr,y0+0.95*yr,r"$\hat{R} = %.4f$" % (BayesInferenceObject.Rhat ( parameter ) ) )

    if warn and BayesInferenceObject.Rhat(parameter)>1.1:
        ax.text(x0+0.5*xr,y0+0.5*yr,"Chains do not seem to sample\nfrom the same distribution!",
                horizontalalignment="center",verticalalignment="center",rotation=45,**(rc.warning+rc.alltext))

    ax.set_xlabel ( "sample #", **(rc.label+rc.alltext) )
    parname = BayesInferenceObject.parnames[parameter]
    if parname in ["alpha","beta","gamma","lambda"]:
        parname = r"$\%s$" % (parname,)
    ax.set_ylabel ( parname, **(rc.label+rc.alltext) )

    return ax

def plotParameterDist ( InferenceObject, parameter=0, ax=None ):
    """Plot the distribution of parameters

    :Parameters:
        *InferenceObject* :
            either a BootstrapInference object or a BayesInference object
            containing the samples of the parameter distribtution
        *parameter* :
            index of the model parameter of interest
        *ax* :
            pylab.axes object where the plot should go
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distribution of parameters requires monte carlo samples. Try to call the sample() method of your inference object."

    if ax is None:
        ax = prepare_axes ( p.axes() )

    samples = InferenceObject.mcestimates[:,parameter]
    h,b,ptch = ax.hist ( samples, bins=20, normed=True, histtype="step", lw=2 )

    if InferenceObject.__repr__().split()[1] in ["BayesInference","ASIRInference"]:
        priorstr = InferenceObject.model["priors"]
        if not priorstr is None:
            priorstr = priorstr[parameter]
            m = re.search (
                r"(\w+)\((-?\d*\.?\d*[eE]?-?\d*),(-?\d*\.?\d*[eE]?-?\d*)\)",
                priorstr )
            if not m is None:
                dist,prm1,prm2 = m.groups()
                prm1,prm2 = float(prm1),float(prm2)
                x = N.mgrid[b.min():b.max():100j]

                if dist.lower () == "gauss":
                    ax.plot(x,stats.norm.pdf(x,prm1,prm2))
                elif dist.lower () == "beta":
                    ax.plot(x,stats.beta.pdf(x,prm1,prm2))
                elif dist.lower () == "gamma":
                    ax.plot(x,stats.gamma.pdf(x,prm1,scale=prm2))
                elif dist.lower () == "ngamma":
                    ax.plot(x,stats.gamma.pdf(-x,prm1,scale=prm2))
                elif dist.lower () == "uniform":
                    ax.plot(x,stats.uniform.pdf(x,prm1,prm2))

    # Highlight estimate and credibility intervals
    prm = InferenceObject.estimate[parameter]
    c25,c975 = p.prctile ( samples, (2.5,97.5) )
    ym = ax.get_ylim()
    ax.plot ( [c25]*2,ym,'b:', [c975]*2,ym,'b:' )
    ax.plot ( [prm]*2,ym,'b' )

    prname = InferenceObject.parnames[parameter]
    if prname in ["alpha","beta","gamma","lambda"]:
        prname = r"\%s" % (prname,)
    message = r"$\hat{%s}"%(prname,)
    message += r"$=%.3f, CI(95)=(%.3f,%.3f)" % ( prm,c25,c975 )
    ax.set_title ( message, **(rc.text+rc.alltext) )
    ax.set_xlabel ( InferenceObject.parnames[parameter], **(rc.label+rc.alltext) )
    ax.set_ylabel ( "density estimate", **(rc.label+rc.alltext) )

def plotThresholdDist ( InferenceObject, cut=0, ax=None ):
    """Plot the distribution of thresholds

    :Parameters:
        *InferenceObjecxt* :
            a BootstrapInference or BayesInference object containing the desired
            data
        *cut* :
            index (!) of the desired cut
        *ax* :
            axes object to place the plot in.
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distributions of thresholds requires monte carlo samples. Try to call the sample() method of your inference object."

    if ax is None:
        ax = prepare_axes ( p.axes() )

    # Plot histogram
    mcthres = InferenceObject.mcthres[:,cut]
    h,b,ptch = ax.hist ( mcthres, bins=20, normed=True, histtype="step", lw=2 )

    # Highlight estimate and credibility intervals
    thres = InferenceObject.getThres ( InferenceObject.cuts[cut] )
    c25,c975 = InferenceObject.getCI ( cut, (0.025,0.975) )

    yl = ax.get_ylim ()
    ax.plot( [c25]*2,yl,'b:', [c975]*2,yl,'b:' )
    ax.plot ( [thres]*2, yl, 'b' )
    ax.set_title ( r"F$^{-1}$(%.2f)=%.3f, CI(95%%)=(%.3f,%.3f)" % (InferenceObject.cuts[cut], thres, c25, c975 ),
            horizontalalignment="center", verticalalignment="bottom", **(rc.text+rc.alltext) )

    ax.set_xlabel ( r"F$^{-1}$(%.2f)" % ( InferenceObject.cuts[cut], ), **(rc.label+rc.alltext) )
    ax.set_ylabel ( "density estimate", **(rc.label+rc.alltext) )

    ax.set_ylim ( yl )

    return ax

def ThresholdPlot ( InferenceObject ):
    """Show distributions and estimates for all thresholds

    :Parameters:
        *InferenceObject*
            a BootstrapInference or BayesInference object containing the
            desired data
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distributions of thresholds requires monte carlo samples. Try to call the sample() method of your inference object."

    nthres = len(InferenceObject.cuts)
    axw = 1./nthres
    fig = p.figure ( figsize=(3*nthres,3) )
    allax = axes_array_h ( fig, nthres, (.8/nthres-.08/(nthres-1),.7), (.1,.2), dist=.1 )

    for k,ax in enumerate ( allax ):
        # ax = p.subplot ( 1,nthres,k+1 )
        ax = plotThresholdDist ( InferenceObject, k, prepare_axes ( ax ) )

    return allax

def ParameterPlot ( InferenceObject ):
    """Show distributions and estimates for all parameters in the model

    :Parameters:
        *InferenceObject* :
            a BootstrapInference or BayesInference object containing the
            desired data
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distributions of parameters requires monte carlo samples. Try to call the sample() method of your inference object."

    nparams = len(InferenceObject.parnames)
    axw = 1./nparams
    fig = p.figure (figsize=(3*nparams,3))
    allax = axes_array_h ( fig, nparams, (.8/nparams-.08/(nparams-1),.65), (.1,.2), dist=.1 )

    for k,ax in enumerate ( allax ):
        # ax = p.subplot ( 1, nparams, k+1 )
        plotParameterDist ( InferenceObject, k, ax )

    return allax

def ConvergenceMCMC ( BayesInferenceObject, parameter=0, warn=True ):
    """Diagram to check convergence of MCMC chains for a single parameter

    :Parameters:
        *BayesInferenceObject* :
            a BayesInference object containing all information about
            the model and the posterior distribution
        *parameter* :
            model parameter of interest. So far, no model derived parameters such as
            thresholds are supported
        *warn* :
            should warnings be displayed if the samples look suspicious?
    """

    if BayesInferenceObject.mcestimates is None:
        raise ValueError, "MCMC convergence diagnostics require monte carlo samples. Try to call the sample() method of your inference object."

    fig = p.figure ( figsize=[9,3] )
    ax_chains,ax_geweke,ax_prm = axes_array_h ( fig, 3, (.2,.65),(.1,.2), dist=.1 )
    plotChains ( BayesInferenceObject, parameter, ax_chains, warn=warn )
    plotGeweke ( BayesInferenceObject, parameter, ax_geweke, warn=warn )
    plotParameterDist ( BayesInferenceObject, parameter, ax_prm )

def plotSensitivity ( BootstrapInferenceObject, ax=None ):
    """Visualize a sensitivity analysis to determine expanded bootstrap confidence intervals

    Sensitivity analysis is used for BootstrapInference objects to expand the confidence intervals
    in order to obtain more realistic coverage. This function calls the sensitivity_analysis() method
    of the BootstrapInferenceObject with default parameters. If other parameters are requested, the
    sensitivity_analysis() method should be called manually

    :Parameters:
        *BootstrapInferenceObject* :
            Inference object to be analyzed
        *ax* :
            pylab axes that should be used for plotting
    """
    if BootstrapInferenceObject.mcestimates is None:
        raise ValueError, "Sensitivity analysis requires monte carlo samples. Try to call the sample() method of your inference object."

    if ax==None:
        ax = p.axes()

    # Determine axes ranges
    prm1 = BootstrapInferenceObject.mcestimates[:,0]
    prm2 = BootstrapInferenceObject.mcestimates[:,1]
    ax.plot(prm1,prm2,'w.',markersize=1)
    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    ax.cla()

    # Plot the density estimate in the background
    x,y = N.mgrid[xmin:xmax:100j,ymin:ymax:100j]
    C = BootstrapInferenceObject.mcdensity(N.c_[N.ravel(x),N.ravel(y)].T)
    C.shape = 100,100
    ax.imshow( C.T,origin="lower",extent=(xmin,xmax,ymin,ymax), cmap=p.cm.gray_r )

    # Get the points and make sure, a sensitivity_analysis has indeed been run
    thres,slope,points = BootstrapInferenceObject.sensitivity_analysis(verbose=False)

    # plot the points
    ax.fill(points[:,0],points[:,1],fill=False,edgecolor="r",linewidth=2)
    ax.plot(prm1,prm2,"b.",markersize=2)
    ax.plot(points[:,0],points[:,1],'rd',markersize=5)
    ax.plot([BootstrapInferenceObject.estimate[0]],[BootstrapInferenceObject.estimate[1]],'ro',markersize=5)

    # plot marginal percentiles
    prm1lims = p.prctile ( BootstrapInferenceObject.mcestimates[:,0], (2.5,25,75,97.5) )
    prm2lims = p.prctile ( BootstrapInferenceObject.mcestimates[:,1], (2.5,25,75,97.5) )
    ax.plot( prm1lims, [ymin-0.05*(ymax-ymin)]*4, 'b-', [xmin-0.05*(xmax-xmin)]*4, prm2lims, 'b-' )
    ax.plot( prm1lims[1:3], [ymin-0.05*(ymax-ymin)]*2, 'b-', [xmin-0.05*(xmax-xmin)]*2, prm2lims[1:3], 'b-', linewidth=5 )

    # Draw axes
    drawaxes ( ax, ax.get_xticks(), "%g", ax.get_yticks(), "%g", BootstrapInferenceObject.parnames[0], BootstrapInferenceObject.parnames[1] )

def plotInfluential ( InferenceObject ):
    """Diagnostic plot for detecting influential observations

    Determining influential observations follows a different logic for bootstrap
    and for bayes inference. A block is labelled an influential observation if
    the fit for a dataset without that point is significantly different from the
    fit including that point. For BootstrapInference objects, this is quantified
    using a normed distance of the maximum likelihood fit including the block and
    withouth that block. This distance is normed in the following way: If the
    maximum likelihood fit for the reduced dataset remains inside the 95% confidence
    limits of the maximum likelihood fit for the full dataset, the influence
    value is below 1. Thus, influence values large than 1 indicate a problem with
    the data set. For BayesInference objects, the influence of a block is simply
    quantified as the Kullbach-Leibler divergence of the posterior for the full
    data set from the posterior for the reduced data set.

    :Parameters:
        *InferenceObject* :
            Data set for which the influential observations are to be plotted
    """
    maxinfl = N.argmax(InferenceObject.infl)
    ind = range ( InferenceObject.data.shape[0] )
    ind.pop(maxinfl)
    # influencedDataset = psignidata.BootstrapInference( InferenceObject.data[ind,:],
    #         sample=False, **(InferenceObject.model))
    # influencedDataset = psignidata.BayesInference ( InferenceObject.data[ind,:], **(InferenceObject.model) )
    est = interface.mapestimate ( InferenceObject.data[ind,:], start=InferenceObject.estimate, **(InferenceObject.model) )[0]
    x = N.mgrid[InferenceObject.data[:,0].min():InferenceObject.data[:,0].max():100j]
    influencedPMF = interface.diagnostics ( x, est,
            nafc = InferenceObject.model["nafc"],
            sigmoid = InferenceObject.model["sigmoid"],
            core = InferenceObject.model["core"] )

    p.figure ( figsize=(6,8) )
    # ax = p.axes ( (0.0,.5,.9,.5) )
    ax = prepare_axes ( p.subplot ( 2,1,1 ) )
    ax.set_ylabel ( r"$\Psi(x)$" )
    if InferenceObject.__repr__().split()[1] in ["BayesInference","ASIRInference"]:
        InferenceObject.drawposteriorexamples ( ax=ax )
    plotPMF ( InferenceObject, ax=ax, showaxes=True, showdesc=False, color="b", linewidth=2 )
    ax.plot ( [InferenceObject.data[maxinfl,0]], [InferenceObject.data[maxinfl,1].astype("d")/InferenceObject.data[maxinfl,2]],
            'rx', markersize=20, markeredgewidth=5 )
    # ax = plotPMF ( influencedDataset, ax=ax, showdesc=False, showaxes=True, color="r", markertype=([(0,0)],0), linewidth=2 )[-1]
    ax.plot ( x, influencedPMF, color="r", linewidth=2 )

    xl = list(ax.get_xlim ())

    # ax = p.axes ( (0.0, 0., .9, .5) )
    ax = p.subplot ( 2,1,2, sharex=ax )
    if InferenceObject.__repr__().split()[1] == "BootstrapInference":
        ax.plot ( [InferenceObject.data[:,0].min(),InferenceObject.data[:,0].max()], [1,1], 'k:' )
        yname = "Influence"
    else:
        yname = "D_KL( full || reduced )"
    ax.plot ( InferenceObject.data[:,0], InferenceObject.infl, 'bo' )
    ax.set_xlim(xl)
    drawaxes ( ax, ax.get_xticks(), "%g", ax.get_yticks(), "%g", r"stimulus intensity $x$", yname )

def plotMultiplePMFs ( *InferenceObjects, **kwargs ):
    """
    Plot multiple psychometric functions

    :Parameters:
        *InferenceObjects* :
            The Inference Objects that should be plotted. If the inference objects contain
            information about themselves, this information is used.
        *ax* :
            the axis object where the plot should go
        *xlabel* :
            text to be written on the y axis
        *ylabel* :
            text to be written on the x axis
        *ci* :
            boolean indicating whether credibility intervals should be drawn
            by default, this is False

    :Example:
    This example shows how to plot multiple psychometric functions

    >>> d0 = [[0, 28, 50], [2, 33, 50], [4, 38, 50], [6, 45, 50], [8, 45, 50], [10, 49, 50]]
    >>> d1 = [[0, 22, 50], [2, 34, 50], [4, 31, 50], [6, 42, 50], [8, 42, 50], [10, 46, 50]]
    >>> d2 = [[0, 26, 50], [2, 31, 50], [4, 38, 50], [6, 47, 50], [8, 49, 50], [10, 49, 50]]
    >>> constraints = ("","","Uniform(0,.1)")
    >>> B0 = BootstrapInference ( d0, priors=constraints,plotprm={"color": "r", "label": "Condition 0"} )
    >>> B1 = BootstrapInference ( d1, priors=constraints, plotprm={"color": "b","label": "Condition 1"} )
    >>> B2 = BootstrapInference ( d2, priors=constraints, plotprm={"color": "b","label": "Condition 2"} )
    >>> plotMultiplePMFs ( B0, B1, B2 )
    """
    ax = kwargs.setdefault ( "ax", None )
    if ax is None:
        ax = p.axes()
    pmflines = []
    pmflabels= []
    pmfdata  = []

    for pmf in InferenceObjects:
        l,d = plotPMF ( pmf, showaxes=False, showdesc=False, ax=ax )[:2]
        pmflines.append(l)
        pmfdata.append(d)
        pmflabels.append(pmf.label)
        if kwargs.setdefault ( "ci", False ):
            plotThres ( pmf, ax, color=pmf.color )

    ylabel_text = kwargs.setdefault("ylabel", None)
    if ylabel_text is None:
        if pmf.model["nafc"] < 2:
            ylabel_text = "P(Yes)"
        else:
            ylabel_text = "P(correct)"

    # Determine tics
    # p.setp(ax,frame_on=False,ylim=(-.05,1.05))
    xtics = p.getp(ax,'xticks')
    ytics = list(p.getp(ax,'yticks'))
    # Clean up ytics
    for k,yt in enumerate(ytics):
        if yt<0 or yt>1:
            ytics.pop(k)
    ytics = N.array(ytics)
    drawaxes ( ax, xtics, "%g", ytics, "%g", kwargs.setdefault("xlabel", "stimulus intensity"), ylabel_text )

    # Draw legend
    p.legend (pmflines,pmflabels)

    return pmflines,pmfdata

def plotJoint ( InferenceObject ):
    """Plot joint distribution of parameters
    
    For ASIRInference objects, also the fitted marginals are shown
    """
    nprm = InferenceObject.nparams
    parnames = InferenceObject.parnames

    fig = p.figure ()
    txt = []
    h,w = .9/nprm,.9/nprm
    for i in xrange ( nprm ):
        ax = prepare_axes ( fig.add_axes ( [.07+i*w,.97-(i+1)*h,.8*w,.8*h] ) )
        ax.xaxis.set_major_locator(MaxNLocator(5))

        if parnames[i] == "lambda":
            pr = r"\lambda"
        elif parnames[i] == "guess":
            pr = r"\gamma"
        else:
            pr = parnames[i]

        if InferenceObject.inference == "ASIR":
            gr  = InferenceObject.grids[i]
            mrg = InferenceObject.margins[i]
            ax.plot ( gr, mrg, 'b-.' )
            gn,gx = InferenceObject.getCI ( parnames[i], conf=(.01,.99) )
            x = p.mgrid[gn:gx:100j]
            ax.plot ( x, InferenceObject.posterior_pdf ( i, x ), 'b-', linewidth=2 )
            ax.plot ( x, InferenceObject.prior_pdf ( i, x ), 'k:' )
            txt.append ( r"$%s\sim%s$" % (pr, InferenceObject.posterior_approximations[i].strip("$")) )
        ax.set_xlabel ( r"$"+pr+r"$" )

        th = InferenceObject.mcestimates[:,i]
        if not (th==0).all():
            hist,b = N.histogram ( th, normed=True )
            ax.bar ( b[:-1], hist, N.diff(b), color=[.8,.8,1], edgecolor=[.8,.8,1] )

    for  i in xrange ( nprm ):
        ti = InferenceObject.mcestimates[:,i]
        for j in xrange ( i+1, nprm ):
            ax = prepare_axes ( fig.add_axes ( [.07+j*w, .97-(i+1)*h, .8*w,.8*h] ) )
            tj = InferenceObject.mcestimates[:,j]
            ax.plot ( tj, ti, '.' )
            ax.xaxis.set_major_locator(MaxNLocator(5))
            a,b,r,pr,se = stats.linregress ( tj, ti )
            x = N.sort(tj)
            ax.plot ( x, a*x+b, 'y' )
            ax.text ( ax.get_xlim()[0], ax.get_ylim()[0],r"$r=%.2f, p=%g$" % ( r,pr ),
                    fontsize=10, horizontalalignment="left", verticalalignment="bottom" )

    if InferenceObject.inference == "ASIR":
        txt.append ( r"duplicates: %g" % (InferenceObject.duplicates,) )
    fig.text ( .1,.1,"\n".join(txt) )


gof = GoodnessOfFit

if __name__ == "__main__":
    import doctest
    doctest.testmod()
