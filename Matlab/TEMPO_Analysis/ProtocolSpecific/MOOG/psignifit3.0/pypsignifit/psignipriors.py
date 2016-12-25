#!/usr/bin/env python

__doc__ = """This module gives default priors for bayesian estimation of psychometric functions"""

import numpy as np
from scipy import stats, optimize

def default(x):
    """ All default priors, with default settings.

    :Parameters:
        *x* :
            array of stimulus intensities used in the experiment

    :Returns:
        tuple of 3 priors: m, w, lapse

    """
    return default_mid(x)[0], default_width(x)[0], default_lapse()[0]

def default_lapse ( observer="normal" ):
    """Default prior for the lapse rate

    :Parameters:
        *observer* :
            Typically, animals have much higher lapse rates than humans.
            The same is true for some clinical patient groups. Based on
            simulation results, we suggest different priors for animals or
            patients and for healthy adult humans. Thus, select either
            observer='normal' for healthy adult humans or observer='lapse'
            for an observer type with unusually high lapse rates.
    """
    if observer=="normal":
        return "Beta(2,20)",0.,0.5
    elif observer=="lapse":
        return "Beta(2.5,12)",0.,0.5
    else:
        raise Exception, "Unknown observer %s" % (observer,)

def default_width ( x, method="moments" ):
    """Default prior for the width of a psychometric function

    :Parameters:
        *x* :
            array of stimulus intensities used in the experiment
        *method* :
            method to determine the range spanned by the prior. This could either be
            'moments' which implies that the maximum width is mean+standard deviation of the
            gamma distribution and the minimum width is set to mean-standard deviation.
            Alternatively method could a sequence of quantiles which implies that the
            quantiles are matched with the maximum and minimum width.
    """
    xx = np.sort(x)
    wmin = np.min(np.diff(xx))
    wmax = 2*(xx[-1]-xx[0])
    if method=='moments':
        wr = wmin/wmax
        k  = ((1+wr)/(1-wr))**2
        th = wmax/(k+np.sqrt(k))
    elif isinstance ( method, (tuple,list) ):
        def error ( prm ):
            e = stats.gamma.cdf ( [wmin,wmax], prm[0], scale=prm[1] ) - np.array(method)
            e = np.sum(e**2)
            if np.isnan ( e ):
                return 10000.
            else:
                return e
        k,th = optimize.fmin ( error, [1.,4.] )

    return "Gamma(%g,%g)" % (k,th),wmin,wmax

def default_mid ( x, method="moments" ):
    """Default prior for the midpoint (threshold) of a psychometric function

    :Parameters:
        *x* :
            array of stimulus intensities used in the experiment
        *method* :
            method to determine the range spanned by the prior. This could either be
            'moments' which implies that the maximum stimulus intensity is mean+standard
            deviation of the normal distribution and the minimum stimulus intensity is
            set to mean-standard deviation. Alternatively method could a sequence of
            quantiles which implies that the quantiles are matched with the maximum and
            minimum stimulus intensities.
    """
    mmin = np.min(x)
    mmax = np.max(x)
    if method == "moments":
        mu = 0.5*(mmin+mmax)
        sg = 0.5*(mmax-mmin)
    elif isinstance ( method, (tuple,list) ):
        zmin,zmax = stats.norm.ppf ( method )
        sg = (mmax-mmin)/(zmax-zmin)
        mu = mmin - sg*zmin

    return "Gauss(%g,%g)" % (mu,sg), mmin, mmax

if __name__ == "__main__":
    import pylab as pl

    x = np.mgrid[0:10:100j]
    xx = np.array([1,2,3,5,6.,4.5,4])

    # g,k,th = default_width ( xx, [.1,.9] )
    # print g

    # # pl.plot ( xx, [0]*len(xx), 'o', x, (x/th)**(k-1)*np.exp(-x/th) )
    # pl.plot ( xx, [0]*len(xx), 'o', x, stats.gamma.pdf ( x, k, scale=th ), x, stats.gamma.cdf ( x, k, scale=th ) )

    # g,mu,sg = default_mid ( xx )
    # print g
    # pl.plot ( xx, [0]*len(xx), 'o', x, stats.norm.pdf ( x, mu, scale=sg ) )
    # pl.show()
