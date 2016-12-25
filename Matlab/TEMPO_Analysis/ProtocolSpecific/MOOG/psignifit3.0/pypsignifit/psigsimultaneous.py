#!/usr/bin/env python

from scipy import stats
from scipy.special import gamma,digamma,polygamma
from numpy import log,mean,var,exp,array
from math import sqrt

__all__ = ["derive_informed_priors"]
__doc__ = """Routines to derive informed priors for quasi-simultaneous bayesian inference

These routines are described in more detail in the document located in

documents/simultaneous.pdf

In short, these routines take a number of fits that were performed with (improper) flat priors
and derive an list of informed priors. These informed priors can be used to perform the bayesian
inference a second time, this time with non flat prios, such that in the end the posterior
distributions of certain parameters are equal.
"""

def normpdf ( x, prm ):
    """normal pdf

    Parameters
    ----------
    x : array
        array on which the normal pdf should be evaluated
    prm : sequence
        pair of mean and variance for the gaussian

    Returns
    -------
    p : array
        an array of densities at the positions of x
    """
    return stats.norm.pdf ( x, prm[0], sqrt(prm[1]) )

def gammapdf ( x, prm ):
    """gamma pdf

    Parameters
    ----------
    x : array
        array on which the gamma pdf should be evaluated
    prm : sequence
        pair of k (shape) and theta (scale) parameters for the gamma distribution

    Returns
    -------
    p : array
        an array of densities at the positions of x
    """
    k,th = prm
    return x**(k-1) * exp ( - x/th )/(gamma(k)*th**k)

def betapdf ( x, prm ):
    """beta pdf

    Parameters
    ----------
    x : array
        array on which the gamma pdf should be evaluated
    prm : sequence
        pair of alpha (prior successes) and beta (prior misses) parameters for the beta distribution

    Returns
    -------
    p : array
        an array of densities at the positions of x
    """
    al,bt = prm
    return stats.beta.pdf(x, al, bt )

def fitgauss ( samples ):
    """fit a normal distribution using maximum likelihood

    Parameters
    ----------
    samples : array
        array of samples on which the distribution should be fitted

    Returns
    -------
    prm : sequence
        pair of mean and variance for the fitted gaussian
    l : double
        likelihood of the data at the fitted parameter values
    """
    m = mean ( samples )
    s = var ( samples )
    l = sum(log(normpdf(samples,(m,s))))
    return (m,s),l

def fitgamma ( samples ):
    """fit a gamma distribution using maximum likelihood

    Parameters
    ----------
    samples : array
        array of samples on which the distribution should be fitted

    Returns
    -------
    prm : sequence
        pair of k (shape) and theta (scale) parameters for the fitted gamma distribution
    l : double
        likelihood of the data at the fitted parameter values
    """
    s = log ( mean(samples) ) - mean ( log(samples) )
    k = 3 - s + sqrt ( (s-3)**2 + 24*s)
    k /= 12 * s
    for i in xrange ( 5 ):
        k -= ( log(k) - digamma ( k ) -s ) / ( 1./k - polygamma( 1, k ) )
    th = mean ( samples ) / k
    l = sum(log(gammapdf ( samples, (k,th) ) ) )
    return (k,th),l

def fitbeta ( samples ):
    """fit a beta distribution using maximum likelihood

    Parameters
    ----------
    samples : array
        array of samples on which the distribution should be fitted

    Returns
    -------
    prm : sequence
        pair of alpha (prior successes) and beta (prior misses) parameters for the beta distribution
    l : double
        likelihood of the data at the fitted parameter values
    """
    m = mean ( samples )
    s = var ( samples )
    al = m * ( m*(1-m)/s - 1 )
    bt = (1-m) * ( m*(1-m)/s - 1 )
    l = sum(log(betapdf(samples, (al, bt) ) ) )
    return (al,bt),l

def combinegauss ( prm ):
    """combine gaussian parameter estimates to give the product prior

    Parameters
    ----------
    prm : sequence
        a sequence of pairs of parameter estimates

    Returns:
    --------
    mubar,varbar : double
        mean an variance of the product prior
    """
    prm = array ( prm )
    varbar = 1./sum(1./prm[:,1])
    mubar = sum(prm[:,0]/prm[:,1]) * varbar
    return mubar,varbar

def combinegamma ( prm ):
    """combine gamma parameter estimates to give the product prior

    Parameters
    ----------
    prm : sequence
        a sequence of pairs of parameter estimates

    Returns:
    --------
    kbar,thetabar : double
        shape and scale parameter for the product prior gamma distribution
    """
    prm = array ( prm )
    kbar = 1 + sum (prm[:,0]-1)
    thbar = 1./sum(1./prm[:,1])
    return kbar,thbar

def combinebeta ( prm ):
    """combine gamma parameter estimates to give the product prior

    Parameters
    ----------
    prm : sequence
        a sequence of pairs of parameter estimates

    Returns:
    --------
    alphabar,betabar : double
        prior successes and prior misses parameters for the product prior beta distribution
    """
    prm = array ( prm )
    albar = 1 + sum ( prm[:,0]-1 )
    btbar = 1 + sum ( prm[:,1]-1 )
    return albar,btbar

def derive_informed_priors ( mcmcobjects, distribution="Gamma", parameter="w" ):
    """Combine mcmc samples to give informed priors for quasi-simultaneous inference

    For further information see the file 'simultaneous' in the documents folder.

    Parameters
    ----------
    mcmcobjects : sequence of BayesInference objects
        BayesInference objects containing samples obtained with flat priors that are used to determine
        the combined informed 'product prior'
    distribution : string
        name of the desired prior. Currently, only 'Gamma', 'Gauss', 'Beta' are supported.
    parameter : string
        name of the parameter to be constrained to be equal across conditions. Currently,
        'm', 'alpha', 'w', 'beta', 'lambda', and 'gamma' are supported

    Returns
    -------
    priors : sequence of strings
        list of priors to be used for performing the quasi simultaneous fit
    """
    # Check that the mcmc objects are comparable!!!
    params_assign = {'m': 0, 'alpha': 0, 'w': 1, 'beta': 1, 'lambda': 2, 'gamma': -1}
    if distribution=="Gamma":
        fit,combine = fitgamma,combinegamma
    elif distribution=="Gauss":
        fit,combine = fitgauss,combinegauss
    elif distribution=="Beta":
        fit,combine = fitbeta,combinebeta
    else:
        raise ArgumentError,"Invalid model to fit posterior distributions"

    fitted_parameters = []
    for mfit in mcmcobjects:
        prm,l = fit ( mfit.mcestimates[:,params_assign[parameter]] )
        fitted_parameters.append(prm)

    priors = []
    for j in xrange(len(fitted_parameters)):
        current = fitted_parameters.pop(0)
        prm = combine ( fitted_parameters )
        priors.append ( distribution+"(%g,%g)"%prm )
        fitted_parameters.append ( current )

    return priors
