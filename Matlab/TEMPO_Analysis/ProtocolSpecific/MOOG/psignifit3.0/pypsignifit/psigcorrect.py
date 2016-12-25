#!/usr/bin/env python

import numpy as N
from scipy.special import polygamma,gamma

__doc__ = """This module implements a correction for overdispersion based on a beta-variance model

The procedure is described in more detail in

documents/correction_of_ci.pdf

In short, a beta distribution with the mode given by the fitted psychometric function,
and the number of trials scaled by a input independent factor nu is fitted to the data.
If a point estimate of the psychometric function has already been obtained, the only
free parameter of this model is nu which can be estimated using newton iterations.
The estimated nu can be used to correct the confidence intervals obtained from inference
that was based on a binomial distribution. Furhtermore, nu can be considered a measure of
binomiality of the residuals.
"""

__all__ = ["estimate_nu"]

def l ( nu, psi, k, n ):
    """log likelihood for psi

    :Parameters:
        *nu* :
            scalar value nu, for which the likelihood should be evaluated
        *psi* :
            m array that gives the psychometric function evaluated at stimulus
            levels x_i, i=1,...,m
        *k* :
            m array that gives the numbers of correct responses (in Yes/No: Yes responses)
            at stimulus levels x_i, i=1,...,m
        *n* :
            m array that gives the numbers of presented trials at stimulus
            levels x_i, i=1,...,m

    :Example:

    >>> psi = [ 0.52370051,  0.58115041,  0.70565915,  0.83343107,  0.89467234,  0.91364765,  0.91867512]
    >>> k   = [28, 29, 35, 41, 46, 46, 45]
    >>> n   = [50]*7
    >>> l ( 1, psi, k, n )
    13.752858759933943
    """
    psi = N.array(psi,'d')
    k = N.array(k,'d')
    n = N.array(n,'d')
    p = k/n
    return (N.log(gamma(nu*n+2))).sum()-(N.log(gamma(psi*nu*n+1))).sum()-(N.log(gamma((1-psi)*nu*n+1))).sum() \
            + (psi*nu*n*N.log(p)).sum() + ((1-psi)*nu*n*N.log(1-p)).sum()

def dl ( nu, psi, k, n ):
    """first derivative of the likelihood function with respect to nu

    :Parameters:
        *nu* :
            scalar value nu, for which the likelihood should be evaluated
        *psi* :
            m array that gives the psychometric function evaluated at stimulus
            levels x_i, i=1,...,m
        *k* :
            m array that gives the numbers of correct responses (in Yes/No: Yes responses)
            at stimulus levels x_i, i=1,...,m
        *n* :
            m array that gives the numbers of presented trials at stimulus
            levels x_i, i=1,...,m
    """
    p = k/n
    return (n*polygamma( 0, nu*n+2 )).sum() \
            - (psi*n*polygamma(0,nu*psi*n+1)).sum() \
            - ((1-psi)*n*polygamma(0,nu*(1-psi)*n+1)).sum() \
            + (psi*n*N.log(p)).sum() \
            + ((1-psi)*n*N.log(1-p)).sum()

def ddl ( nu, psi, k, n ):
    """second derivative of the likelihood function with respect to nu

    :Parameters:
        *nu* :
            scalar value nu, for which the likelihood should be evaluated
        *psi* :
            m array that gives the psychometric function evaluated at stimulus
            levels x_i, i=1,...,m
        *k* :
            m array that gives the numbers of correct responses (in Yes/No: Yes responses)
            at stimulus levels x_i, i=1,...,m
        *n* :
            m array that gives the numbers of presented trials at stimulus
            levels x_i, i=1,...,m
    """
    return (n**2*polygamma( 1, nu*n+2)).sum() \
            - (psi**2*n**2*polygamma(1,nu*psi*n+1)).sum() \
            - ((1-psi)**2*n**2*polygamma(1,nu*(1-psi)*n+1)).sum()

def estimate_nu ( InferenceObject ):
    """Perform a couple of newton iterations to estimate nu

    :Parameters:
        *InferenceObject* :
            An Inference object (either Bayes- or Bootstrap) for which the nu parameter should
            be estimated

    :Return:
        nu,nu_i,l_i
        *nu* :
            estimated nu parameter
        *nu_i* :
            sequence of nu values during optimization
        *l_i* :
            sequence of likelihoods associated with the nu values
    """
    psi = InferenceObject.evaluate ( InferenceObject.data[:,0] )
    k = InferenceObject.data[:,1].astype('d')
    n = InferenceObject.data[:,2].astype('d')
    k = N.where ( k==n, k-.01, k )
    k = N.where ( k==0, .01, k )

    nu = 1.
    nu_i = [nu]
    l_i = [l(nu,psi,k,n)]

    for i in xrange(10):
        nu -= dl(nu,psi,k,n)/ddl(nu,psi,k,n)
        if nu>1:
            nu=1
        elif nu<0:
            nu=0
        nu_i.append ( nu )
        l_i.append ( l(nu,psi,k,n) )

    return nu, nu_i, l_i

def main ( ):
    # An Example of usage
    from pypsignifit import BootstrapInference
    from pypsignifit.psigobservers import BetaBinomialObserver
    import pylab as p

    O = BetaBinomialObserver ( 4, 1, .02, 10 )

    d = N.array ( O.DoAnExperiment( [1,2,3,4,5,6,7] ) )
    B = BootstrapInference ( d, priors=("","","Uniform(0,.1)"))

    psi = B.evaluate ( B.data[:,0] )
    k = B.data[:,1].astype('d')
    n = B.data[:,2].astype('d')

    x = N.mgrid[0.001:0.999:100j]
    ll = []
    for xx in x:
        ll.append ( l(xx,psi,k,n) )
        print ll[-1]

    p.plot(x,ll,'b-')

    nu,nu_i,l_i = estimate_nu ( B )

    p.plot ( nu_i, l_i, 'ro' )
    p.show()

if __name__ == "__main__":
    # main()
    import doctest
    doctest.testmod()
