#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__docformat__ = "restructuredtext"

__doc__ = """This module implements the Raftery & Lewis gibbsit program in python."""

from numpy import *
import pylab as p
from scipy import stats

__all__ = ["gibbsit"]

def recode ( U, q ):
    """recode a series of the parameter of interest to a binary series

    :Parameters:
        U       chain of samples of the parameter of interest (1d-array!)
        q       quantile of interest

    :Output:
        Z       a series that is 1 for all U[t] less than the q-quantile
                of U and that is 0 otherwise
    """
    u = p.prctile ( U, 100*q )
    return (U<=u).astype("d")

def mctest ( Z ):
    """Compare 1st order and 2nd order Markov Chain models

    :Parameters:
        Z       a binary series of 0 and 1

    :Output:
        BIC,G2
        BIC     Bayesian information criterion for the comparison between
                1st order and 2nd order Markov Chain models. If BIC<=0,
                a 1st order model should be preferred.
        G2      the corresponding likelihood ratio statistic
    """
    # Determine Transition array
    T = zeros ( (2,2,2), 'd' )
    for i in xrange (2):
        for j in xrange (2):
            for k in xrange (2):
                T[i,j,k] = ( (Z[:-2]==i) & (Z[1:-1]==j) & (Z[2:]==k) ).sum()

    # Fitted model
    fitted = zeros ( (2,2,2), 'd' )
    for i in xrange (2):
        for j in xrange (2):
            for k in xrange (2):
                fitted[i,j,k] = T[i,j,:].sum() * T[:,j,k].sum() / T[:,j,:].sum()

    # Statistics
    G2 = 2 * sum ( T*log(T/fitted) )
    BIC = G2 - 2*log(Z.shape[0]-2)
    return BIC,G2

def indtest ( Z ):
    """Compare 1st order Markov Chain and independence chain

    :Parameters:
        Z       a binary series of 0 and 1

    :Output:
        BIC,G2
        BIC     Bayesian information criterion for the comparison between
                1st order Markov Chain and an independence chain, i.e. a
                chain in which all samples are independent of each others.
                If BIC<=0, the independence chain is to be preferred.
        G2      the corresponding likelihood ratio statistic
    """
    # Determine Transition matrix
    T = zeros ( (2,2), 'd' )
    for i in xrange (2):
        for j in xrange (2):
            T[i,j] = ( (Z[:-1]==i) & (Z[1:]==j) ).sum()

    # Fitted model
    fitted = zeros ( (2,2), 'd' )
    for i in xrange (2):
        for j in xrange (2):
            fitted[i,j] = T[i,:].sum() * T[j,:].sum() / (Z.shape[0]-1)

    # Statistics
    G2 = 2 * sum ( T*log(T/fitted) )
    BIC = G2 - 2*log(Z.shape[0]-1)
    return BIC,G2

def mcest ( Z ):
    """Estimate parameters of a 1st order Markov Chain

    :Parameters:
        Z       a binary series of 0 and 1

    :Output:
        alpha,beta
        alpha   Probability to move from state 0 to state 1
        beta    Probability to move from state 1 to state 0
    """
    # Determine Transition matrix
    T = zeros ( (2,2), 'd' )
    for i in xrange (2):
        for j in xrange (2):
            T[i,j] = ( (Z[:-1]==i) & (Z[1:]==j) ).sum()

    # Calculate transition probability
    alpha = T[0,1]/T[0,:].sum()
    beta  = T[1,0]/T[1,:].sum()
    return alpha,beta

def find_index ( Z, test=mctest, kstart=1 ):
    """Find an index based on a test function

    Find a thinning k from which a particular model is preferred.
    The function starts with k=kstart and performs model comparisons
    parameterized by the parameter test for increasing k. An IndexError is
    raised if k>200.

    :Parameters:
        Z       a binary series of 0 and 1
        test    a function, that returns the BIC for the test of interest
                a first argument. Useful functions for this task are
                mctest and indtest
        kstart  starting value for thinning to be used

    :Output:
        kstop   the first k value at which BIC<=0
    """
    k = kstart
    while True:
        BIC = test ( Z[::k] )[0]
        if BIC<=0:
            return k
        if k>200:
            raise IndexError
        k += 1

def gibbsit ( D=None, q=0.025, r=0.0125, s=0.95, eps=0.001 ):
    """The real gibbsit routine

    This function resembles the functionality of the gibbsit program by
    Raftery & Lewis (1995). It estimates parameters for a Markov Chain that
    should result in reasonable estimates of quantiles. These parameters
    are probabilistic in nature. Thus, they are to be taken as a good
    guess, not as exact values.

    :Parameters:
        D       a chain from a testrun of the sampling routine. If D is
                None, only the length that would be required for this chain
                is returned
        q       quantile of interest
        r       desired accuracy of the quantile of interest
        s       probability with which the estimated quantile is within
                q+/-r in the end
        eps     accuracy of estimation of the stationary density around the
                quantile

    :Output:
        params  a dictionary with the parameters for an MCMC run
    """
    # A correction parameter that is used multiple times
    phiterm = (stats.norm.ppf ( 0.5*(s+1) ) / r)**2

    class MCMCpar (dict):
        def __init__ (self, d):
            dict.__init__(self,d)
        burnin = property ( fget=lambda self: self.setdefault("M",0) )
        thin   = property ( fget=lambda self: self.setdefault("kthin",1) )
        Nsamples = property ( fget=lambda self: self.setdefault("N",self.setdefault("Nmin",1000) ) )
    params = MCMCpar({})

    if not D is None:
        Z = recode ( D, q )

        # Determine thinning
        kthin = find_index ( Z, mctest, 1 )       # Thinning required for 1st order markov chain
        kind  = find_index ( Z, indtest, kthin )  # Thinning required for independence chain

        # Estimate parameters of the 1st order markov chain ...
        alpha,beta = mcest ( Z[::kthin] )

        # Derive Quantities that have an explicit form
        ms = log ( (alpha+beta)*eps/max(alpha,beta) ) / log (1-alpha-beta)
        ns = (2-alpha-beta)*alpha*beta/(alpha+beta)**3 * phiterm
        params = MCMCpar({ "kthin": kthin, "kind": kind, "ms": ms, "M": int(ceil(ms*kthin)), "ns": ns, "N": int(ceil(ns*kthin)) })
    Nmin = q*(1-q) * phiterm
    params["Nmin"] = int(ceil(Nmin))

    return params

def main ( ):
    # Check whether we obtain similar results as the original
    D = fromfile ( "/home/ingo/tmp/mcmc_short.dat", sep=" " )
    D.shape = (-1,3)
    par = gibbsit(D[:,0])
    print par.thin,par.burnin,par.Nsamples

if __name__ == "__main__":
    main()
