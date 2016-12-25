#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__docformat__ = "restructuredtext"

import sys
import subprocess

# This is the interface to psi++
import swignifit as interface

# This is "real" psignifit
from psignidata import *
from psigniplot import *

# Methods to set default priors
import psignipriors

# This is to enable display of graphics
from pylab import show

try:
    from __version__ import version
except ImportError:
    __version__ = 'Fatal: no version found!'

interface.set_seed( 0 )

def set_seed(value):
    interface.set_seed(value)

def dump_info():
    """
    Print some basic system info.
    
    NOTE: This will be extended to a more sophisticated scheme.
    """
    print("psignifit version: \t" + version)
    print("python version: \t" + sys.version)

def __test__ ( ):
    "If we call the file directly, we perform a test run"
    import numpy as N
    import pylab as p
    import sys

    if len(sys.argv) == 1 or sys.argv[1] == "bootstrap":
        bootstrap = True
    elif sys.argv[1] == "bayes":
        bootstrap = False

    x = [float(2*k) for k in xrange(6)]
    k = [34,32,40,48,50,48]
    n = [50]*6
    d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    d = N.array(zip(x,k,n))
    priors = ("flat","flat","Uniform(0,0.1)")

    if bootstrap:
        b = BootstrapInference ( d, sample=2000, priors=priors )
        GoodnessOfFit(b)
        ParameterPlot(b)
    else:
        priors = ("Gauss(0,4)","Gamma(1,3)","Beta(2,30)")
        mcmc = BayesInference ( d, sigmoid="cauchy", priors=priors )
        mcmc.sample(start=(6,4,.3))
        mcmc.sample(start=(1,1,.1))
        print "Posterior Intervals",mcmc.getCI()
        print "Model Evidence", mcmc.evidence
        print "Rhat (m):",mcmc.Rhat ()
        print "Nsamples:",mcmc.nsamples
        print "DIC:",mcmc.DIC
        print "pD:", mcmc.pD

        GoodnessOfFit(mcmc)
        for prm in xrange(3):
            ConvergenceMCMC ( mcmc, parameter=prm )
        print mcmc

    p.show()

if __name__ == "__main__":
    __test__()
