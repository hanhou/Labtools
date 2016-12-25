#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__doc__ = """This file is intended to be run with

valgrind --tool=memcheck --leak-check=full /usr/bin/python memory_test.py

to detect memory leaks in the python interface."""

from pypsignifit.psigobservers import Observer
import swignifit.interface_methods as interface
import os
import pylab

O = Observer ( 4,2,.02 )
d = O.DoAnExperiment ( [0.,1.,2.,3,4.,5.,6.,7.,8.,9.], 20 )
priors = ("Gauss(0,100)","Gamma(1,4)","Beta(2,50)")

def sample ():
    boots = interface.bootstrap ( d, priors=priors, nsamples=1500-2*k )
    mape = interface.mapestimate ( d, priors=priors )
    mcmc = interface.mcmc ( d, start=(4,2,.02), priors=priors, nsamples = 1500-2*k )
    diag = interface.diagnostics ( d, (4,1,.02) )
    return float(os.popen ( "ps -C python -o rss" ).readlines()[1])/1024

n = []
n.append ( float(os.popen ( "ps -C python -o rss" ).readlines()[1])/1024 )
for k in xrange ( 200 ):
    print k
    n.append ( sample() )

pylab.plot(n)
pylab.ylabel("kB")
pylab.show()
