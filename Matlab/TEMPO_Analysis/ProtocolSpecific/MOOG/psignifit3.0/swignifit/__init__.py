#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

from interface_methods import bootstrap, mcmc, mapestimate, diagnostics, asir

def set_seed(value):
    swignifit_raw.setSeed(value)
