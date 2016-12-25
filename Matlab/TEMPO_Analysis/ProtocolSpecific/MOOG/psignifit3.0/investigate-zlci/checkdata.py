#!/usr/bin/env python

import reader
import sys
import pypsignifit as pf
import pypsignifit.psignipriors as pfp
import pylab as pl
import numpy as np
import swignifit.swignifit_raw as sfr
import integrate as ig
# import pypsignifit.psigniplot as pp

d,s = reader.read_data ( sys.argv[1] )
d = np.array(d)
# stimulus_intensities = [0.0,2.0,4.0,6.0,8.0,10.0]
# number_of_correct = [34,32,40,48,50,48]
# number_of_trials  = [50]*len(stimulus_intensities)
# data = zip(stimulus_intensities,number_of_correct,number_of_trials)
# d = np.array ( data )

model = {'nafc':1, 'sigmoid':"logistic", 'core':'mw0.1'}
m = 4.0
w = 4.0
l = 0.05
g = 0.02

priors = ["Gauss(%f,%f)" % (m, m), "Gauss(%f,%f)" % (w, w*2), "Beta(2,50)", "Beta(1,50)"]
# priors = (pfp.default_mid(d[:,0])[0],"Gamma(2,4)",pfp.default_lapse(),pfp.default_lapse())
# priors = ("Gauss(4,.1)","Gauss(4,.1)","Beta(2,50)","Beta(1,50)")
print priors
x,fx,priors = ig.integration_grid ( d )
post_approx = ig.fit_posterior(fx,x)

print priors
print post_approx

# for i,prm in enumerate ( ["m","w","lm","gm"] ):
#     pl.subplot(221+i)
#     xx = np.mgrid[x[i].min():x[i].max():100j]
#     f = ig.sfu.get_prior ( post_approx[i] )
#     lf = np.array ([ f.pdf ( y ) for y in xx ])
#     pl.plot ( xx, lf, '-' )
#     pl.plot ( x[i], lf.max()*fx[i], 'o' )
#     pl.title ( prm )
# pl.show()


mcmc = pf.BayesInference(d, sigmoid=model['sigmoid'], core=model['core'],
        nafc=model['nafc'], priors=priors, verbose=True, sampler="independence",
        stepwidths=post_approx, automatic=False, sample=False )
mcmc.sample ( 8000 )

# print mcmc.approx
print post_approx

pf.ConvergenceMCMC(mcmc,0)
pf.ConvergenceMCMC(mcmc,1)
pf.ConvergenceMCMC(mcmc,2)
pf.ConvergenceMCMC(mcmc,3)
# pl.subplot(414)
# pl.plot ( mcmc.debug_samples[0][:,3] )
pl.show()

# pf.ConvergenceMCMC(mcmc,0)
# pf.ConvergenceMCMC(mcmc,1)
# pf.ConvergenceMCMC(mcmc,2)
# pf.show()
