#!/usr/bin/env python

# This file illustrates the analysis of influential observations
# using constrained ML and bootstrap as well as bayesian inference
#
# The analysis is explained in more detail in the "Quick start to psignifit"
# that can be found at http://psignifit.sourceforge.net/TUTORIAL.html

from pypsignifit import *
from pylab import figure

# The data are form a 2afc task
nafc = 2

# Now we get the data
stimulus_intensities = [0.0,2.0,4.0,6.0,8.0,10.0]
number_of_correct = [34,32,40,48,50,48]
number_of_trials  = [50]*len(stimulus_intensities)
data = zip(stimulus_intensities,number_of_correct,number_of_trials)

# Constraints for Bootstrap
constraints = ( 'unconstrained', 'unconstrained', 'Uniform(0,0.1)' )

# Select priors
# These priors are clearly informative priors. This means they impose considerable
# prior information on the inference. This is useful in this case to illustrate
# bayesian analysis in general. However, in many cases, other, less informative priors
# might be of interest.
priors = ( 'Gauss(0,5)', 'Gamma(1,3)', 'Beta(2,30)' )

# Perform the actual inference
B    = BootstrapInference ( data, priors=constraints, nafc=nafc, sample=True )
mcmc = BayesInference ( data, priors=priors, nafc=nafc )

# plot BayesInfluence
plotInfluential ( B )
plotInfluential ( mcmc )

# for i in xrange (3 ):
#     ConvergenceMCMC ( mcmc, i )

# show results
show()
