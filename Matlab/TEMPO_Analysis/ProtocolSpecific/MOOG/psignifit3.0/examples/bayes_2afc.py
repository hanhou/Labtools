#!/usr/bin/env python

# This file illustrates the analysis of 2afc data using bayesian
# inference and markov chain monte carlo sampling.
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

# Select priors
# These priors are clearly informative priors. This means they impose considerable
# prior information on the inference. This is useful in this case to illustrate
# bayesian analysis in general. However, in many cases, other, less informative priors
# might be of interest.
priors = ( 'Gauss(0,5)', 'Gamma(1,3)', 'Beta(2,30)' )

# Perform the actual inference
mcmc = BayesInference ( data, priors=priors, nafc=nafc )

# Add some more chains for convergence diagnostic
mcmc.sample ( start = (0,1,0.01) )
mcmc.sample ( start = (6,11,0.3) )

# Generate convergence plots for all three paramters
for i in xrange ( 3 ):
    ConvergenceMCMC ( mcmc, i )

# Assess goodness of fit
GoodnessOfFit ( mcmc )

# See parameter plot
ParameterPlot ( mcmc )

# Show everything
show()
