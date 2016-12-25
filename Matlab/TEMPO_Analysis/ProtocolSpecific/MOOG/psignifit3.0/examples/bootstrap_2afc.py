#!/usr/bin/env python

# This file illustrates the analysis of 2afc data using constrained
# maximum likelihood and bootstrap analyzes.
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

# Constraints for the parameters
constraints = ( 'unconstrained', 'unconstrained', 'Uniform(0,0.1)' )

# Determine point estimate
# this uses the default values, i.e. a logistic sigmoid and the ab-core
# resulting in a parameterization of the psychometric function of the
# form
#
# psi ( x ) = gamma + (1-gamma-lambda) / ( 1 + exp ( - (x-alpha)/beta ) )
#
# With a parameter vector (alpha,beta,lambda) and gamma=1/2.
B = BootstrapInference ( data, priors=constraints, nafc=nafc )

# Now we perform bootstrap sampling to obtain confidence regions and goodness
# of fit statistics.
#
# Again the default values are used which is: 2000 bootstrap samples are generated
# by parametric bootstrap
B.sample ()

# We generate a summary of the goodness of fit statistics
GoodnessOfFit ( B )

# We plot information about the parameters and their distributions
ParameterPlot ( B )

# information about thresholds and their distributions
ThresholdPlot ( B )

# Now we print the confidence intervals of the 0.5-threshold before and after
# sensitivity analysis (where we perform the sensitivity analsis implicitely
# by calling plotSensitivity)
print "CI_0 =",B.getCI(1)
fig = figure()
ax = fig.add_axes ( [.1,.1,.8,.8] )
plotSensitivity ( B, ax )
print "CI_1 =",B.getCI(1)

# Show all figures
show()
