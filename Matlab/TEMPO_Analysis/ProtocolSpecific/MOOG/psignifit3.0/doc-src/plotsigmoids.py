#!/usr/bin/env python

from pylab import *
from scipy import stats


x = mgrid[-5:5:100j]

sigmoids = [ lambda x: 1./(1+exp(-x)),
        lambda x: stats.norm.cdf(x),
        lambda x: arctan ( x )/pi + 0.5,
        lambda x: 1-exp(-exp(x)),
        lambda x: exp(-exp(-x)),
        lambda x: where ( x>0, 1-exp(-x), 0 ) ]
names = ["logistic","gauss","cauchy","gumbel_l", "gumbel_r", "exponential"]

for k,f in enumerate ( sigmoids ):
    subplot(231+k)
    plot(x,f(x))
    title(names[k])
savefig ( "sigmoids.png" )
