#!/usr/bin/env python

from pylab import *

cores = [ lambda x,a,b: (x-a)/b,
        lambda x,a,b: 2*log(9)*(x-a)/b,
        lambda x,a,b: a*x+b,
        lambda x,a,b: where(x>0, a*log(x)+b, 0.),
        lambda x,a,b: 2./log(2)*a*b*(log(x)-log(a))+log(log(2)),
        lambda x,a,b: (x/a)**b ]
names = ["ab", "mw0.1 (for logistic)", "linear", "log", "weibull", "poly"]

x = mgrid[-3:3:100j]

for k,g in enumerate ( cores ):
    subplot(231+k)
    l = []
    l.append(plot ( x, g( x, 0, 1 ), label="g(x,0,1)" ))
    l.append(plot ( x, g( x, 2, 1 ), label="g(x,2,1)" ))
    l.append(plot ( x, g( x, 2, 2 ), label="g(x,2,2)" ))
    l.append(plot ( x, g( x, 0, 2 ), label="g(x,0,2)" ))
    title(names[k])
figlegend(l,("g(x,0,1)","g(x,2,1)","g(x,2,2)","g(x,0,2)"),loc="lower right")

savefig("cores.png")
