from numpy import *
import pypsignifit.psignidata as D
import pypsignifit.psigniplot as P
from pylab import show, savefig, axes, plot, setp, fill, figure, rc

x = [0.,2.,4.,6.,8.,10.]
k = [34,32,40,50,48,50]
n = [50]*len(x)
d = zip(x,k,n)

x0 = -25

priors = ("Gauss(0,1e10)","Gamma(1,6)","Beta(2,50)")

figure ( figsize=(2.3,.5) )
ax = axes([0,0,1,1], frame_on=False,xticks=(),yticks=())

mcmc = D.BayesInference ( d, priors=priors)

print mcmc.burnin,mcmc.thin

est = mcmc.getsamples()
dev = mcmc.getmcdeviance()
dev -= dev.min()
dev /= dev.max()
dev = clip(.4+dev,0,1)

print est.shape,dev.shape,dev[mcmc.burnin::mcmc.thin].shape

x = mgrid[x0:10:100j]
for k in xrange ( 50 ):
    ind = random.randint ( dev.shape[0] )
    psi = mcmc.evaluate ( x, est[ind] )
    ax.plot(x,psi,color=[dev[ind]]*2+[1])

psi = mcmc.evaluate ( x )

ax.plot ( mcmc.data[:,0], mcmc.data[:,1]/mcmc.data[:,2], 'bo' )
ax.plot ( x, psi, 'b-', linewidth=2 )

# fill ( [x0-.2,1.7,1.7,x0-.2],[.53,.53,.64,.64], facecolor=[1,1,1], edgecolor=[1,1,1], alpha=.5, zorder=51 )

ax.text ( x0, .71, "psignifit", color="k", fontsize=28, horizontalalignment="left", verticalalignment="center", zorder=52 )

ax.text ( 6.4, .55, r"3.0", color="k", fontsize=17, horizontalalignment="center", verticalalignment="center" )
rc("text",usetex=True)
ax.text ( 9.6,.5, r"$\beta$", color="k", fontsize=14, horizontalalignment="center",verticalalignment="center" )

setp ( ax, xlim=(x0-.5,10+.5), ylim=(.4,1.1) )

# savefig ( "psignifit_logo.svg" )
savefig ( "psignifit_logo.png" )
savefig ( "psignifit_logo.eps" )
savefig ( "psignifit_logo.pdf" )
show()
