#!/usr/bin/env python

import pylab as p
import swignifit.swignifit_raw as sft

def prepare_axes ( ax, haveon=("bottom","left") ):
    """set spines correctly"""
    for loc,spine in ax.spines.iteritems():
        if loc in haveon:
            spine.set_position ( ("outward",10) )
        else:
            spine.set_color ( "none" )
    if "bottom" in haveon:
        ax.xaxis.set_ticks_position ( "bottom" )
    else:
        ax.xaxis.set_ticks_position ( "none" )
        ax.xaxis.set_ticklabels ( "" )
    if "left" in haveon:
        ax.yaxis.set_ticks_position ( "left" )
    else:
        ax.yaxis.set_ticks_position ( "none" )
        ax.yaxis.set_ticklabels ( "" )
    return ax

def axes_array_h ( fig, naxes, axsize, lowerleft=(0.1,0.1), dist=0.05, showally=False, nox=False ):
    """draw a horizontal array of axes"""
    xsize,ysize = axsize
    xdist,ydist = lowerleft
    step = xsize+dist
    
    if nox:
        axs = [prepare_axes ( fig.add_axes ( [xdist,ydist,xsize,ysize] ), haveon="left" )]
    else:
        axs = [prepare_axes ( fig.add_axes ( [xdist,ydist,xsize,ysize] ) )]
    for n in xrange(1,naxes):
        xdist += step
        if nox:
            axs.append ( prepare_axes ( fig.add_axes ( [xdist,ydist,xsize,ysize] ), haveon=[] ) )
        else:
            axs.append ( prepare_axes ( fig.add_axes ( [xdist,ydist,xsize,ysize] ), haveon=("bottom",) ) )

    return axs


fig = p.figure ( figsize=(6,6) )

ax_ab,ax_weibull = axes_array_h ( fig, 2, [.2,.2], [.5,.65], nox=True )
ax_logist,ax_logistab,ax_logistweib = axes_array_h ( fig, 3, [.2,.2], [.25,.37], nox=True )
ax_lgumbel,ax_gumbelab,ax_gumbelweib = axes_array_h ( fig, 3, [.2,.2], [.25,.1], nox=False )

fontspecs = {"horizontalalignment": "center", "verticalalignment": "center"}

fig.text ( .05, .47, r"logistic", rotation=90, **fontspecs)
fig.text ( .1, .47, r"$f(x)=\frac{1}{1+\exp(-x)}$", rotation=90, **fontspecs )
fig.text ( .05, .2, r"exponential", rotation=90, **fontspecs)
fig.text ( .1, .2, r"$f(x)=1-\exp(-x)$", rotation=90, **fontspecs )

fig.text ( .6, .95, r"ab", **fontspecs )
fig.text ( .6, .9, r"$g(x,a,b)=\frac{x-a}{b}$", **fontspecs )
fig.text ( .85, .95, r"poly", **fontspecs )
fig.text ( .85, .9, r"$g(x,a,b)=(\frac{x}{a})^b$", **fontspecs )

axes_ab = [ax_ab,ax_logist,ax_logistab,ax_lgumbel,ax_gumbelab]
axes_wei = [ax_weibull,ax_logistweib,ax_gumbelweib]

abcore   = sft.abCore ( )
weicore  = sft.polyCore ( sft.PsiData ( [1,2,3],[2,2,2],[2,2,2],2) )
logistic = sft.PsiLogistic ()
gumbel   = sft.PsiExponential ()

x = p.mgrid[-7:7:1000j]

ax_logist.plot  ( x, [logistic.f(xx) for xx in x] , color="k")
ax_lgumbel.plot ( x, [gumbel.f(xx)   for xx in x] , color="k")

prm = [(1,2,.02),(2,2,.02),(1,1,.02),(2,1,.02)]
col = ['b','r','g','c']

for pr,c in zip(prm,col):
    dab = [abcore.g  ( xx, pr ) for xx in x]
    ax_ab.plot       ( x, dab, color=c )
    ax_logistab.plot ( x, [logistic.f(dd) for dd in dab], color=c )
    ax_gumbelab.plot ( x, [gumbel.f(dd)   for dd in dab], color=c )

p.setp ( axes_ab, xlim = (x.min(),x.max()) )

x = p.mgrid[0:7:1000j]

for pr,c in zip(prm,col):
    dwei = [weicore.g ( xx, pr ) for xx in x]
    ax_weibull.plot  ( x, dwei, color=c )
    ax_logistweib.plot ( x, [logistic.f(dd) for dd in dwei], color=c )
    ax_gumbelweib.plot ( x, [gumbel.f(dd)   for dd in dwei], color=c )

p.setp ( axes_wei, xlim = (x.min(),x.max()) )
p.setp ( [ax_ab,ax_weibull], ylim=(-6,6) )
p.setp ( [ax_logist,ax_lgumbel,ax_logistab,ax_logistweib,ax_gumbelab,ax_gumbelweib], ylim=(0,1) )

k = 1
for pr,c in zip ( prm,col ):
    fig.text ( .2, .6+.05*k, "a=%g, b=%g" % pr[:2], color=c )
    k += 1

p.savefig ( "coreandsigmoid.png" )

p.show()
