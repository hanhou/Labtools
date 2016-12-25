#!/usr/bin/env python

import numpy as np
from scipy.optimize import fmin,leastsq, brentq
import pypsignifit as pf
import pypsignifit.psignipriors as pfp
import swignifit.swignifit_raw as sfr
import swignifit.utility as sfu
import pypsignifit.psigobservers as pfo
from scipy import stats

# Error functions based on KL-Divergence?

def bounds ( mapest, pdata, ppmf, parameter="m" ):
    prm = mapest.copy()

    if mapest[2] < 0:
        prm[2] = 1e-5
    if len(mapest)>3:
        if mapest[3] < 0:
            prm[3] = 1e-5

    maxpost = -ppmf.neglpost ( prm, pdata )

    x = np.array ( [ pdata.getIntensity ( i ) for i in xrange ( pdata.getNblocks() ) ] )
    if parameter=="m":
        index = 0
        pmin = x.min ()
        pmax = x.max ()
    elif parameter=="w":
        index = 1
        x.sort ()
        pmin = 0.3*np.min ( np.diff(x) )
        pmax = 3*(x.max()-x.min())
    elif parameter=="lm":
        index = 2
        pmin = 0
        pmax = .5
    elif parameter=="gm":
        index = 3
        pmin = 0
        pmax = .5

    def error ( x, target=0.1*maxpost ):
        prm[index] = x
        lpost = -ppmf.neglpost ( prm, pdata )
        return lpost - target

    lmax = pmax
    lmin = pmin
    print "Maxpost =",maxpost
    for r in [.1,.3,.5,.8,.9]:
        trg = maxpost+np.log(r)
        if error ( mapest[index], trg )*error ( pmax, trg ) < 0:
            lmax = brentq ( error, mapest[index], pmax, args=(trg,) )
            print parameter,r,"max"
            break

    for r in [.1,.3,.5,.8,.9]:
        trg = maxpost+np.log(r)
        if error ( mapest[index], trg ) * error ( pmin, trg ) < 0:
            lmin = brentq ( error, pmin, mapest[index], args=(trg,) )
            print parameter,r,"min"
            break

    if lmin==pmin:
        print "WARNING: did not optimize lower bound for parameter",parameter
    if lmax==pmax:
        print "WARNING: did not optimize upper bound for parameter",parameter

    if lmin>lmax:
        lmin,lmax = lmax,lmin
    if parameter in ["lm","gm"]:
        if lmax<0.1:
            lmax = 0.1
        if lmin>0.02:
            lmin = 0.02

    return lmin,lmax

def dist2class ( dist ):
    name,prm = dist.split("(")
    prm = prm.strip(") ")
    p1,p2 = [float(p) for p in prm.split(",")]
    if name=="Gauss":
        f = stats.norm ( p1, p2 )
    elif name=="Gamma":
        f = stats.gamma ( p1, scale=p2 )
    elif name=="Beta":
        f = stats.beta ( p1, p2 )
    else:
        raise ValueError, "Unknown distribution: %s" % (name,)
    return f

def integration_grid ( data, run=1, dists=None ):
    data = np.array(data)
    mprior,mmin,mmax = pfp.default_mid ( data[:,0] )
    wprior,wmin,wmax = pfp.default_width ( data[:,0] )
    lprior,lmin,lmax = pfp.default_lapse ( )
    gprior,gmin,gmax = pfp.default_lapse ( )
    priors = (mprior,wprior,lprior,gprior)

    pdata,ppmf,pn = sfu.make_dataset_and_pmf ( data, 1, "logistic", "mw0.1", priors )
    mapest = pf.BootstrapInference ( data, priors, nafc=1 ).estimate

    if run==1:
        gridsize = 9
        mmin,mmax = bounds ( mapest, pdata, ppmf, "m" )
        wmin,wmax = bounds ( mapest, pdata, ppmf, "w" )
        lmin,lmax = bounds ( mapest, pdata, ppmf, "lm" )
        gmin,gmax = bounds ( mapest, pdata, ppmf, "gm" )
        grid = np.reshape ( np.mgrid[
            mmin:mmax:1j*gridsize,
            wmin:wmax:1j*gridsize,
            lmin:lmax:1j*gridsize,
            gmin:gmax:1j*gridsize
            ], (4,-1) )
    elif run==2:
        f_m = dist2class ( dists[0] )
        f_w = dist2class ( dists[1] )
        f_l = dist2class ( dists[2] )
        f_g = dist2class ( dists[3] )
        grid = np.mgrid[.025:.975:7j,.025:.975:7j,.025:.975:7j,.025:.975:7j]
        grid[0] = f_m.ppf ( grid[0] )
        grid[1] = f_w.ppf ( grid[1] )
        grid[2] = f_l.ppf ( grid[2] )
        grid[3] = f_g.ppf ( grid[3] )
        grid = np.reshape ( grid, (4,-1) )
        gridsize = 7

    post = np.reshape ( np.array ( map ( lambda prm: ppmf.neglpost ( prm, pdata ), grid.T ) ), [gridsize]*pn ) # negative log posterior
    post = np.exp ( -post ) # posterior
    grid = np.reshape ( grid, [pn]+[gridsize]*pn )

    d = [ np.diff ( grid[i], axis=i )[0].max() for i in xrange ( pn ) ]

    fx = [ marginalize ( post, d, i ) for i in xrange ( pn ) ]
    x = []
    for i in xrange ( pn ):
        s =  [ i ] + [0]*pn
        s[i+1] = slice(0,grid.shape[i+1])
        x.append ( grid[s] )
    return x,fx,priors

def marginalize ( post, d, i ):
    p = reduce ( lambda x,y: x*y, d )
    fx = np.zeros ( post.shape[i] )
    for j in xrange ( post.shape[i] ):
        s = [slice(0,post.shape[0]),slice(0,post.shape[1]),slice(0,post.shape[2]),slice(0,post.shape[3])]
        s[i] = j
        fx[j] = post[s].sum()
    # return fx * p / d[i]
    return fx

def error_gauss ( prm, fx, x ):
    Z,mu,sg = prm
    sg = sg*sg
    # return np.sum ( ( Z**2*np.exp ( -0.5*((x-mu)/sg)**2 ) - fx )**2 )
    # return Z*Z*np.exp ( -0.5*((x-mu)/sg**2 ) ) - fx
    return np.sum ( np.log ( Z**2*np.exp ( -0.5*((x-mu)/sg)**2 ) / fx )**2 )

def error_gamma ( prm, fx, x ):
    Z,k,th = prm
    k = k*k
    th = th*th
    return np.sum ( ( Z**2*x**(k-1)*np.exp(-x/th) - fx )**2 )
    # return np.sum ( np.log ( Z**2*x**(k-1)*np.exp(-x/th) / fx )**2 )

def error_beta ( prm, fx, x ):
    Z,al,bt = prm
    al = al**2
    bt = bt**2
    # return np.sum ( ( Z**2*x**(al-1)*(1-x)**(bt-1) - fx )**2 )
    return np.sum ( np.log( Z**2*x**(al-1)*(1-x)**(bt-1) / fx )**2 )

def fit_posterior ( fx, x ):
    post = []
    I = 10000
    N = 10000

    mu = x[0][np.argmax(fx[0])]
    fx[0] /= fx[0].max()
    mprm = fmin ( error_gauss, [1.,mu,1.5], args=(fx[0],x[0]), maxfun=N, maxiter=I )
    print mprm
    post.append ( "Gauss(%g,%g)" % ( mprm[1],mprm[2]**2 ) )

    fx[1] /= fx[1].max()
    wprm = fmin ( error_gamma, [1.,2,4], args=(fx[1],x[1]), maxfun=N, maxiter=I )
    post.append ( "Gamma(%g,%g)" % ( wprm[1]**2,wprm[2]**2 ) )

    fx[2] /= fx[2].max()
    lprm = fmin ( error_beta, [1.,2,20], args=(fx[2],x[2]), maxfun=N, maxiter=I )
    post.append ( "Beta(%g,%g)" % ( lprm[1]**2,lprm[2]**2 ) )

    if len(fx)>3:
        fx[3] /=  fx[3].max()
        gprm = fmin ( error_beta, [1.,2,20], args=(fx[3],x[3]), maxfun=N, maxiter=I )
        post.append ( "Beta(%g,%g)" % ( gprm[1]**2,gprm[2]**2 ) )

    return post

def sample_importance_resample ( post, pdata, ppmf, nresample=600, nsamples=6000 ):
    f = [ dist2class ( dist ) for dist in post ]
    print "Sample"
    presamples = []
    for dist in f:
        presamples.append ( dist.rvs ( nsamples ) )
    presamples = np.array(presamples)

    print "Calculate importance weights"
    w = []
    for prm in presamples.T:
        q = 1.
        for th,dist in zip ( prm, f ):
            q *= dist.pdf(th)
        p = np.exp ( -ppmf.neglpost ( prm, pdata ) )
        w.append ( p/q )
    w = np.array ( w )
    w /= w.sum()
    P = np.cumsum(w)
    u = np.sort ( np.random.rand ( nresample ) )

    print "Resample"
    samples = []
    ind = 0
    k = 0
    while k<nresample:
        n = 0
        while u[k] <= P[ind]:
            k += 1
            n += 1
            if k>=nresample:
                break
        samples += [presamples[:,ind]]*n
        ind += 1
    return np.array(samples)

if __name__ == "__main__":
    import pylab as pl
    O = pfo.Observer ( 5,3,.05,.05, core="mw0.1", sigmoid="logistic", nafc=1 )
    data = O.DoAnExperiment ( [1,2,3,4,5,6,7,8,12], 30 )
    # data = [[1, 2, 30], [2, 2, 30], [3, 2, 30], [4, 5, 30], [5, 16, 30], [6, 22, 30], [7, 26, 30], [8, 27, 30], [12, 29, 30]]
    # data = [[1, 1, 30], [2, 4, 30], [3, 3, 30], [4, 7, 30], [5, 11, 30], [6, 26, 30], [7, 27, 30], [8, 29, 30], [12, 30, 30]]
    # data = [[1, 1, 30], [2, 2, 30], [3, 4, 30], [4, 9, 30], [5, 19, 30], [6, 27, 30], [7, 25, 30], [8, 27, 30], [12, 30, 30]]  # First optimization does not converge
    # data = [[1, 1, 30], [2, 2, 30], [3, 5, 30], [4, 9, 30], [5, 20, 30], [6, 25, 30], [7, 29, 30], [8, 24, 30], [12, 28, 30]]
    # data = [[1, 3, 30], [2, 0, 30], [3, 0, 30], [4, 5, 30], [5, 17, 30], [6, 22, 30], [7, 24, 30], [8, 27, 30], [12, 29, 30]] # Bad initial fit ~> log fitting?
    # data = [[1, 1, 30], [2, 2, 30], [3, 1, 30], [4, 5, 30], [5, 13, 30], [6, 28, 30], [7, 26, 30], [8, 29, 30], [12, 27, 30]] # Takes many refinements to converge
    # data = [[1, 0, 30], [2, 2, 30], [3, 1, 30], [4, 9, 30], [5, 17, 30], [6, 24, 30], [7, 27, 30], [8, 28, 30], [12, 29, 30]] # Bad initial fit for w
    # data = [[1, 1, 30], [2, 2, 30], [3, 1, 30], [4, 3, 30], [5, 14, 30], [6, 26, 30], [7, 28, 30], [8, 28, 30], [12, 28, 30]] # Bad initial fit for w
    # data = [[1, 0, 30], [2, 2, 30], [3, 2, 30], [4, 9, 30], [5, 17, 30], [6, 26, 30], [7, 25, 30], [8, 28, 30], [12, 29, 30]] # Bad initial fit for w
    # data = [[1, 2, 30], [2, 4, 30], [3, 2, 30], [4, 9, 30], [5, 16, 30], [6, 17, 30], [7, 27, 30], [8, 27, 30], [12, 30, 30]] # Bad initial fit for w
    nrefine = 2
    print data

    x,fx,priors = integration_grid ( data )
    print "x1 =",x
    print "f1 =",fx
    post = fit_posterior(fx,x)
    for i in xrange ( nrefine ):
        x,fx,priors = integration_grid ( data, 2, post )
        print post
        print "x%d =" % (i+1,),x
        print "f%d =" % (i+1),fx
        post = fit_posterior (fx, x)
    f = [ sfu.get_prior ( p ) for p in post ]

    mapest = pf.BootstrapInference ( data, priors, core="mw0.1", nafc=1 ).estimate
    pdata,ppmf,pn = sfu.make_dataset_and_pmf ( data, 1, "logistic", "mw0.1", priors )
    samples = sample_importance_resample ( post, pdata, ppmf, nresample=600, nsamples=6000 )

    rng = [(3,7),(0,6),(0,.5),(0,.5)]

    hist_ax = [pl.axes ( [.15+.2*i,.75-.2*i,.15,.15] ) for i in xrange ( 4 )]
    labels = ["m","w","lm","gm"]

    for i,prm in enumerate ( labels ):
        ax = hist_ax[i]
        h,b = np.histogram ( samples[:,i], normed=True )
        ax.step ( np.convolve ( [.5,.5], b, 'valid' ), h, where='mid' )
        xx = np.mgrid[rng[i][0]:rng[i][1]:1000j]
        g = np.array ( [f[i].pdf(x_) for x_ in xx] )
        ax.plot ( xx, g, '-', linewidth=2 )
        mxind = np.argmax(fx[i])
        r = fx[i][mxind]/f[i].pdf(x[i][mxind])
        ax.plot ( x[i], fx[i]/r, 'o' )
        ax.plot ( [O.params[i]]*2,[0,f[i].pdf(O.params[i])], 'k' )
        ax.plot ( [mapest[i]]*2,[0,f[i].pdf(mapest[i])], 'r' )
        ax.set_title ( prm  )
        print fx[i]/r
        ax.set_ylim ( 0, 1.5 * np.max(fx[i]/r) )
        ax.set_xlim ( rng[i] )

    for i in xrange ( 4 ):
        for j in xrange ( i+1, 4 ):
            ax = pl.axes ( [.15+.2*j,.75-.2*i,.15,.15 ] )
            ax.plot ( samples[:,j], samples[:,i], '.' )

            m,b,r,pr,se = stats.linregress ( samples[:,j], samples[:,i] )
            x = np.sort ( samples[:,j] )
            ax.plot ( x, m*x+b )
            ax.set_title ( r"$r=%.2f\pm%.2f$" % (r,se) )

            ax.set_xlim ( rng[j] )
            ax.set_ylim ( rng[i] )

            pl.setp ( ax, xticklabels=() )
            if j-i > 1:
                pl.setp ( ax, yticklabels=() )

    ax = pl.axes ( [.1,.1,.3,.3] )
    data = np.array(data, dtype="d")
    x = np.mgrid[data[:,0].min():data[:,0].max():100j]
    for k in xrange ( 20 ):
        ax.plot ( x, [ppmf.evaluate ( xx, samples[k,:] ) for xx in x], color=[.7,.7,1.] )
    ax.plot ( x, [ppmf.evaluate ( xx, samples.mean(0) ) for xx in x], color='b', linewidth=2 )
    ax.plot ( data[:,0], data[:,1]/data[:,2], 'ko' )
    ax.set_xlabel ( "stimulus intensity" )
    ax.set_ylabel ( r"$\Psi(x|\theta)$" )

    print "Number of duplicates:",1-float(len(np.unique(samples[:,0])))/samples.shape[0]

    print mapest

    pl.show()

    print "data =",data
