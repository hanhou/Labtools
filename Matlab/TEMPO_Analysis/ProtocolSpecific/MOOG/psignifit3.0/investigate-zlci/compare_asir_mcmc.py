#!/usr/bin/env python

import pypsignifit as psi
import pylab as pl
import numpy as np
import scipy.stats
import pypsignifit.psigobservers as Obs
from pypsignifit.psigniplot import plotJoint

sigmoid = "logistic"
core    = "mw0.1"
priors  = ("Gauss(0,100)","invGamma(.5,.5)","Beta(2,30)")
Beta = scipy.stats.beta
nafc = 2

blocks = [4,6,8,12]
blocksizes = [20,40,60]
wgen = [1,2,4]
nafc = [1,2]

def compare_on_dataset ( d, nafc=2, nruns=3 ):
    """Perform a comparison on a single dataset

    :Parameters:
        d
            dataset
        nafc
            number of alternatives
        nruns
            number of repetitions of the analysis

    TODO:
    This needs to be revised internally!
    """
    if nafc==2:
        priors = ("Gauss(0,100)","invGamma(.5,.5)","Beta(2,30)")
    else:
        priors = ("Gauss(0,100)","invGamma(.5,.5)","Beta(2,30)","Beta(2,30)")

    results = []
    for k in xrange ( nruns ):
        A = psi.ASIRInference ( d, nafc=nafc, priors=priors, sigmoid=sigmoid, core=core )

        Am025,Am975 = pl.prctile ( A.mcestimates[:,0], (2.5,97.5) )
        Aw025,Aw975 = pl.prctile ( A.mcestimates[:,1], (2.5,97.5) )
        Amm,Awm = A.mcestimates[:,:2].mean(0)
        Ams,Aws = A.mcestimates[:,:2].std(0)

        M = psi.BayesInference ( d, nafc=nafc, priors=priors, sigmoid=sigmoid, core=core, maxnsamples=1000000 )
        M.sample( start=M.farstart )
        M.sample( start=M.farstart )

        Mm025,Mm975 = pl.prctile ( M.mcestimates[:,0], (2.5,97.5) )
        Mw025,Mw975 = pl.prctile ( M.mcestimates[:,1], (2.5,97.5) )
        Mmm,Mwm = M.mcestimates[:,:2].mean(0)
        Mms,Mws = M.mcestimates[:,:2].std(0)
        MR = M.Rhat ()
        print MR

        result = "%g %g %g %g "  % (Am025,Am975,Amm,Ams)
        result += "%g %g %g %g " % (Aw025,Aw975,Awm,Aws)
        result += "%g %g %g %g " % (Mm025,Mm975,Mmm,Mms)
        result += "%g %g %g %g " % (Mw025,Mw975,Mwm,Mws)
        result += "%g %g %g " % (M.Rhat(0),M.Rhat(1),M.Rhat(2))
        result += "%g %g " % (A._ASIRInference__inference["resampling-entropy"],
                A._ASIRInference__inference["duplicates"])

        results.append ( result )

    return results

def sampling_scheme ( observer, nblocks ):
    if observer.model["nafc"] == 2:
        B = Beta ( 3.5, 1.5 )
        # B = Beta ( 1.5,.6 )
    elif observer.model["nafc"] == 1:
        B = Beta ( .5,.5 )
    Fx = B.ppf ( np.mgrid[.025:.975:nblocks*1j] )
    return observer.getlevels ( Fx )

def assemble_conditions ( blocks, blocksizes, wgen, nafc ):
    prm = []
    for b in blocks:
        for bs in blocksizes:
            for w in wgen:
                for n in nafc:
                    prm.append ( (b,bs,w,n) )
    return prm

def run_single_simulation ( prm, outf ):
    nblocks,blocksize,wgen,nafc = prm
    if nafc==2:
        O = Obs.Observer ( 4, wgen, .03, sigmoid=sigmoid, core=core, nafc=nafc )
    else:
        O = Obs.Observer ( 4, wgen, .03, .03, sigmoid=sigmoid, core=core, nafc=nafc )
    x = sampling_scheme ( O, nblocks )
    for di in xrange ( 3 ):
        d = O.DoAnExperiment ( x, blocksize )
        print d
        results = compare_on_dataset ( d, nafc )
        outf.write ( "%(nblocks)s %(ntrials)s %(nafc)s %(wgen)s " % {"nblocks": nblocks, "ntrials": blocksize, "nafc": nafc, "wgen": wgen} )
        outf.write ( " ".join ( results ) + "\n" )
    outf.flush()

if __name__ == "__main__":
    conditions = assemble_conditions ( blocks, blocksizes, wgen, nafc )
    outf = open ( "test", "w" )
    outf.write ( "nblocks ntrials nafc wgen Am025_1 Am975_1 Amm_1 Ams_1 Aw025_1 Aw975_1 Awm_1 Aws_1 Mm025_1 Mm975_1 Mmm_1 Mms_1 Mw025_1 Mw975_1 Mwm_1 Mws_1 Rm_1 Rw_1 Rl_1 H_1 dup_1 Am025_2 Am975_2 Amm_2 Ams_2 Aw025_2 Aw975_2 Awm_2 Aws_2 Mm025_2 Mm975_2 Mmm_2 Mms_2 Mw025_2 Mw975_2 Mwm_2 Mws_2 Rm_2 Rw_2 Rl_2 H_2 dup_2 Am025_3 Am975_3 Amm_3 Ams_3 Aw025_3 Aw975_3 Awm_3 Aws_3 Mm025_3 Mm975_3 Mmm_3 Mms_3 Mw025_3 Mw975_3 Mwm_3 Mws_3 Rm_3 Rw_3 Rl_3 H_3 dup_3\n" )

    for cond in conditions:
        run_single_simulation ( cond, outf )
    outf.close()
