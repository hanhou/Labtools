#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" Unit tests for interface methods provided by swignifit. """

import numpy as np
import numpy.testing as npt
import unittest as ut
import os
import swignifit.swignifit_raw as sfr
import swignifit.utility as sfu
import swignifit.interface_methods as interface

################################################################################
# original swignifit tests

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

class TestBootstrap(ut.TestCase):

    def test_basic(self):
        interface.bootstrap(data)

    def test_old_doctest(self):

        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('flat','flat','Uniform(0,0.1)')
        sfr.setSeed(1)
        samples,est,D,thres,thbias,thacc,slope,slbias,slacc,Rkd,Rpd,out,influ = interface.bootstrap(d,nsamples=2000,priors=priors)
        self.assertAlmostEqual( np.mean(est[:,0]), 2.7537742610139397, places=2)
        self.assertAlmostEqual( np.mean(est[:,1]), 1.4072288392075412, places=2)

    def test_start(self):
        interface.bootstrap(data, nsamples=25, start=[0.1, 0.2, 0.3])

    def test_nsamples(self):
        interface.bootstrap(data, nsamples=666)

    def test_nafc(self):
        interface.bootstrap(data, nafc=23)

    def test_sigmoid(self):
        interface.bootstrap(data, nsamples=25, sigmoid='gumbel_l')

    def test_core(self):
        interface.bootstrap(data, nsamples=25, core='linear')

    def test_prior(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        interface.bootstrap(data, nsamples=25, priors=priors)

    def test_cuts(self):
        interface.bootstrap(data, nsamples=25, cuts=[0.5,0.6,0.75])

    def test_parameteric(self):
        interface.bootstrap(data, nsamples=25, parametric=False)

class TestMCMC(ut.TestCase):

    def test_basic(self):
        interface.mcmc(data)

    def test_old_doctest(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')
        stepwidths = (1.,1.,0.01)
        sfr.setSeed(1)
        (estimates, deviance, posterior_predictive_data,
        posterior_predictive_deviances, posterior_predictive_Rpd,
        posterior_predictive_Rkd, logposterior_ratios, accept_rate) = interface.mcmc(d,nsamples=10000,priors=priors,stepwidths=stepwidths)
        self.assertAlmostEqual( np.mean(estimates[:,0]), 2.5463976926832483)
        self.assertAlmostEqual( np.mean(estimates[:,1]), 7.335732619111738)

    def test_start(self):
        interface.mcmc(data,nsamples=25, start=[0.1,0.2,0.3])

    def test_nsamples(self):
        interface.mcmc(data,nsamples=666)

    def test_nafc(self):
        interface.mcmc(data,nsamples=25, nafc=23)

    def test_sigmoid(self):
        interface.mcmc(data,nsamples=25, sigmoid='gumbel_r')

    def test_core(self):
        interface.mcmc(data, nsamples=25, core='ab')

    def test_prior(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        interface.mcmc(data, nsamples=25, priors=priors)

    def test_stepwidth(self):
        interface.mcmc(data, nsamples=25, stepwidths=[0.1, 0.2, 0.3])

    def test_alternate_samplers(self):
        # this used to fail for psipy, since it does not support alternative
        # samplers
        interface.mcmc(data, nsamples=25, sampler="MetropolisHastings")
        interface.mcmc(data, nsamples=25, sampler="GenericMetropolis")
        self.assertRaises(sfu.PsignifitException, interface.mcmc, data,
                sampler="DoesNotExist")

class TestMapestimate(ut.TestCase):

    def test_basic(self):
        interface.mapestimate(data)

    def test_old_doctest(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('flat','flat','Uniform(0,0.1)')
        estimate, fisher, thres, slope, deviance = interface.mapestimate ( d, priors=priors )
        npt.assert_almost_equal(np.array([ 2.7517523 ,  1.45728454,  0.01555464]), estimate)
        self.assertAlmostEqual(2.7517523000189286, thres[0])
        self.assertAlmostEqual(8.07133136423, deviance)

    def test_nafc(self):
        interface.mapestimate(data, nafc=23)

    def test_sigmoid(self):
        interface.mapestimate(data, sigmoid='gauss')

    def test_core(self):
        interface.mapestimate(data, core='mw0.2')

    def test_priors(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        interface.mapestimate(data, priors=priors)

    def test_cuts(self):
        interface.mapestimate(data, cuts=[0.5, 0.75, 0.85])

    def test_start(self):
        estimate, fisher, thres, slope, deviance = interface.mapestimate (data,
                start=[0.1, 0.2, 0.3])

class TestDiagnostics(ut.TestCase):

    prm = [2.75, 1.45, 0.015]

    def test_basic(self):
        interface.diagnostics(data, TestDiagnostics.prm)

    def test_old_doctest(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        prm = [2.75, 1.45, 0.015]
        pred,di,D,thres,slope,Rpd,Rkd = interface.diagnostics(d,prm)
        self.assertAlmostEqual(8.07484858608, D)
        self.assertAlmostEqual(1.68932796526, di[0])
        self.assertAlmostEqual(-0.19344675783032761, Rpd)

    def test_nafc(self):
        interface.diagnostics(data, TestDiagnostics.prm, nafc=23)

    def test_sigmoid(self):
        interface.diagnostics(data, TestDiagnostics.prm, sigmoid='logistic')

    def test_core(self):
        interface.diagnostics(data, TestDiagnostics.prm, core='linear')

    def test_cuts(self):
        interface.diagnostics(data, TestDiagnostics.prm, cuts=[0.5, 0.75, 0.85])

    def test_intensities_only(self):
        predicted = interface.diagnostics(x, TestDiagnostics.prm)

    def test_empty_data(self):
        # if an empty sequence is passed we only obtain the threshold
        result = interface.diagnostics([], TestDiagnostics.prm)

################################################################################
# old psipy tests, now also useful for swignifit


#################### 2afc #########################

def makedata ( neg=False ):
    if neg:
        return zip([1,2,3,4,5,6],[49,50,45,29,26,25],[50]*6)
    else:
        return zip([1,2,3,4,5,6],[25,26,29,45,50,49],[50]*6)

class TestPsipy_2afc ( ut.TestCase ):
    def test_bootstrap ( self ):
        """Call bootstrap"""
        priors = ("","","Uniform(0,.1)")
        interface.bootstrap ( makedata(), priors=priors )
    def test_mcmc ( self ):
        """Call mcmc"""
        priors = ("Gauss(0,100)","Gamma(1,4)","Beta(2,50)")
        interface.mcmc ( makedata(), start=(4,1,.02), priors=priors )
    def test_map ( self ):
        """Call mapestimator"""
        priors = ("","","Uniform(0,.1)")
        interface.mapestimate ( makedata(), priors=priors )
    def test_diagnostics ( self ):
        """Call diagnostics"""
        interface.diagnostics ( makedata(), (4,1,.02) )
    def test_sigmoids ( self ):
        """Try different sigmoids (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)")
        for s in ["logistic","cauchy","gauss","gumbel_r","gumbel_l","exponential"]:
            interface.mapestimate ( makedata(), priors=priors, sigmoid=s )
    def test_cores ( self ):
        """Try different cores (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)")
        for c in ["ab","mw0.1","poly","weibull","log","linear"]:
            interface.mapestimate ( makedata(), priors=priors, core=c )
    def test_weibull ( self ):
        """Try different parameterizations of the weibull (assumes working mapestimator)"""
        interface.mapestimate ( makedata(), priors = ("","","Uniform(0,.1)"), sigmoid="exponential", core="poly" )
        interface.mapestimate ( makedata(), priors = ("","","Uniform(0,.1)"), sigmoid="gumbel_r", core="weibull" )
        interface.mapestimate ( makedata(), priors = ("","","Uniform(0,.1)"), sigmoid="gumbel_r", core="log" )
    def test_priors ( self ):
        """Try different combinations of priors (assumes working mapestimator)"""
        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","Gamma(1,4)","Uniform(0,5)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    interface.mapestimate ( makedata(), priors = (mprior,wprior,lprior), core="mw0.1" )

        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","nGamma(1,4)","Uniform(-5,0)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    interface.mapestimate ( makedata(True), priors = (mprior,wprior,lprior), core="mw0.1" )

##################### yes/no ########################

def makedata_yn ( neg=False ):
    if neg:
        return zip ( [1,2,3,4,5,6],[40,46,27,25,14,1],[50]*6)
    else:
        return zip ( [1,2,3,4,5,6],[1,14,25,27,46,49],[50]*6)

class TestPsipy_yn ( ut.TestCase ):
    def test_bootstrap ( self ):
        """Call bootstrap"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        interface.bootstrap ( makedata_yn(), priors=priors, nafc=1 )
    def test_mcmc ( self ):
        """Call mcmc"""
        priors = ("Gauss(0,100)","Gamma(1,4)","Beta(2,50)","Beta(2,50)")
        interface.mcmc ( makedata_yn(), start=(4,1,.02,.02), priors=priors, nafc=1 )
    def test_map ( self ):
        """Call mapestimator"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        interface.mapestimate ( makedata_yn(), priors=priors, nafc=1 )
    def test_diagnostics ( self ):
        """Call diagnostics"""
        interface.diagnostics ( makedata_yn(), (4,1,.02,.02), nafc=1 )
    def test_sigmoids ( self ):
        """Try different sigmoids (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        for s in ["logistic","cauchy","gauss","gumbel_r","gumbel_l","exponential"]:
            interface.mapestimate ( makedata_yn(), priors=priors, sigmoid=s, nafc=1 )
    def test_cores ( self ):
        """Try different cores (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        for c in ["ab","mw0.1","poly","weibull","log","linear"]:
            interface.mapestimate ( makedata_yn(), priors=priors, core=c, nafc=1 )
    def test_weibull ( self ):
        """Try different parameterizations of the weibull (assumes working mapestimator)"""
        interface.mapestimate ( makedata_yn(), priors = ("","","Uniform(0,.1)","Uniform(0,.1)"), nafc=1, sigmoid="exponential", core="poly" )
        interface.mapestimate ( makedata_yn(), priors = ("","","Uniform(0,.1)","Uniform(0,.1)"), nafc=1, sigmoid="gumbel_r", core="weibull" )
        interface.mapestimate ( makedata_yn(), priors = ("","","Uniform(0,.1)","Uniform(0,.1)"), nafc=1, sigmoid="gumbel_r", core="log" )
    def test_priors ( self ):
        """Try different combinations of priors (assumes working mapestimator)"""
        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","Gamma(1,4)","Uniform(0,5)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    for gprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                        interface.mapestimate ( makedata_yn(), priors = (mprior,wprior,lprior,gprior), core="mw0.1", nafc=1 )

        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","nGamma(1,4)","Uniform(-5,0)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    for gprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                        interface.mapestimate ( makedata_yn(), priors = (mprior,wprior,lprior,gprior), core="mw0.1", nafc=1 )


if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromName("__main__")
    ut.TextTestRunner(verbosity=2).run(suite)

