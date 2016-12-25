#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" Unit Tests for raw swig wrapper """

import numpy, pylab
import unittest as ut
import swignifit.swignifit_raw as sfr
import swignifit.utility as sfu
import pypsignifit.psignipriors as pfp

class TestData(ut.TestCase):

    @staticmethod
    def generate_test_dataset():
        x = sfr.vector_double([0.,2.,4.,6., 8., 10.])
        k = sfr.vector_int([24, 32, 40,48, 50,48])
        n = sfr.vector_int(6*[50])
        return sfr.PsiData(x, n, k, 2)

    def test_data(self):
        data = TestData.generate_test_dataset()
        data.setNcorrect(sfr.vector_int([24, 32, 40,48, 50,48]))

        data.getIntensities()
        data.getNtrials()
        data.getNcorrect()
        data.getPcorrect()

        blocks = data.getNblocks()

        for i in range(blocks):
            data.getIntensity(i)
            data.getNtrials(i)
            data.getNcorrect(i)
            data.getPcorrect(i)
            data.getNoverK(i)

        data.getNalternatives()
        data.nonasymptotic()

class TestSigmoid(ut.TestCase):
    """ test that all sigmoids have been wrapped and can be executed """

    def all_methods(self, sigmoid):
        s = sigmoid()
        s.f(0.0)
        s.df(0.0)
        s.ddf(0.0)
        s.inv(0.1)
        s.clone()
        s2 = sigmoid(s)

    def test_cauchy(self):
        self.all_methods(sfr.PsiCauchy)

    def test_exponential(self):
        self.all_methods(sfr.PsiExponential)

    def test_gauss(self):
        self.all_methods(sfr.PsiGauss)

    def test_gumbell(self):
        self.all_methods(sfr.PsiGumbelL)

    def test_gumbelr(self):
        self.all_methods(sfr.PsiGumbelR)

    def test_logistic(self):
        self.all_methods(sfr.PsiLogistic)

    def test_exponential_exception(self):
        s = sfr.PsiExponential()
        self.assertRaises(ValueError, s.inv, 0)
        self.assertRaises(ValueError, s.inv, 1)

class TestCore(ut.TestCase):

    data = TestData.generate_test_dataset()

    def all_methods(self, core):
        c = core(TestCore.data, 1, 0.1)
        params = sfr.vector_double([1.0,1.0])
        c.g(0.0, params)
        c.dg(0.0,params,0)
        c.dg(0.0,params,1)
        c.ddg(0.0,params,0,0)
        c.ddg(0.0,params,0,1)
        c.ddg(0.0,params,1,0)
        c.ddg(0.0,params,1,1)
        c.inv(0.0,params)
        c.dinv(0.0,params,0)
        c.dinv(0.0,params,1)
        c.transform(2,1.0,1.0)
        c.clone()
        c2 = core(c)

    def test_ab_core(self):
        self.all_methods(sfr.abCore)

    def test_linear_core(self):
        self.all_methods(sfr.linearCore)

    def test_log_core(self):
        self.all_methods(sfr.logCore)

    def test_mw_core(self):
        # mwCore constructor is a bit different than the rest
        self.all_methods(sfr.mwCore)

    def test_poly_core(self):
        self.all_methods(sfr.polyCore)

    def test_weibull_core(self):
        self.all_methods(sfr.weibullCore)

    def test_exceptions(self):
        c = sfr.logCore(TestCore.data)
        params = sfr.vector_double([1.0,1.0])
        self.assertRaises(ValueError, c.g, -1.0, params)
        c = sfr.weibullCore(TestCore.data)
        self.assertRaises(ValueError, c.dg, -1.0, params, 0)
        self.assertRaises(ValueError, c.ddg, -1.0, params, 0, 1)

class TestPriors(ut.TestCase):

    def all_methods(self, prior):
        p = prior(1.5, 3)
        p.pdf(0.0)
        p.dpdf(0.0)
        p.rand()
        p.clone()
        p2 = prior(p)
        p.get_code()

    def test_beta_prior(self):
        self.all_methods(sfr.BetaPrior)

    def test_gamma_prior(self):
        self.all_methods(sfr.GammaPrior)

    def test_ngamma_prior(self):
        self.all_methods(sfr.nGammaPrior)

    def test_gauss_prior(self):
        self.all_methods(sfr.GaussPrior)

    def test_uniform_prior(self):
        self.all_methods(sfr.UniformPrior)

class TestPsychometric(ut.TestCase):

    @staticmethod
    def generate_test_model():
        # IMPORTANT: here we can use the fact that PsiPsychometic manages its
        # own memory, and we don't need to hang on to th sigmoid, core, and
        # prior.
        return sfr.PsiPsychometric(2, sfr.abCore(), sfr.PsiLogistic())

    def test_pschometric(self):
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        params = sfr.vector_double([0.5,0.5,0.01])

        pr = sfr.UniformPrior(0,1)
        psi.setPrior(0,pr)
        psi.setPrior(1,pr)
        psi.setPrior(2,pr)

        psi.evaluate(0.0,params)
        psi.negllikeli(params,data)
        psi.neglpost(params, data)
        psi.leastfavourable(params, data, 0.0)
        psi.deviance(params, data)
        psi.ddnegllikeli(params, data)
        psi.dnegllikeli(params, data)
        psi.getCore()
        psi.getSigmoid()

    def test_memory_management(self):
        core = sfr.abCore()
        sigmoid = sfr.PsiLogistic()
        psi = sfr.PsiPsychometric(2,core,sigmoid)

    def test_exceptions(self):
        psi = TestPsychometric.generate_test_model()
        # for 2AFC we have 3 paramters with indices [0,1,2]
        self.assertRaises(ValueError, psi.setPrior,3, sfr.UniformPrior(0,1))

class TestOptimizer(ut.TestCase):

    def test_optimize(self):
        model = TestPsychometric.generate_test_model()
        data = TestData.generate_test_dataset()
        opt = sfr.PsiOptimizer(model, data)
        opt.optimize(model, data, None)

class TestBootstrap(ut.TestCase):

    @staticmethod
    def generate_test_bootstrap_list():
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        cuts = sfr.vector_double([1, 0.5])
        return sfr.bootstrap(999, data, psi, cuts)

    @staticmethod
    def generate_test_jackknife_list():
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        return sfr.jackknifedata(data, psi)

    def test_bootstrap(self):
        TestBootstrap.generate_test_bootstrap_list()

    def test_jackknifedata(self):
        TestBootstrap.generate_test_jackknife_list()

class TestMCMC(ut.TestCase):

    def all_sampler_methods(self, sampler):
        sampler.draw()
        theta = sampler.getTheta()
        sampler.setTheta(theta)
        sampler.setStepSize(sfr.vector_double([0.1,0.2,0.3]))
        sampler.getDeviance()
        sampler.sample(25)
        sampler.getModel()
        sampler.getData()

    def test_metropolis_hastings(self):
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        sampler = sfr.MetropolisHastings(psi, data, sfr.GaussRandom())
        self.all_sampler_methods(sampler)
        new_theta = sfr.vector_double(sampler.getNparams())
        sampler.proposePoint(sfr.vector_double(sampler.getTheta()),
                              sfr.vector_double(sampler.getStepsize()),
                              sfr.GaussRandom(),
                              new_theta)

    def test_generic_metropolis(self):
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        sampler = sfr.GenericMetropolis(psi, data, sfr.GaussRandom())
        self.all_sampler_methods(sampler)
        new_theta = sfr.vector_double(sampler.getNparams())
        sampler.proposePoint(sfr.vector_double(sampler.getTheta()),
                              sfr.vector_double(sampler.getStepsize()),
                              sfr.GaussRandom(),
                              new_theta)
        mclist = sampler.sample(4)
        sampler.findOptimalStepwidth(mclist)

    def test_independence_mcmc(self):
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        # just test initialization
        sampler = sfr.DefaultMCMC(psi, data, sfr.GaussRandom())
        # default priors has length 4
        priors = pfp.default(data.getIntensities())
        for i,p in enumerate(priors):
            sampler.set_proposal(i, sfu.get_prior(p))

    def test_hybrid_mcmc(self):
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        sampler = sfr.HybridMCMC(psi, data, 5)
        self.all_sampler_methods(sampler)

    def test_not_enough_samples(self):
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        sampler = sfr.GenericMetropolis(psi, data, sfr.GaussRandom())
        pilot = sampler.sample(2)
        self.assertRaises(ValueError, sampler.findOptimalStepwidth, pilot)

class TestMCList(ut.TestCase):

    def test_psi_mclist(self):
        bs_list = TestBootstrap.generate_test_bootstrap_list()
        bs_list.getEst(0)
        bs_list.getEst(0,0)
        bs_list.getPercentile(0.95, 0)
        bs_list.getMean(0)
        bs_list.getdeviance(0)
        bs_list.getNsamples()
        bs_list.getNparams()
        bs_list.getDeviancePercentile(0.95)

        bs_list.setEst(0, sfr.vector_double([0.1,0.1,0.1]), 0.95)
        bs_list.setdeviance(0,0.95)

    def test_bootstrap_list(self):
        bs_list = TestBootstrap.generate_test_bootstrap_list()
        bs_list.getData(0)
        # segmentation fault?
        #bs_list.getThres(0.95, 0)
        bs_list.getThres_byPos(0,0)
        bs_list.getNblocks()
        bs_list.getCut(0)
        bs_list.getAcc_t(0)
        bs_list.getAcc_s(0)
        bs_list.getBias_t(0)
        bs_list.getBias_s(0)
        bs_list.getRpd(0)
        # should this not throw a BadIndexError
        bs_list.percRpd(0)
        bs_list.getRkd(0)
        # should this not thow a BadIndexError?
        bs_list.percRkd(0)

        bs_list.setBCa_t(0, 0.1, 0.1)
        bs_list.setBCa_s(0, 0.1, 0.1)
        bs_list.setData(0, sfr.vector_int([24, 32, 40,48, 50,48]))
        bs_list.setThres(0.5, 0, 0)
        bs_list.setRpd(0, 0.5)
        bs_list.setRkd(0, 0.5)

    def test_jackknifedata(self):
        jk_list = TestBootstrap.generate_test_jackknife_list()
        jk_list.getNblocks()
        jk_list.influential(0, sfr.vector_double([0.0, 0.0, 0.0]),
                sfr.vector_double([0.0, 0.0, 0.0]))
        jk_list.outlier(0)

class TestRNG(ut.TestCase):

    def all_methods(self, random):
        random.draw()
        random.clone()

    def test_gauss_random(self):
        self.all_methods(sfr.GaussRandom())
        self.all_methods(sfr.GaussRandom(mean=5, standarddeviation=0.1))

    def test_uniform_random(self):
        self.all_methods(sfr.UniformRandom())
        self.all_methods(sfr.UniformRandom(low=-1, up=2))

    def test_binomial_random(self):
        binomial = sfr.BinomialRandom(6, 0.25)
        self.all_methods(binomial)
        binomial.setprm(10, 0.9)
        self.all_methods(binomial)

    def test_set_seed(self):
        sfr.setSeed(1)

class TestLinalg(ut.TestCase):

    def test_matrix(self):
        rows = columns = 5
        matrix = sfr.Matrix(rows, columns)
        for (i,j) in ((i,j) for i in xrange(rows) for j in xrange(columns)):
            # here we use the typemap magic from cpointer.i
            sfr.doublep_assign(matrix(i,j), 1 if i ==j else 0)
        # print is a python keyword, hence it was wrapped as _print
        # do not test print! 
        #matrix._print()
        matrix.getnrows()
        matrix.getncols()
        matrix.cholesky_dec()
        # TODO some of these fail due to segementation faults
        matrix.lu_dec()
        matrix.qr_dec()
        matrix.inverse_qr()
        matrix.regularized_inverse(0.1)
        matrix.solve(sfr.vector_double([0.1]*5))
        matrix.inverse()
        matrix * sfr.vector_double([0.5] * 5)
        matrix.scale(0.1)
        matrix.symmetric()
        # TODO use a dependent matrix (all ones) to test that decompositions,
        # solve and inverse don't work

if __name__ == "__main__":
    #ut.main()
    suite = ut.TestLoader().loadTestsFromName("__main__")
    ut.TextTestRunner(verbosity=2).run(suite)
