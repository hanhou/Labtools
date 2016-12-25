#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

import unittest as ut
import numpy as np
import swignifit.swignifit_raw as sfr
import swignifit.utility as sfu
from interface_test import k, n, x, data


class TestUtility(ut.TestCase):

    def test_make_dataset(self):
        dataset = sfu.make_dataset(data, 1)
        self.assertTrue((np.array(x) == np.array(dataset.getIntensities())).all())
        self.assertTrue((np.array(k) == np.array(dataset.getNcorrect())).all())
        self.assertTrue((np.array(n) == np.array(dataset.getNtrials())).all())
        self.assertEqual(1, dataset.getNalternatives())

    def test_get_sigmoid(self):
        logistic = sfr.PsiLogistic()
        self.assertEqual("PsiLogistic", logistic.__class__.__name__)
        gauss = sfr.PsiGauss()
        self.assertEqual("PsiGauss", gauss.__class__.__name__)
        gumbel_l = sfr.PsiGumbelL()
        self.assertEqual("PsiGumbelL", gumbel_l.__class__.__name__)
        gumbel_r = sfr.PsiGumbelR()
        self.assertEqual("PsiGumbelR", gumbel_r.__class__.__name__)
        cauchy = sfr.PsiCauchy()
        self.assertEqual("PsiCauchy", cauchy.__class__.__name__)
        exponential = sfr.PsiExponential()
        self.assertEqual("PsiExponential", exponential.__class__.__name__)


    def test_get_core(self):
        sigmoid = sfr.PsiLogistic()
        dataset = sfu.make_dataset(data, 1)
        ab = sfu.get_core("ab", dataset, sigmoid)
        self.assertEqual("abCore", ab.__class__.__name__)
        mw = sfu.get_core("mw", dataset, sigmoid)
        self.assertEqual("mwCore", mw.__class__.__name__)
        self.assertEqual(0.1, mw.getAlpha())
        mw = sfu.get_core("mw0.2", dataset, sigmoid)
        self.assertEqual("mwCore", mw.__class__.__name__)
        self.assertEqual(0.2, mw.getAlpha())
        linear = sfu.get_core("linear", dataset, sigmoid)
        self.assertEqual("linearCore", linear.__class__.__name__)
        log = sfu.get_core("log", dataset, sigmoid)
        self.assertEqual("logCore", log.__class__.__name__)
        weibull = sfu.get_core("weibull", dataset, sigmoid)
        self.assertEqual("weibullCore", weibull.__class__.__name__)
        poly = sfu.get_core("poly", dataset, sigmoid)
        self.assertEqual("polyCore", poly.__class__.__name__)

    def test_get_prior(self):
        uniform = sfu.get_prior("Uniform(1,2)")
        self.assertEqual("UniformPrior", uniform.__class__.__name__)
        gauss = sfu.get_prior("Gauss(0,1)")
        self.assertEqual("GaussPrior", gauss.__class__.__name__)
        beta = sfu.get_prior("Beta(1.5, 3)")
        self.assertEqual("BetaPrior", beta.__class__.__name__)
        gamma = sfu.get_prior("Gamma(1.5, 3)")
        self.assertEqual("GammaPrior", gamma.__class__.__name__)
        ngamma = sfu.get_prior("nGamma(1.5,3)")
        self.assertEqual("nGammaPrior", ngamma.__class__.__name__)
        flat = sfu.get_prior("flat")
        self.assertEqual(None, flat)
        unconstrained = sfu.get_prior("unconstrained")
        self.assertEqual(None, unconstrained)

    def test_get_cuts(self):
        # this used to cause an error since
        # operator.isNumberType() on an ndarry is always true
        cuts = np.array([1.0, 2.0, 3.0])
        sfu.get_cuts(cuts)

if __name__ == "__main__":
    ut.main()

