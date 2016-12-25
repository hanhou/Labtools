#!/usr/bin/env python

import unittest as ut
import pypsignifit.psignidata as pd
import swignifit.swignifit_raw as sft

import numpy as np
import sys

def approximatedly_equal ( x, y, eps=1e-4 ):
    e = abs ( x-y ) > eps
    if not e==0:
        sys.stderr.write ( "%g != %g\n" % (x,y) )
    return e


class TestPsiInference ( ut.TestCase ):
    def setUp ( self ):
        self.pinf = pd.PsiInference ()
        self.pinf.estimate = [2.,1.,.02]
    def test_evaluate ( self ):
        evaluated = self.pinf.evaluate ( [1,2],[2.,1.,.02] )
        self.assertEqual ( approximatedly_equal ( evaluated[0], 0.62909188225759771), 0 )
        self.assertEqual ( approximatedly_equal ( evaluated[1], 0.74), 0 )
        evaluated = self.pinf.evaluate ( [1,2] )
        self.assertEqual ( approximatedly_equal ( evaluated[0], 0.62909188225759771), 0 )
        self.assertEqual ( approximatedly_equal ( evaluated[1], 0.74), 0 )
    def test_getThres ( self ):
        self.assertRaises ( NotImplementedError, self.pinf.getThres, .5 )
        self.pinf.data = [[1,2,3]]
        evaluated = self.pinf.getThres ( .5 )
        self.assertEqual ( approximatedly_equal ( evaluated, 2.), 0 )
        evaluated = self.pinf.getThres ( .3 )
        self.assertEqual ( approximatedly_equal ( evaluated, 1.1527021396127963), 0 )
    def test_repr ( self ):
        self.assertEqual ( self.pinf.__repr__(), "< PsiInference object >" )
    def test_properties ( self ):
        self.assertEqual ( self.pinf.desc, 'sigmoid: logistic\ncore: ab\nnAFC: 2' )
        self.assertEqual ( self.pinf.label, "Psychometric function fit" )
        self.assertEqual ( self.pinf.color, "b" )
        self.assertEqual ( self.pinf.linestyle, "-" )
        self.assertEqual ( self.pinf.marker, "o" )
        self.assertEqual ( self.pinf.linewidth, 1 )

class TestBootstrapInference ( ut.TestCase ):
    def setUp ( self ):
        sft.setSeed(0)
        nafc = 2
        stimulus_intensities = [0.0,2.0,4.0,6.0,8.0,10.0]
        number_of_correct = [34,32,40,48,50,48]
        number_of_trials  = [50]*len(stimulus_intensities)
        data = zip(stimulus_intensities,number_of_correct,number_of_trials)
        self.parametric    = pd.BootstrapInference ( data, priors=("","","Beta(2,30)"), parametric=True )
        self.nonparametric = pd.BootstrapInference ( data, priors=("","","Beta(2,30)"), parametric=False )
    def test_map ( self ):
        map1 = self.parametric.estimate
        map2 = self.nonparametric.estimate
        should_be = [ 2.7373, 1.40406, 0.020320093764199146 ]
        for val1,val2,val3 in zip(map1,map2,should_be):
            self.assertEqual ( approximatedly_equal ( val1, val2 ), 0 )
            self.assertEqual ( approximatedly_equal ( val1, val3 ), 0 )

    def test_boots ( self ):
        self.parametric.sample ()
        self.nonparametric.sample ()
        parci = self.parametric.getCI(1)
        nprci = self.nonparametric.getCI(1)
        self.assertEqual ( approximatedly_equal ( parci[0], 1.69 ), 0 )
        self.assertEqual ( approximatedly_equal ( parci[1], 3.8539 ), 0 )
        self.assertEqual ( approximatedly_equal ( nprci[0], 1.11463 ), 0 )
        self.assertEqual ( approximatedly_equal ( nprci[1], 4.05597 ), 0 )

        self.assertEqual ( self.parametric.nsamples, 2000 )
        self.assertEqual ( self.nonparametric.nsamples, 2000 )

        self.assertEqual ( approximatedly_equal ( self.parametric.deviance,  8.1689126711025022 ), 0 )
        self.assertEqual ( approximatedly_equal ( self.nonparametric.deviance,  8.1689126711025022 ), 0 )

    def test_sensitivity ( self ):
        self.parametric.sensitivity_analysis (verbose=False)
        parci = [ 1.5905, 3.87779 ]
        extci = self.parametric.getCI(1)
        for par,ext in zip(parci,extci):
            self.assertEqual ( approximatedly_equal ( par, ext ), 0 )

    def test_keywordhandling ( self ):
        self.assertRaises ( ValueError, pd.BootstrapInference, self.parametric.data, shape="logistic" )

class TestBayesInference ( ut.TestCase ):
    def setUp ( self ):
        sft.setSeed(0)
        nafc = 2
        stimulus_intensities = [0.0,2.0,4.0,6.0,8.0,10.0]
        number_of_correct = [34,32,40,48,50,48]
        number_of_trials  = [50]*len(stimulus_intensities)
        data = zip(stimulus_intensities,number_of_correct,number_of_trials)
        self.mcmc = pd.BayesInference ( data, priors=("Gauss(0,100)","Gamma(1.01,200)","Beta(2,30)") )
    def test_all ( self ):
        mapest = self.mcmc.mapestimate
        meanest = self.mcmc.estimate
        map_target = [ 2.73973931, 6.15554732, 0.02034599]
        mean_target =[ 2.64938, 6.44707, 0.027297]
        steps_made = self.mcmc._steps
        steps = [ 0.726551, 2.45564, 0.013264]
        burnin = 0
        thinning = 1
        nsamples = 600

        for k in xrange ( 3 ):
            self.assertEqual ( approximatedly_equal ( mapest[k],     map_target[k] ),  0 )
            self.assertEqual ( approximatedly_equal ( meanest[k],    mean_target[k] ), 0 )
            self.assertEqual ( approximatedly_equal ( steps_made[k], steps[k] ),       0 )

        self.assertEqual ( approximatedly_equal ( self.mcmc.bayesian_p(),0.126667 ), 0 )
        self.assertEqual ( burnin, self.mcmc.burnin )
        self.assertEqual ( thinning, self.mcmc.thin )
        self.assertEqual ( nsamples, self.mcmc.nsamples )

        target_rpd = -0.0244598
        target_rkd = -0.362064
        self.assertEqual ( approximatedly_equal ( self.mcmc.Rpd, target_rpd ), 0 )
        self.assertEqual ( approximatedly_equal ( self.mcmc.Rkd, target_rkd ), 0 )

        target_dr = (1.64122, -0.675137, -0.709666, 0.925372, 2.00248, -0.376286)
        target_thres =  [ 1.03761, 2.64938, 4.26115]
        target_deviance = 8.66087
        for dr,tdr in zip ( self.mcmc.devianceresiduals, target_dr ):
            self.assertEqual ( approximatedly_equal ( dr, tdr ), 0 )
        for th,tth in zip ( self.mcmc.thres, target_thres ):
            self.assertEqual ( approximatedly_equal ( th, tth ), 0 )
        self.assertEqual ( approximatedly_equal ( self.mcmc.deviance, target_deviance ), 0 )

        # Randomly check single samples
        target_mcRpd = [ -0.162011225773 , -0.658780099748 , 0.142264200236 ]
        target_mcRkd = [ -0.519220509587 , -0.969883465483 , -0.199933951214 ]
        target_mcthres = [ [ 1.20085248224 , 2.73973931207 , 4.27862614189 ], [ 2.87682033438 , 3.73674348349 , 4.5966666326 ], [ 0.70405560915 , 2.68052230561 , 4.65698900207 ] ]
        indices = [10,50,100]
        for k in xrange ( 3 ):
            self.assertEqual ( approximatedly_equal ( self.mcmc.mcRpd[indices[k]], target_mcRpd[k] ), 0 )
            self.assertEqual ( approximatedly_equal ( self.mcmc.mcRkd[indices[k]], target_mcRkd[k] ), 0 )
            for l in xrange ( 3 ):
                self.assertEqual ( approximatedly_equal ( self.mcmc.mcthres[indices[k]][l], target_mcthres[k][l] ), 0 )
    def test_keywordhandling ( self ):
        self.assertRaises ( ValueError, pd.BayesInference, self.mcmc.data, shape="logistic" )

class Testcheck_kwargs ( ut.TestCase ):
    def test_checking ( self ):
        self.assertRaises ( ValueError, pd.check_kwargs, {"test": 1}, "Some text" )
        docstr = """:Parameters:
        *test* :
            useless documentation
        *anotherprm* :
            dummy parameter
        *prm0* :
            should work with numbers, too
        *prmCamelCase* :
            should work with capitals, too
        *prm_with_underscores* :
            and should work with underscores
        *prmtype* : float
            and should work with type specifications
        """
        self.assertEqual ( 0,               pd.check_kwargs ( {"test": 1},                           docstr ) )
        self.assertEqual ( "notavailable",  pd.check_kwargs ( {"notavailable": 1},                   docstr ) )
        self.assertEqual ( 0,               pd.check_kwargs ( {"prm0": 1},                           docstr ) )
        self.assertEqual ( 0,               pd.check_kwargs ( {"prmCamelCase": 1},                   docstr ) )
        self.assertEqual ( 0,               pd.check_kwargs ( {"prm_with_underscores": 1},           docstr ) )
        self.assertEqual ( 0,               pd.check_kwargs ( {"prmtype": 1},                        docstr ) )
        self.assertEqual ( "notin1",        pd.check_kwargs ( {"notin1": 1, "test": 1},              docstr ) )
        self.assertEqual ( "notin1",        pd.check_kwargs ( {"test": 1, "notin1": 1},              docstr ) )
        self.assertEqual ( "notin1",        pd.check_kwargs ( {"test": 1, "notin1": 1, "notin2": 1}, docstr ) )

if __name__ == "__main__":
    ut.main()
    # suite = ut.TestLoader().loadTestsFromTestCase(TestBayesInference)
    # ut.TextTestRunner ().run(suite)
