#!/usr/bin/env python

# Unit tests for the command line interface
# These tests mainly assume that the fitting procedures work. So the tests just make sure that
# all commands can be called, that the results can be read again, and that the results are complete.

import os
import unittest as ut
import re
import numpy as np

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data2afc = zip ( x, k, n )

x = [float(2*k) for k in xrange(6)]
k = [2,21,29,42,50,49]
n = [50]*6
data1afc = zip ( x, k, n )

def writedata ( f, d ):
    for l in d:
        f.write ( "%g %d %d\n" % l )
    return None

def getvariable ( varname, datastring ):
    pattern = r"# %s\n(.*?)\n\n" % (varname,)
    m = re.search ( pattern, datastring, re.DOTALL )
    v = np.fromstring ( m.group(1), sep=' ' )
    return v

class TestCLImapestimate ( ut.TestCase ):
    def test_2afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data2afc )
        f.close()
        cmd = "psignifit-mapestimate -nafc 2 .testdata -o .testresults2"
        os.system ( cmd )
        f = open ( ".testresults2" )
        results = f.read()
        f.close()

        th = getvariable ( "thetahat", results )
        self.assertAlmostEqual ( th[0], 2.751832 )
        self.assertAlmostEqual ( th[1], 6.403953 )
        self.assertAlmostEqual ( th[2], 0.015555 )

        fisher = np.reshape( getvariable ( "fisher_info", results ), (3,3) )
        self.assertAlmostEqual ( fisher[0,0], -4.339672 )
        self.assertAlmostEqual ( fisher[0,1], -0.826837 )
        self.assertAlmostEqual ( fisher[0,2], -58.007662 )
        self.assertAlmostEqual ( fisher[1,0], -0.826837 )
        self.assertAlmostEqual ( fisher[1,1], -0.384938 )
        self.assertAlmostEqual ( fisher[1,2], -25.336736 )
        self.assertAlmostEqual ( fisher[2,0], -58.007662 )
        self.assertAlmostEqual ( fisher[2,1], -25.336736 )
        self.assertAlmostEqual ( fisher[2,2], -6349.783561 )

        thres = getvariable ( "thres", results )
        self.assertAlmostEqual ( thres[0], 1.150844 )
        self.assertAlmostEqual ( thres[1], 2.751832 )
        self.assertAlmostEqual ( thres[2], 4.352820 )

        D = float ( getvariable ( "deviance", results ) )
        self.assertAlmostEqual ( D, 8.071331 )

        os.remove ( ".testdata" )
        os.remove ( ".testresults2" )

    def test_1afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-mapestimate -nafc 1 .testdata -o .testresults1"
        os.system ( cmd )
        f = open ( ".testresults1" )
        results = f.read()
        f.close()

        th = getvariable ( "thetahat", results )
        self.assertAlmostEqual ( th[0], 3.268666 )
        self.assertAlmostEqual ( th[1], 6.050049 )
        self.assertAlmostEqual ( th[2], 0.004847 )
        self.assertAlmostEqual ( th[3], 0.021376 )

        fisher = np.reshape( getvariable ( "fisher_info", results ), (4,4) )
        self.assertAlmostEqual ( fisher[0,0], -15.996856 )
        self.assertAlmostEqual ( fisher[0,1], -0.544973 )
        self.assertAlmostEqual ( fisher[0,2], -112.664222 )
        self.assertAlmostEqual ( fisher[0,3],  62.230531 )
        self.assertAlmostEqual ( fisher[1,0], -0.544973 )
        self.assertAlmostEqual ( fisher[1,1], -1.933384 )
        self.assertAlmostEqual ( fisher[1,2], -57.822245 )
        self.assertAlmostEqual ( fisher[1,3], -8.574196 )
        self.assertAlmostEqual ( fisher[2,0], -112.664222 )
        self.assertAlmostEqual ( fisher[2,1], -57.822245 )
        self.assertAlmostEqual ( fisher[2,2], -7375.069598 )
        self.assertAlmostEqual ( fisher[2,3],  247.270163 )
        self.assertAlmostEqual ( fisher[3,0],  62.230531 )
        self.assertAlmostEqual ( fisher[3,1], -8.574196 )
        self.assertAlmostEqual ( fisher[3,2],  247.270163 )
        self.assertAlmostEqual ( fisher[3,3],  -395.088389 )

        thres = getvariable ( "thres", results )
        self.assertAlmostEqual ( thres[0], 1.756153 )
        self.assertAlmostEqual ( thres[1], 3.268666 )
        self.assertAlmostEqual ( thres[2], 4.781178 )

        D = float ( getvariable ( "deviance", results ) )
        self.assertAlmostEqual ( D, 11.162991 )

        os.remove ( ".testdata" )
        os.remove ( ".testresults1" )

    def test_1afc_e ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-mapestimate -nafc 1 -e .testdata -o .testresults1e"
        os.system ( cmd )
        f = open ( ".testresults1e" )
        results = f.read()
        f.close()

        th = getvariable ( "thetahat", results )
        self.assertAlmostEqual ( th[0], 3.210313 )
        self.assertAlmostEqual ( th[1], 6.261821)
        self.assertAlmostEqual ( th[2], 0.001765 )

        fisher = np.reshape( getvariable ( "fisher_info", results ), (3,3) )
        self.assertAlmostEqual ( fisher[0,0], -16.421624 )
        self.assertAlmostEqual ( fisher[0,1], -0.571775 )
        self.assertAlmostEqual ( fisher[0,2], -133.113719 )
        self.assertAlmostEqual ( fisher[1,0], -0.571775 )
        self.assertAlmostEqual ( fisher[1,1], -1.977617 )
        self.assertAlmostEqual ( fisher[1,2], -78.395747 )
        self.assertAlmostEqual ( fisher[2,0], -133.113719 )
        self.assertAlmostEqual ( fisher[2,1], -78.395747 )
        self.assertAlmostEqual ( fisher[2,2], -10130.614861 )

        thres = getvariable ( "thres", results )
        self.assertAlmostEqual ( thres[0], 1.644858 )
        self.assertAlmostEqual ( thres[1], 3.210313 )
        self.assertAlmostEqual ( thres[2], 4.775768 )

        D = float ( getvariable ( "deviance", results ) )
        self.assertAlmostEqual ( D, 10.616 )

        os.remove ( ".testdata" )
        os.remove ( ".testresults1e" )

class TestCLIdiagnostics ( ut.TestCase ):
    def test_2afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data2afc )
        f.close()
        cmd = "psignifit-diagnostics -nafc 2 -params '4,2,.02' .testdata -o .testdiag2"
        os.system ( cmd )
        f = open ( ".testdiag2" )
        diagnostics = f.read()
        f.close()

        prediction = np.reshape ( getvariable ( "prediction", diagnostics ), (-1,2) )
        self.assertAlmostEqual ( prediction[0,1], 0.500073 )
        self.assertAlmostEqual ( prediction[1,1], 0.505854 )
        self.assertAlmostEqual ( prediction[2,1], 0.74     )
        self.assertAlmostEqual ( prediction[3,1], 0.974146 )
        self.assertAlmostEqual ( prediction[4,1], 0.979927 )
        self.assertAlmostEqual ( prediction[5,1], 0.979999 )

        evaluated = np.reshape ( getvariable ( "pmf", diagnostics ), (-1,2) )
        self.assertAlmostEqual ( evaluated[0,1], 0.500073 )
        self.assertAlmostEqual ( evaluated[-1,1], 0.979999 )

        d = getvariable ( "devianceresiduals", diagnostics )
        d_ = [ 2.573423, 1.911003, 0.994806, -0.584292, 1.423986, -0.890531]
        for found,should in zip ( d, d_ ):
            self.assertAlmostEqual ( found, should )

        thres = getvariable ( "thres", diagnostics )
        thres_ = np.array([ 3.5,  4. ,  4.5])
        for found, should in zip ( thres, thres_ ):
            self.assertAlmostEqual ( found, should )

        D = getvariable ( "deviance", diagnostics )
        self.assertAlmostEqual ( D, 14.426254 )

        rpd = getvariable ( "rpd", diagnostics )
        self.assertAlmostEqual ( rpd, -0.801282 )
        rkd = getvariable ( "rkd", diagnostics )
        self.assertAlmostEqual ( rkd, -0.982893 )

        os.remove ( ".testdata" )
        os.remove ( ".testdiag2" )

    def test_1afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-diagnostics -nafc 1 -params '4,2,.02,.01' .testdata -o .testdiag1"
        os.system ( cmd )
        f = open ( ".testdiag1" )
        diagnostics = f.read ()
        f.close()

        prediction = np.reshape ( getvariable ( "prediction", diagnostics ), (-1,2) )
        self.assertAlmostEqual ( prediction[0,1], 0.010148 )
        self.assertAlmostEqual ( prediction[1,1], 0.021829 )
        self.assertAlmostEqual ( prediction[2,1], 0.495000 )
        self.assertAlmostEqual ( prediction[3,1], 0.968171 )
        self.assertAlmostEqual ( prediction[4,1], 0.979852 )
        self.assertAlmostEqual ( prediction[5,1], 0.979998 )

        evaluated = np.reshape ( getvariable ( "pmf", diagnostics ), (-1,2) )
        self.assertAlmostEqual ( evaluated[0,1], 0.010148 )
        self.assertAlmostEqual ( evaluated[-1,1], 0.979998 )

        d = getvariable ( "devianceresiduals", diagnostics )
        d_ =  [1.595850, 9.689173, 1.204377, -3.729350, 1.426659, 0.000092]
        for found,should in zip ( d, d_ ):
            self.assertAlmostEqual ( found, should )

        thres = getvariable ( "thres", diagnostics )
        thres_ = np.array([ 3.5,  4. ,  4.5])
        for found, should in zip ( thres, thres_ ):
            self.assertAlmostEqual ( found, should )

        D = getvariable ( "deviance", diagnostics )
        self.assertAlmostEqual ( D, 113.820741 )

        rpd = getvariable ( "rpd", diagnostics )
        self.assertAlmostEqual ( rpd, -0.699618 )
        rkd = getvariable ( "rkd", diagnostics )
        self.assertAlmostEqual ( rkd, -0.534870 )

        os.remove ( ".testdata" )
        os.remove ( ".testdiag1" )

    def test_1afc_e ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-diagnostics -nafc 1 -e -params '4,2,.02' .testdata -o .testdiag1e"
        os.system ( cmd )
        f = open ( ".testdiag1e" )
        diagnostics = f.read ()
        f.close()

        prediction = np.reshape ( getvariable ( "prediction", diagnostics ), (-1,2) )
        self.assertAlmostEqual ( prediction[0,1], 0.020146 )
        self.assertAlmostEqual ( prediction[1,1], 0.031707 )
        self.assertAlmostEqual ( prediction[2,1], 0.500000 )
        self.assertAlmostEqual ( prediction[3,1], 0.968293 )
        self.assertAlmostEqual ( prediction[4,1], 0.979854 )
        self.assertAlmostEqual ( prediction[5,1], 0.979998 )

        evaluated = np.reshape ( getvariable ( "pmf", diagnostics ), (-1,2) )
        self.assertAlmostEqual ( evaluated[0,1], 0.020146 )
        self.assertAlmostEqual ( evaluated[-1,1], 0.979998 )

        d = getvariable ( "devianceresiduals", diagnostics )
        d_ =  [ 0.882222, 8.876392, 1.133807, -3.736160, 1.426604, 0.000091]
        for found,should in zip ( d, d_ ):
            self.assertAlmostEqual ( found, should )

        thres = getvariable ( "thres", diagnostics )
        thres_ = np.array([ 3.5,  4. ,  4.5])
        for found, should in zip ( thres, thres_ ):
            self.assertAlmostEqual ( found, should )

        D = getvariable ( "deviance", diagnostics )
        self.assertAlmostEqual ( D, 96.848264 )

        rpd = getvariable ( "rpd", diagnostics )
        self.assertAlmostEqual ( rpd, -0.659141 )
        rkd = getvariable ( "rkd", diagnostics )
        self.assertAlmostEqual ( rkd, -0.494372 )

        os.remove ( ".testdata" )
        os.remove ( ".testdiag1e" )

class TestCLIbootstrap ( ut.TestCase ):
    def test_2afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data2afc )
        f.close()
        cmd = "psignifit-bootstrap -nafc 2 -nsamples 100 .testdata -o .testboots2"
        os.system ( cmd )
        f = open ( ".testboots2" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,3) )
        mcdeviance  = getvariable ( "mcdeviance",  results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslopes    = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd",       results )
        mcRkd       = getvariable ( "mcRkd",       results )
        influential = getvariable ( "influential", results )
        outliers    = getvariable ( "outliers",    results )
        bias_thres  = getvariable ( "bias_thres",  results )
        bias_slope  = getvariable ( "bias_slope",  results )
        acc_thres   = getvariable ( "acc_thres",   results )
        acc_slope   = getvariable ( "acc_slope",   results )

        data_ = np.array([[ 31.,  33.,  37.,  49.,  49.,  50.],
               [ 28.,  35.,  39.,  49.,  49.,  50.],
               [ 22.,  33.,  45.,  46.,  48.,  50.],
               [ 24.,  31.,  40.,  50.,  48.,  49.]])
        est_ = np.array([[  3.33819500e+00,   5.60458100e+00,   2.68500000e-03],
               [  3.03828900e+00,   5.70957800e+00,   5.40700000e-03],
               [  2.71095100e+00,   4.19333400e+00,   2.83910000e-02],
               [  3.02637000e+00,   6.40395300e+00,   1.55550000e-02]])
        dev_ = np.array([ 5.590425,  3.67402 ,  6.205791,  9.107808])
        thres_ = np.array([[ 1.93705 ,  3.338195,  4.739341],
               [ 1.610894,  3.038289,  4.465683],
               [ 1.662617,  2.710951,  3.759284],
               [ 1.425382,  3.02637 ,  4.627359]])
        slopes_ = np.array([[ 0.147015,  0.19602 ,  0.147015],
               [ 0.144312,  0.192416,  0.144312],
               [ 0.196493,  0.26199 ,  0.196493],
               [ 0.128664,  0.171552,  0.128664]])
        rpd_ = np.array([-0.183562,  0.07488 ,  0.388136,  0.515715])
        rkd_ = np.array([-0.314035, -0.087545,  0.077809,  0.957769])

        for k,b in enumerate([0,13,52,99]):
            self.assertAlmostEqual ( mcdeviance[b], dev_[k] )
            self.assertAlmostEqual ( mcRpd[b], rpd_[k] )
            self.assertAlmostEqual ( mcRkd[b], rkd_[k] )
            for i in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[b,i],      data_[k,i] )
            for i in xrange ( 3 ):
                self.assertAlmostEqual ( mcestimates[b,i], est_[k,i] )
                self.assertAlmostEqual ( mcthres[b,i], thres_[k,i] )
                self.assertAlmostEqual ( mcslopes[b,i], slopes_[k,i] )

        influential_ = np.array([ 0.48018 ,  0.862859,  0.356132,  0.434502,  0.357766,  0.847707])

        for i in xrange ( 6 ):
            self.assertAlmostEqual ( influential[i], influential_[i] )

        bias_thres_ = np.array([-0.368003, -0.162024,  0.086973])
        bias_slope_ = np.array([-0.476577, -0.476577, -0.476577])
        acc_thres_ = np.array([ -6.78200000e-03,  -4.49780000e-02,  -2.60209610e+01])
        acc_slope_ = np.array([ -6.78200000e-03,  -4.49780000e-02,  -2.60209610e+01])

        for i in xrange ( 3 ):
            self.assertAlmostEqual ( bias_thres[i], bias_thres_[i] )
            self.assertAlmostEqual ( bias_slope[i], bias_slope_[i] )
            self.assertAlmostEqual ( acc_thres[i], acc_thres_[i] )
            self.assertAlmostEqual ( acc_slope[i], acc_slope_[i] )

        os.remove ( ".testdata" )
        os.remove ( ".testboots2" )

    def test_1afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-bootstrap -nafc 1 -nsamples 100 .testdata -o .testboots1"
        os.system ( cmd )
        f = open ( ".testboots1" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,4) )
        mcdeviance  = getvariable ( "mcdeviance",  results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslopes    = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd",       results )
        mcRkd       = getvariable ( "mcRkd",       results )
        influential = getvariable ( "influential", results )
        outliers    = getvariable ( "outliers",    results )
        bias_thres  = getvariable ( "bias_thres",  results )
        bias_slope  = getvariable ( "bias_slope",  results )
        acc_thres   = getvariable ( "acc_thres",   results )
        acc_slope   = getvariable ( "acc_slope",   results )

        data_ = np.array([[ 11.,  14.,  26.,  47.,  49.,  50.],
               [  8.,  15.,  31.,  47.,  47.,  50.],
               [  4.,  19.,  33.,  44.,  48.,  50.],
               [  5.,  12.,  33.,  50.,  48.,  49.]])
        est_ = np.array([[  3.39125900e+00,   6.38257200e+00,   8.30000000e-05,
                  6.53660000e-02],
               [  3.26135500e+00,   6.12111300e+00,   9.10000000e-05,
                  5.42920000e-02],
               [  3.04542000e+00,   6.19016800e+00,   2.00000000e-06,
                  1.33770000e-02],
               [  3.42222900e+00,   3.14740200e+00,   2.68020000e-02,
                  9.93700000e-02]])
        dev_ = np.array([ 9.116104,  4.294998,  2.197431,  6.289036])
        thres_ = np.array([[ 1.795616,  3.391259,  4.986902],
               [ 1.731077,  3.261355,  4.791633],
               [ 1.497878,  3.04542 ,  4.592962],
               [ 2.635379,  3.422229,  4.20908 ]])
        slopes_ = np.array([[ 0.129095,  0.172127,  0.129095],
               [ 0.134609,  0.179479,  0.134609],
               [ 0.133108,  0.177477,  0.133108],
               [ 0.26179 ,  0.349054,  0.26179 ]])
        rpd_ = np.array([ 0.219986,  0.097852,  0.233106,  0.190961])
        rkd_ = np.array([ 0.123775, -0.218641, -0.076982, -0.04086 ])


        for k,b in enumerate([0,13,52,99]):
            self.assertAlmostEqual ( mcdeviance[b], dev_[k] )
            self.assertAlmostEqual ( mcRpd[b], rpd_[k] )
            self.assertAlmostEqual ( mcRkd[b], rkd_[k] )
            for i in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[b,i],      data_[k,i] )
            for i in xrange ( 3 ):
                self.assertAlmostEqual ( mcestimates[b,i], est_[k,i] )
                self.assertAlmostEqual ( mcthres[b,i], thres_[k,i] )
                self.assertAlmostEqual ( mcslopes[b,i], slopes_[k,i] )

        influential_ = np.array([ 0.827883,  0.85784 ,  0.538477,  0.302408,  0.96394 ,  0.965507])

        for i in xrange ( 6 ):
            self.assertAlmostEqual ( influential[i], influential_[i] )

        bias_thres_ = np.array([-0.368003, -0.162024,  0.062085])
        bias_slope_ = np.array([-0.561792, -0.561792, -0.561792])
        acc_thres_  = np.array([ 0.012912,  0.014571,  0.016806])
        acc_slope_  = np.array([ 0.012912,  0.014571,  0.016806])

        for i in xrange ( 3 ):
            self.assertAlmostEqual ( bias_thres[i], bias_thres_[i] )
            self.assertAlmostEqual ( bias_slope[i], bias_slope_[i] )
            self.assertAlmostEqual ( acc_thres[i], acc_thres_[i] )
            self.assertAlmostEqual ( acc_slope[i], acc_slope_[i] )

        os.remove ( ".testdata" )
        os.remove ( ".testboots1" )


    def test_1afc_e ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-bootstrap -nafc 1 -e -nsamples 100 .testdata -o .testboots1e"
        os.system ( cmd )
        f = open ( ".testboots1e" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,3) )
        mcdeviance  = getvariable ( "mcdeviance",  results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslopes    = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd",       results )
        mcRkd       = getvariable ( "mcRkd",       results )
        influential = getvariable ( "influential", results )
        outliers    = getvariable ( "outliers",    results )
        bias_thres  = getvariable ( "bias_thres",  results )
        bias_slope  = getvariable ( "bias_slope",  results )
        acc_thres   = getvariable ( "acc_thres",   results )
        acc_slope   = getvariable ( "acc_slope",   results )

        data_  = np.array([[ 10.,  14.,  26.,  47.,  49.,  50.],
               [  7.,  15.,  31.,  47.,  47.,  50.],
               [  4.,  19.,  33.,  44.,  48.,  50.],
               [  5.,  12.,  33.,  50.,  48.,  49.]])
        est_   = np.array([[  3.07800800e+00,   6.57922400e+00,   3.90000000e-05],
               [  3.03700700e+00,   6.54174800e+00,   3.70000000e-05],
               [  3.00542200e+00,   6.28976500e+00,   6.10000000e-05],
               [  3.07839800e+00,   4.12339800e+00,   2.40690000e-02]])
        dev_   = np.array([ 10.87362 ,   4.61029 ,   1.987836,   9.019021])
        thres_ = np.array([[ 1.433202,  3.078008,  4.722814],
               [ 1.40157 ,  3.037007,  4.672444],
               [ 1.432981,  3.005422,  4.577864],
               [ 2.047548,  3.078398,  4.109247]])
        slopes_= np.array([[ 0.125237,  0.166982,  0.125237],
               [ 0.125954,  0.167939,  0.125954],
               [ 0.131   ,  0.174667,  0.131   ],
               [ 0.199825,  0.266434,  0.199825]])
        rpd_   = np.array([ 0.065197,  0.137176,  0.218727,  0.006982])
        rkd_   = np.array([ 0.003441, -0.155417, -0.114861, -0.426397])

        for k,b in enumerate([0,13,52,99]):
            self.assertAlmostEqual ( mcdeviance[b], dev_[k] )
            self.assertAlmostEqual ( mcRpd[b], rpd_[k] )
            self.assertAlmostEqual ( mcRkd[b], rkd_[k] )
            for i in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[b,i],      data_[k,i] )
            for i in xrange ( 3 ):
                self.assertAlmostEqual ( mcestimates[b,i], est_[k,i] )
                self.assertAlmostEqual ( mcthres[b,i], thres_[k,i] )
                self.assertAlmostEqual ( mcslopes[b,i], slopes_[k,i] )

        influential_ = np.array([ 0.678465,  1.106128,  0.484737,  0.289692,  0.382489,  0.16962 ])

        for i in xrange ( 6 ):
            self.assertAlmostEqual ( influential[i], influential_[i] )

        bias_thres_ = np.array([-0.263612,  0.037236,  0.289397])
        bias_slope_ = np.array([-0.421668, -0.421668, -0.421668])
        acc_thres_ = np.array([-0.001037, -0.000433,  0.004147])
        acc_slope_ = np.array([-0.001037, -0.000433,  0.004147])

        for i in xrange ( 3 ):
            self.assertAlmostEqual ( bias_thres[i], bias_thres_[i] )
            self.assertAlmostEqual ( bias_slope[i], bias_slope_[i] )
            self.assertAlmostEqual ( acc_thres[i], acc_thres_[i] )
            self.assertAlmostEqual ( acc_slope[i], acc_slope_[i] )

        os.remove ( ".testdata" )
        os.remove ( ".testboots1e" )

class TestCLImcmc ( ut.TestCase ):
    def test_2afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data2afc )
        f.close()
        cmd = "psignifit-mcmc -nafc 2 -nsamples 100 .testdata -o .testmcmc2"
        os.system ( cmd )
        f = open ( ".testmcmc2" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,3) )
        mcdeviance  = getvariable ( "mcdeviance", results )
        ppdeviance  = getvariable ( "ppdeviance", results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslope     = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd", results )
        mcRkd       = getvariable ( "mcRkd", results )
        ppRpd       = getvariable ( "ppRpd", results )
        ppRkd       = getvariable ( "ppRkd", results )
        influential = getvariable ( "influential", results )

        ts = [0,13,52,99]

        # In 2afc
        # print "mcdata_ = np.%s" % (repr(mcdata[ts,:]),)
        # print "mcestimates_ = np.%s" % (repr(mcestimates[ts,:]),)
        # print "mcdeviance_ = np.%s" % (repr(mcdeviance[ts]),)
        # print "ppdeviance_ = np.%s" % (repr(ppdeviance[ts]),)
        # print "mcthres_ = np.%s" % (repr(mcthres[ts,:]),)
        # print "mcslope_ = np.%s" % (repr(mcslope[ts,:]),)
        # print "mcRpd_ = np.%s" % (repr(mcRpd[ts]),)
        # print "ppRpd_ = np.%s" % (repr(ppRpd[ts]),)
        # print "mcRkd_ = np.%s" % (repr(mcRkd[ts]),)
        # print "ppRkd_ = np.%s" % (repr(ppRkd[ts]),)
        # print "influential_ = np.%s" % (repr(influential),)

        mcdata_ = np.array([[ 29.,  32.,  40.,  49.,  49.,  50.],
               [ 29.,  35.,  41.,  50.,  50.,  49.],
               [ 30.,  32.,  43.,  46.,  50.,  50.],
               [ 25.,  28.,  43.,  49.,  47.,  50.]])
        mcestimates_ = np.array([[  2.64190300e+00,   6.35792300e+00,   1.06100000e-03],
               [  2.92947800e+00,   6.50634300e+00,   7.14700000e-03],
               [  2.92993800e+00,   6.27773900e+00,   1.49940000e-02],
               [  2.53914700e+00,   6.02901700e+00,   1.30110000e-02]])
        mcdeviance_ = np.array([ 11.469979,   8.596646,   8.173506,   8.652527])
        ppdeviance_ = np.array([ 3.522112,  9.31121 ,  5.82404 ,  9.986468])
        mcthres_ = np.array([[ 1.052422,  2.641903,  4.231384],
               [ 1.302892,  2.929478,  4.556064],
               [ 1.360503,  2.929938,  4.499373],
               [ 1.031893,  2.539147,  4.046401]])
        mcslope_ = np.array([[ 0.09321 ,  0.10437 ,  0.115816],
               [ 0.081644,  0.091848,  0.102528],
               [ 0.080642,  0.091363,  0.102685],
               [ 0.097277,  0.109657,  0.122384]])
        mcRpd_ = np.array([-0.387785, -0.291975, -0.226538, -0.255566])
        ppRpd_ = np.array([ 0.291323,  0.269911,  0.41512 ,  0.544918])
        mcRkd_ = np.array([-0.732192, -0.668121, -0.610352, -0.580557])
        ppRkd_ = np.array([ 0.131472, -0.897023, -0.470602,  0.249781])
        influential_ = np.array([ 187.606016,  113.069619,  109.988082,   83.446715,   85.701597,
                 87.319936])

        for k,s in enumerate ( ts ):
            for block in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[s,block], mcdata_[k,block] )
            for prm in xrange ( 3 ):
                self.assertAlmostEqual ( mcestimates[s,prm], mcestimates_[k,prm] )
            for cut in xrange ( 3 ):
                self.assertAlmostEqual ( mcthres[s,cut], mcthres_[k,cut] )
                self.assertAlmostEqual ( mcslope[s,cut], mcslope_[k,cut] )
            self.assertAlmostEqual ( mcdeviance[s],mcdeviance_[k] )
            self.assertAlmostEqual ( mcRpd[s],mcRpd_[k] )
            self.assertAlmostEqual ( mcRkd[s],mcRkd_[k] )
            self.assertAlmostEqual ( ppRpd[s],ppRpd_[k] )
            self.assertAlmostEqual ( ppRkd[s],ppRkd_[k] )

        for block in xrange ( 6 ):
            self.assertAlmostEqual ( influential[block], influential_[block] )

        os.remove ( ".testdata" )
        os.remove ( ".testmcmc2" )

    def test_1afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-mcmc -nafc 1 -nsamples 100 .testdata -o .testmcmc1"
        os.system ( cmd )
        f = open ( ".testmcmc1" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,4) )
        mcdeviance  = getvariable ( "mcdeviance", results )
        ppdeviance  = getvariable ( "ppdeviance", results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslope     = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd", results )
        mcRkd       = getvariable ( "mcRkd", results )
        ppRpd       = getvariable ( "ppRpd", results )
        ppRkd       = getvariable ( "ppRkd", results )
        influential = getvariable ( "influential", results )

        ts = [10,13,52,99]

        # print "In 1afc"
        # print "mcdata_ = np.%s" % (repr(mcdata[ts,:]),)
        # print "mcestimates_ = np.%s" % (repr(mcestimates[ts,:]),)
        # print "mcdeviance_ = np.%s" % (repr(mcdeviance[ts]),)
        # print "ppdeviance_ = np.%s" % (repr(ppdeviance[ts]),)
        # print "mcthres_ = np.%s" % (repr(mcthres[ts,:]),)
        # print "mcslope_ = np.%s" % (repr(mcslope[ts,:]),)
        # print "mcRpd_ = np.%s" % (repr(mcRpd[ts]),)
        # print "ppRpd_ = np.%s" % (repr(ppRpd[ts]),)
        # print "mcRkd_ = np.%s" % (repr(mcRkd[ts]),)
        # print "ppRkd_ = np.%s" % (repr(ppRkd[ts]),)
        # print "influential_ = np.%s" % (repr(influential),)

        mcdata_ = np.array([[  4.,  15.,  36.,  39.,  50.,  50.],
               [  4.,  14.,  33.,  47.,  48.,  50.],
               [  2.,  16.,  29.,  45.,  50.,  50.],
               [  2.,  14.,  32.,  46.,  49.,  50.]])
        mcestimates_ = np.array([[  3.17255800e+00,   5.94703500e+00,   6.81800000e-03,
                  2.13760000e-02],
               [  3.38775600e+00,   5.98046000e+00,   1.83400000e-03,
                  2.13760000e-02],
               [  3.54890100e+00,   5.51406200e+00,   4.38100000e-03,
                  2.13760000e-02],
               [  3.42993300e+00,   5.25737800e+00,   9.49900000e-03,
                  2.13760000e-02]])
        mcdeviance_ = np.array([ 11.314784,  11.404113,  12.97937 ,  12.815567])
        ppdeviance_ = np.array([ 10.751542,   3.978192,   7.274981,   3.450069])
        mcthres_ = np.array([[ 1.685799,  3.172558,  4.659317],
               [ 1.892641,  3.387756,  4.882871],
               [ 2.170386,  3.548901,  4.927417],
               [ 2.115589,  3.429933,  4.744278]])
        mcslope_ = np.array([[ 0.068529,  0.079079,  0.090588],
               [ 0.060576,  0.070202,  0.080833],
               [ 0.050019,  0.059275,  0.069829],
               [ 0.051162,  0.061175,  0.072682]])
        mcRpd_ = np.array([ 0.213388,  0.099306, -0.08982 , -0.09154 ])
        ppRpd_ = np.array([ 0.320027,  0.580915,  0.600779,  0.768739])
        mcRkd_ = np.array([ 0.00158 , -0.129608, -0.268318, -0.25366 ])
        ppRkd_ = np.array([-0.327134,  0.372361,  0.431476,  0.707143])
        influential_ = np.array([  95.035529,  273.533598,  118.391833,  117.362097,   83.772145,
                 54.976046])

        for k,s in enumerate ( ts ):
            for block in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[s,block], mcdata_[k,block] )
            for prm in xrange ( 4 ):
                self.assertAlmostEqual ( mcestimates[s,prm], mcestimates_[k,prm] )
            for cut in xrange ( 3 ):
                self.assertAlmostEqual ( mcthres[s,cut], mcthres_[k,cut] )
                self.assertAlmostEqual ( mcslope[s,cut], mcslope_[k,cut] )
            self.assertAlmostEqual ( mcdeviance[s],mcdeviance_[k] )
            self.assertAlmostEqual ( mcRpd[s],mcRpd_[k] )
            self.assertAlmostEqual ( mcRkd[s],mcRkd_[k] )
            self.assertAlmostEqual ( ppRpd[s],ppRpd_[k] )
            self.assertAlmostEqual ( ppRkd[s],ppRkd_[k] )

        for block in xrange ( 6 ):
            self.assertAlmostEqual ( influential[block], influential_[block] )

        os.remove ( ".testdata" )
        os.remove ( ".testmcmc1" )

    def test_1afc_e ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data2afc )
        f.close()
        cmd = "psignifit-mcmc -nafc 1 -e -nsamples 100 .testdata -o .testmcmc1e"
        os.system ( cmd )
        f = open ( ".testmcmc1e" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,3) )
        mcdeviance  = getvariable ( "mcdeviance", results )
        ppdeviance  = getvariable ( "ppdeviance", results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslope     = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd", results )
        mcRkd       = getvariable ( "mcRkd", results )
        ppRpd       = getvariable ( "ppRpd", results )
        ppRkd       = getvariable ( "ppRkd", results )
        influential = getvariable ( "influential", results )

        ts = [11,13,52,99]

        # print "In 1afc e"
        # print "mcdata_ = np.%s" % (repr(mcdata[ts,:]),)
        # print "mcestimates_ = np.%s" % (repr(mcestimates[ts,:]),)
        # print "mcdeviance_ = np.%s" % (repr(mcdeviance[ts]),)
        # print "ppdeviance_ = np.%s" % (repr(ppdeviance[ts]),)
        # print "mcthres_ = np.%s" % (repr(mcthres[ts,:]),)
        # print "mcslope_ = np.%s" % (repr(mcslope[ts,:]),)
        # print "mcRpd_ = np.%s" % (repr(mcRpd[ts]),)
        # print "ppRpd_ = np.%s" % (repr(ppRpd[ts]),)
        # print "mcRkd_ = np.%s" % (repr(mcRkd[ts]),)
        # print "ppRkd_ = np.%s" % (repr(ppRkd[ts]),)
        # print "influential_ = np.%s" % (repr(influential),)

        mcdata_ = np.array([[ 28.,  41.,  43.,  46.,  48.,  49.],
               [ 26.,  37.,  42.,  49.,  50.,  46.],
               [ 30.,  38.,  46.,  46.,  50.,  50.],
               [ 30.,  34.,  45.,  47.,  47.,  50.]])
        mcestimates_ = np.array([[ -0.132204,  10.175802,   0.0116  ],
               [ -0.132204,  10.175802,   0.0116  ],
               [ -0.374314,   9.790229,   0.019928],
               [ -0.765105,   9.541507,   0.017945]])
        mcdeviance_ = np.array([ 13.295048,  13.295048,  12.981131,  13.662284])
        ppdeviance_ = np.array([  3.740625,  11.692917,   9.792765,   5.63142 ])
        mcthres_ = np.array([[-2.676155, -0.132204,  2.411746],
               [-2.676155, -0.132204,  2.411746],
               [-2.821871, -0.374314,  2.073243],
               [-3.150482, -0.765105,  1.620272]])
        mcslope_ = np.array([[ 0.107231,  0.105976,  0.104138],
               [ 0.107231,  0.105976,  0.104138],
               [ 0.110041,  0.108003,  0.105363],
               [ 0.109071,  0.105896,  0.102188]])
        mcRpd_ = np.array([-0.256815, -0.256815, -0.144461, -0.059835])
        ppRpd_ = np.array([-0.599741,  0.04276 ,  0.318483,  0.19923 ])
        mcRkd_ = np.array([-0.420435, -0.420435, -0.267914, -0.14653 ])
        ppRkd_ = np.array([-0.615887, -0.344268, -0.475148, -0.160335])
        influential_ = np.array([ 196.257643,  170.316221,  139.005295,   91.174052,  119.386247,
                 71.220669])

        for k,s in enumerate ( ts ):
            for block in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[s,block], mcdata_[k,block] )
            for prm in xrange ( 3 ):
                self.assertAlmostEqual ( mcestimates[s,prm], mcestimates_[k,prm] )
            for cut in xrange ( 3 ):
                self.assertAlmostEqual ( mcthres[s,cut], mcthres_[k,cut] )
                self.assertAlmostEqual ( mcslope[s,cut], mcslope_[k,cut] )
            self.assertAlmostEqual ( mcdeviance[s],mcdeviance_[k] )
            self.assertAlmostEqual ( mcRpd[s],mcRpd_[k] )
            self.assertAlmostEqual ( mcRkd[s],mcRkd_[k] )
            self.assertAlmostEqual ( ppRpd[s],ppRpd_[k] )
            self.assertAlmostEqual ( ppRkd[s],ppRkd_[k] )

        for block in xrange ( 6 ):
            self.assertAlmostEqual ( influential[block], influential_[block] )

        os.remove ( ".testdata" )
        os.remove ( ".testmcmc1e" )

class TestCLIgmcmc ( ut.TestCase ):
    def test_2afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data2afc )
        f.close()
        cmd = "psignifit-bootstrap -nafc 2 -nsamples 100 .testdata -o .pilot2"
        os.system ( cmd )
        cmd = "psignifit-mcmc -nafc 2 -nsamples 100 -generic -proposal .pilot2 .testdata -o .testmcmc2"
        os.system ( cmd )
        f = open ( ".testmcmc2" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,3) )
        mcdeviance  = getvariable ( "mcdeviance", results )
        ppdeviance  = getvariable ( "ppdeviance", results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslope     = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd", results )
        mcRkd       = getvariable ( "mcRkd", results )
        ppRpd       = getvariable ( "ppRpd", results )
        ppRkd       = getvariable ( "ppRkd", results )
        influential = getvariable ( "influential", results )

        ts = [0,13,52,99]

        # print "In 2afc"
        # print "mcdata_ = np.%s" % (repr(mcdata[ts,:]),)
        # print "mcestimates_ = np.%s" % (repr(mcestimates[ts,:]),)
        # print "mcdeviance_ = np.%s" % (repr(mcdeviance[ts]),)
        # print "ppdeviance_ = np.%s" % (repr(ppdeviance[ts]),)
        # print "mcthres_ = np.%s" % (repr(mcthres[ts,:]),)
        # print "mcslope_ = np.%s" % (repr(mcslope[ts,:]),)
        # print "mcRpd_ = np.%s" % (repr(mcRpd[ts]),)
        # print "ppRpd_ = np.%s" % (repr(ppRpd[ts]),)
        # print "mcRkd_ = np.%s" % (repr(mcRkd[ts]),)
        # print "ppRkd_ = np.%s" % (repr(ppRkd[ts]),)
        # print "influential_ = np.%s" % (repr(influential),)

        mcdata_ = np.array([[ 31.,  36.,  41.,  49.,  49.,  50.],
               [ 31.,  36.,  40.,  46.,  50.,  50.],
               [ 31.,  37.,  43.,  42.,  47.,  48.],
               [ 32.,  38.,  45.,  46.,  48.,  48.]])
        mcestimates_ = np.array([[  1.90600100e+00,   6.40395300e+00,   1.55550000e-02],
               [  2.64465800e+00,   8.67314900e+00,   4.01100000e-03],
               [  1.74393300e+00,   1.58297130e+01,   2.56720000e-02],
               [  1.12192900e+00,   2.03367120e+01,   1.18100000e-03]])
        mcdeviance_ = np.array([ 11.222675,   8.888714,  18.997799,  20.212403])
        ppdeviance_ = np.array([ 5.122639,  5.528634,  3.606824,  7.732305])
        mcthres_ = np.array([[ 0.305013,  1.906001,  3.506989],
               [ 0.476371,  2.644658,  4.812945],
               [-2.213495,  1.743933,  5.701361],
               [-3.962249,  1.121929,  6.206107]])
        mcslope_ = np.array([[ 0.126225,  0.137095,  0.147162],
               [ 0.08949 ,  0.095569,  0.101446],
               [ 0.066501,  0.067373,  0.068098],
               [ 0.053545,  0.053778,  0.053934]])
        mcRpd_ = np.array([-0.087845,  0.04226 ,  0.683976,  0.7236  ])
        ppRpd_ = np.array([ 0.270714,  0.522033,  0.70528 ,  0.788184])
        mcRkd_ = np.array([-0.309592, -0.327452,  0.62385 ,  0.692216])
        ppRkd_ = np.array([ 0.203812, -0.507867,  0.696117,  0.722177])
        influential_ = np.array([ 111.878579,  208.20494 ,  104.901806,  181.676591,  238.259331,
                 79.115505])

        for k,s in enumerate ( ts ):
            for block in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[s,block], mcdata_[k,block] )
            for prm in xrange ( 3 ):
                self.assertAlmostEqual ( mcestimates[s,prm], mcestimates_[k,prm] )
            for cut in xrange ( 3 ):
                self.assertAlmostEqual ( mcthres[s,cut], mcthres_[k,cut] )
                self.assertAlmostEqual ( mcslope[s,cut], mcslope_[k,cut] )
            self.assertAlmostEqual ( mcdeviance[s],mcdeviance_[k] )
            self.assertAlmostEqual ( mcRpd[s],mcRpd_[k] )
            self.assertAlmostEqual ( mcRkd[s],mcRkd_[k] )
            self.assertAlmostEqual ( ppRpd[s],ppRpd_[k] )
            self.assertAlmostEqual ( ppRkd[s],ppRkd_[k] )

        for block in xrange ( 6 ):
            self.assertAlmostEqual ( influential[block], influential_[block] )

        os.remove ( ".testdata" )
        os.remove ( ".testmcmc2" )
        os.remove ( ".pilot2" )

    def test_1afc ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-bootstrap -nafc 1 -nsamples 100 .testdata -o .pilot1"
        os.system ( cmd )
        cmd = "psignifit-mcmc -nafc 1 -nsamples 100 .testdata -generic -proposal .pilot1 -o .testmcmc1"
        os.system ( cmd )
        f = open ( ".testmcmc1" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,4) )
        mcdeviance  = getvariable ( "mcdeviance", results )
        ppdeviance  = getvariable ( "ppdeviance", results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslope     = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd", results )
        mcRkd       = getvariable ( "mcRkd", results )
        ppRpd       = getvariable ( "ppRpd", results )
        ppRkd       = getvariable ( "ppRkd", results )
        influential = getvariable ( "influential", results )

        ts = [10,13,52,99]

        # print "In 1afc"
        # print "mcdata_ = np.%s" % (repr(mcdata[ts,:]),)
        # print "mcestimates_ = np.%s" % (repr(mcestimates[ts,:]),)
        # print "mcdeviance_ = np.%s" % (repr(mcdeviance[ts]),)
        # print "ppdeviance_ = np.%s" % (repr(ppdeviance[ts]),)
        # print "mcthres_ = np.%s" % (repr(mcthres[ts,:]),)
        # print "mcslope_ = np.%s" % (repr(mcslope[ts,:]),)
        # print "mcRpd_ = np.%s" % (repr(mcRpd[ts]),)
        # print "ppRpd_ = np.%s" % (repr(ppRpd[ts]),)
        # print "mcRkd_ = np.%s" % (repr(mcRkd[ts]),)
        # print "ppRkd_ = np.%s" % (repr(ppRkd[ts]),)
        # print "influential_ = np.%s" % (repr(influential),)

        mcdata_ = np.array([[  5.,  17.,  34.,  48.,  47.,  50.],
               [  7.,  17.,  35.,  46.,  50.,  50.],
               [  6.,  18.,  34.,  41.,  47.,  49.],
               [  3.,   8.,  26.,  47.,  48.,  50.]])
        mcestimates_ = np.array([[  2.72517100e+00,   5.81558800e+00,   2.85860000e-02,
                  2.13760000e-02],
               [  2.72517100e+00,   5.81558800e+00,   2.85860000e-02,
                  2.07710000e-02],
               [  3.42634800e+00,   5.59988900e+00,   4.34700000e-02,
                  2.00650000e-02],
               [  3.17301000e+00,   5.52024300e+00,   1.48230000e-02,
                  2.48100000e-03]])
        mcdeviance_ = np.array([ 15.504518,  15.458316,  16.423193,  11.446102])
        ppdeviance_ = np.array([  6.982574,   8.55622 ,   7.033939,  11.095793])
        mcthres_ = np.array([[ 1.271274,  2.725171,  4.179068],
               [ 1.271274,  2.725171,  4.179068],
               [ 2.026376,  3.426348,  4.82632 ],
               [ 1.792949,  3.17301 ,  4.553071]])
        mcslope_ = np.array([[ 0.087412,  0.099961,  0.113237],
               [ 0.087412,  0.099961,  0.113237],
               [ 0.055359,  0.065182,  0.076258],
               [ 0.064491,  0.075699,  0.088187]])
        mcRpd_ = np.array([ 0.453156,  0.451473,  0.295966,  0.045298])
        ppRpd_ = np.array([ 0.617999,  0.687053, -0.561595,  0.597242])
        mcRkd_ = np.array([ 0.392086,  0.389713,  0.182683, -0.144049])
        ppRkd_ = np.array([ 0.400548,  0.359163, -0.386448,  0.440935])
        influential_ = np.array([ 126.033895,  226.040904,  131.002927,  113.542458,  128.385688,
                 60.451446])

        for k,s in enumerate ( ts ):
            for block in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[s,block], mcdata_[k,block] )
            for prm in xrange ( 4 ):
                self.assertAlmostEqual ( mcestimates[s,prm], mcestimates_[k,prm] )
            for cut in xrange ( 3 ):
                self.assertAlmostEqual ( mcthres[s,cut], mcthres_[k,cut] )
                self.assertAlmostEqual ( mcslope[s,cut], mcslope_[k,cut] )
            self.assertAlmostEqual ( mcdeviance[s],mcdeviance_[k] )
            self.assertAlmostEqual ( mcRpd[s],mcRpd_[k] )
            self.assertAlmostEqual ( mcRkd[s],mcRkd_[k] )
            self.assertAlmostEqual ( ppRpd[s],ppRpd_[k] )
            self.assertAlmostEqual ( ppRkd[s],ppRkd_[k] )

        for block in xrange ( 6 ):
            self.assertAlmostEqual ( influential[block], influential_[block] )

        os.remove ( ".testdata" )
        os.remove ( ".testmcmc1" )
        os.remove ( ".pilot1" )

    def test_1afc_e ( self ):
        f = open ( ".testdata", "w" )
        writedata ( f, data1afc )
        f.close()
        f = open ( ".testd", "w" )
        writedata ( f, data1afc )
        f.close()
        cmd = "psignifit-bootstrap -nafc 1 -e -nsamples 100 .testdata -o .pilot1e"
        os.system ( cmd )
        cmd = "psignifit-mcmc -nafc 1 -e -nsamples 100 -generic -proposal .pilot1e .testdata -o .testmcmc1e"
        os.system ( cmd )
        f = open ( ".testmcmc1e" )
        results = f.read()
        f.close()

        mcdata      = np.reshape ( getvariable ( "mcdata",      results ), (-1,6) )
        mcestimates = np.reshape ( getvariable ( "mcestimates", results ), (-1,3) )
        mcdeviance  = getvariable ( "mcdeviance", results )
        ppdeviance  = getvariable ( "ppdeviance", results )
        mcthres     = np.reshape ( getvariable ( "mcthres",     results ), (-1,3) )
        mcslope     = np.reshape ( getvariable ( "mcslopes",    results ), (-1,3) )
        mcRpd       = getvariable ( "mcRpd", results )
        mcRkd       = getvariable ( "mcRkd", results )
        ppRpd       = getvariable ( "ppRpd", results )
        ppRkd       = getvariable ( "ppRkd", results )
        influential = getvariable ( "influential", results )

        ts = [11,13,52,99]

        # print "In 1afc e"
        # print "mcdata_ = np.%s" % (repr(mcdata[ts,:]),)
        # print "mcestimates_ = np.%s" % (repr(mcestimates[ts,:]),)
        # print "mcdeviance_ = np.%s" % (repr(mcdeviance[ts]),)
        # print "ppdeviance_ = np.%s" % (repr(ppdeviance[ts]),)
        # print "mcthres_ = np.%s" % (repr(mcthres[ts,:]),)
        # print "mcslope_ = np.%s" % (repr(mcslope[ts,:]),)
        # print "mcRpd_ = np.%s" % (repr(mcRpd[ts]),)
        # print "ppRpd_ = np.%s" % (repr(ppRpd[ts]),)
        # print "mcRkd_ = np.%s" % (repr(mcRkd[ts]),)
        # print "ppRkd_ = np.%s" % (repr(ppRkd[ts]),)
        # print "influential_ = np.%s" % (repr(influential),)

        mcdata_ = np.array([[  7.,  13.,  40.,  42.,  49.,  48.],
               [  7.,  15.,  31.,  44.,  50.,  50.],
               [  9.,  22.,  34.,  41.,  47.,  49.],
               [  6.,  12.,  26.,  44.,  48.,  50.]])
        mcestimates_ = np.array([[  3.26083400e+00,   6.78055500e+00,   3.95100000e-03],
               [  3.02569000e+00,   6.78055500e+00,   3.95100000e-03],
               [  3.06885200e+00,   6.86530000e+00,   2.35320000e-02],
               [  3.04677000e+00,   7.72253400e+00,   4.85300000e-03]])
        mcdeviance_ = np.array([ 11.364161,  11.641574,  14.460014,  14.385695])
        ppdeviance_ = np.array([ 11.153502,   6.482785,   3.708924,   9.635855])
        mcthres_ = np.array([[ 1.565695,  3.260834,  4.955972],
               [ 1.330551,  3.02569 ,  4.720829],
               [ 1.352527,  3.068852,  4.785177],
               [ 1.116137,  3.04677 ,  4.977404]])
        mcslope_ = np.array([[ 0.070599,  0.079499,  0.088946],
               [ 0.078954,  0.088372,  0.098207],
               [ 0.077677,  0.086845,  0.096422],
               [ 0.079982,  0.087619,  0.095388]])
        mcRpd_ = np.array([ 0.368229,  0.428669,  0.612444,  0.63752 ])
        ppRpd_ = np.array([-0.098848,  0.615304, -0.630584,  0.753689])
        mcRkd_ = np.array([ 0.167808,  0.24518 ,  0.537823,  0.543597])
        ppRkd_ = np.array([-0.239527, -0.020323, -0.475612,  0.682009])
        influential_ = np.array([ 205.008381,  147.835769,  124.361312,   97.577949,  158.19065 ,
                 54.382314])

        for k,s in enumerate ( ts ):
            for block in xrange ( 6 ):
                self.assertAlmostEqual ( mcdata[s,block], mcdata_[k,block] )
            for prm in xrange ( 3 ):
                self.assertAlmostEqual ( mcestimates[s,prm], mcestimates_[k,prm] )
            for cut in xrange ( 3 ):
                self.assertAlmostEqual ( mcthres[s,cut], mcthres_[k,cut] )
                self.assertAlmostEqual ( mcslope[s,cut], mcslope_[k,cut] )
            self.assertAlmostEqual ( mcdeviance[s],mcdeviance_[k] )
            self.assertAlmostEqual ( mcRpd[s],mcRpd_[k] )
            self.assertAlmostEqual ( mcRkd[s],mcRkd_[k] )
            self.assertAlmostEqual ( ppRpd[s],ppRpd_[k] )
            self.assertAlmostEqual ( ppRkd[s],ppRkd_[k] )

        for block in xrange ( 6 ):
            self.assertAlmostEqual ( influential[block], influential_[block] )

        os.remove ( ".testdata" )
        os.remove ( ".testmcmc1e" )
        # os.remove ( ".pilot1e" )



if __name__ == "__main__":
    ut.main()
