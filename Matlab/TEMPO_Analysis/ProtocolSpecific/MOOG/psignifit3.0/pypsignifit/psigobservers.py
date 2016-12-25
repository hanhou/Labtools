#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__docformat__ = "restructuredtext"

import numpy as N
from scipy import stats
import pypsignifit
interface = pypsignifit.interface

from psignidata import Property

__doc__ = """
If we want to know how well a psychometric function can describe an observers behavior, we
may want to simulate an observer. This module implements a number of simulated observers.
The basic observer does not violate any assumptions. However more elaborated observers
violate some of the assumptions that are typical when fitting psychometric functions.
"""

__all__ = ["Observer","LinearSystemLearner","CriterionSettingObserver","BetaBinomialObserver"]

class Observer ( object ):
    def __init__ ( self, *params, **model ):
        """A stationary binomial observer

        This is the observer, we all want: No interdependencies between trials, no learning,
        no fluctuations in attention or motivation. Perfectly binomial responses in accordance
        with the psychometric function shape you supply.

        :Parameters:
            *params* :
                a list of parameters in the model. For nAFC tasks the parameters are a,b,lapse
                for a Yes/No task the parameters are a,b,lapse,guess
            *model* :
                a list of keyword arguments specifying the model. These are the same as in
                the psignidata module

        :Example:
        >>> O = Observer ( 4, 1, .02 )
        >>> O.seed ( 0 )
        >>> O.DoATrial ( 3 )
        1
        >>> O.DoABlock ( 4, 40 )
        28
        >>> O.DoABlock ( 6, 40 )
        37
        >>> O.DoAnExperiment ( [2,4,6], 50 )
        [[2, 27, 50], [4, 38, 50], [6, 46, 50]]
        >>> O.data
        [[3, 1, 1], [4, 28, 40], [6, 37, 40], [2, 27, 50], [4, 38, 50], [6, 46, 50]]

        :Example:
        >>> model ={"sigmoid" : "gumbel_r", "core" : "mw01", "nafc" : 2}
        >>> correct = [0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
        >>> ob = Observer(4, 4, 0.05, **model)
        >>> levels = ob.getlevels(correct)
        >>> data = ob.DoAnExperiment(levels, ntrials=50)
        """
        if model.setdefault( "nafc", 2 ) == 1:
            self.a,self.b,self.lapse,self.guess = params
        else:
            self.a,self.b,self.lapse = params
            self.guess = 1./model["nafc"]

        # Make sure, a,b,lapse and guess are floats
        self.a = float(self.a)
        self.b = float(self.b)
        self.lapse = float(self.lapse)
        self.guess = float(self.guess)

        self.model = {
                "sigmoid": model.setdefault ( "sigmoid", "logistic" ),
                "core":    model.setdefault ( "core",    "ab" ),
                "nafc":    model.setdefault ( "nafc",    2 )
                }
        self.data = []

    def DoATrial ( self, stimulus_intensity=1 ):
        """Simulate a single trial

        :Parameters:
            *stimulus_intensity* :
                stimulus intensity to be presented

        :Output:
            The response in the trail (1/0 coding for Yes/No in Yes/No-Tasks or for
            Correct/Incorrect in nAFC tasks)
        """
        prob = float( interface.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

        resp = int ( N.random.rand()<prob )

        if len(self.data) == 0 or not stimulus_intensity==self.data[-1][0]:
            self.data.append ( [stimulus_intensity,resp,1] )
        else:
            self.data[-1][1] += resp
            self.data[-1][2] += 1

        return resp

    def DoABlock ( self, stimulus_intensity=1, ntrials=50 ):
        """Simulate a block of trials

        :Parameters:
            *stimulus_intensity* :
                stimulus intensity to be presented
            *ntrials* :
                number of trials in the block

        :Output:
            The number of Yes-responses (in Yes/No) or the number of correct responses (nAFC)
        """
        prob = float( interface.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

        resp = N.random.binomial ( ntrials, prob )

        self.data.append ( [stimulus_intensity, resp, ntrials] )

        return resp

    def DoAnExperiment ( self, stimulus_intensities, ntrials=50 ):
        """Simulate a whole experiment

        :Parameters:
            *stimulus_intensities* :
                a sequence of stimulus intensities in the order they should be presented
            *ntrials* :
                Either an integer or a sequence. If this is an integer, it is interpreted as a constant
                number of trials that is presented in each trial, otherwise, the sequence is expected to
                contain a number of trials for each trials.

        :Output:
            A list of lists. Each element of the list has the structure
            [stimulus_intensity, number_of_correct_or_yes_responses,number_of_trials]
        """
        if isinstance ( ntrials, int ):
            ntrials = [ntrials] * len(stimulus_intensities)

        data = []

        for s,n in zip ( stimulus_intensities, ntrials ):
            data.append ( [s, self.DoABlock ( s, n ), n] )

        return data

    def seed ( self, seed ):
        """Seed the underlying random number generator to a defined value"""
        N.random.seed ( seed )

    def __str__ ( self ):
        if self.model["nafc"] < 2:
            return "< Observer a=%g,b=%g,lapse=%g,guess=%g,core=%s,sigmoid=%s >" % (self.a,self.b,self.lapse,self.guess,self.model["core"],self.model["sigmoid"])
        else:
            return "< Observer a=%g,b=%g,lapse=%g,core=%s,sigmoid=%s >" % (self.a,self.b,self.lapse,self.model["core"],self.model["sigmoid"])

    def getlevels ( self, cuts ):
        """Determine stimulus levels that correspond to predefinde levels of performance"""
        return interface.diagnostics ( [[1,2,3]], self.params, self.model["nafc"], self.model["sigmoid"], self.model["core"], cuts )[3]

    def evaluate ( self, stimulus_intensities ):
        """Evaluate the psychometric function

        :Parameters:
            *stimulus_intensities* :
                stimulus intensities at which the psychometric function
                should be evaluated.
        """
        return interface.diagnostics ( stimulus_intensities, self.params,
                sigmoid=self.model["sigmoid"], core = self.model["core"], nafc=self.model["nafc"] )

    @Property
    def params ():
        "parameters of the model"
        def fget ( self ):
            if self.model["nafc"] < 2:
                return [self.a,self.b,self.lapse,self.guess]
            else:
                return [self.a,self.b,self.lapse]

    @Property
    def thres ():
        "determine 50%% of the model"
        def fget ( self ):
            return float(interface.diagnostics ( [], self.params, cuts=0.5, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3])

class LinearSystemLearner ( Observer ):
    def __init__( self, *params, **model ):
        """A nonstationary observer that learns like a linear system in one or more parameters

        For this observer, the parameters of the psychometric function change: The parameters given
        to the constructor are initial values. In addition, a dictionary is given containing the learning
        parameters.

        The constructor requires a dictionary that describes the changes of the parameters. The dictionary
        has keys 'a', 'b', and/or 'lapse'. For each key, there are two parameters:
        the rate of decay and the asymptote where the parameter converges. Thus
        {'a' : ( 40, 3)}
        means that after every trial, the a parameter is modified according to
        a -= (a-3)/4

        :Parameters:
            *params* :
                parameters of the model a, b, lapse and in Yes/No tasks also guess. In addition, a dictionary
                is required to describe the changes of the parameters (see above)
            *model* :
                a list of keyword arguments to describe the model. These are the same as in the psignidata
                module.

        :Example:
        >>> O = LinearSystemLearner ( 7, 2, .02, {'a': (40,3)} )
        >>> O.seed ( 0 )
        >>> O.a
        7.0
        >>> O.DoATrial ( 3 )
        1
        >>> O.DoABlock ( 4, 50 )
        34
        >>> O.DoAnExperiment ( [4,2,8,10,6], 50 )
        [[4, 43, 50], [2, 32, 50], [8, 48, 50], [10, 49, 50], [6, 38, 50]]
        >>> O.data
        [[3, 1, 1], [4, 34, 50], [4, 43, 50], [2, 32, 50], [8, 48, 50], [10, 49, 50], [6, 38, 50]]
        >>> O.a
        3.0019608723226945
        """
        # TODO: An other example for sharper discriminations (b) and for fatigue (lapse)
        Observer.__init__ ( self, *(params[:-1]), **model )
        self.__parstring = "a=%g,b=%g,lapse=%g" % (self.a,self.b,self.lapse)
        if self.model["nafc"] < 2:
            self.__parstring += ",guess=%g" % (self.guess,)
        self.__parstring += ",core=%s,sigmoid=%s" % (self.model["core"],self.model["sigmoid"])
        self.learn = params[-1]

    def DoATrial ( self, stimulus_intensity=1 ):
        """Simulate a single trial with learning

        :Parameters:
            *stimulus_intensity* :
                intensity of the presented stimulus

        :Output:
            either 1 or 0 indicating Yes/No in Yes/No-Tasks or Correct/Incorrect in nAFC
        """
        # Get the response:
        resp = Observer.DoATrial ( self, stimulus_intensity )

        # Modify parameters:
        for prm,vals in self.learn.iteritems():
            if prm in ["a","m"]:
                self.a -= (self.a-vals[1])/vals[0]
            elif prm in ["b","w"]:
                self.b -= (self.b-vals[1])/vals[0]
            elif prm=="lapse":
                self.lapse -= (self.lapse-vals[1])/vals[0]
            else:
                # This should issue a warning
                raise ValueError, "Trying to modify parameter "+prm+" which does not make sense"

        return resp

    def DoABlock ( self, stimulus_intensity=1, ntrials=50 ):
        """Simulate a block of trials with learning

        :Parameters:
            *stimulus_intenstiy* :
                intensity of the presented stimulus
            *ntrials* :
                number of repetitions of the stimulus

        :Output:
            the number of Yes-responses in a Yes/No-task or the number of correct responses in nAFC
        """
        self.data.append ( [stimulus_intensity, 0, 0] )
        resp = 0
        for k in xrange(ntrials):
            resp += self.DoATrial ( stimulus_intensity )
        return resp

    def __str__ ( self ):
        return "< LinearSystemLearner %s,learning: %s >" % (self.__parstring,self.learn)

class CriterionSettingObserver ( Observer ):
    def __init__ ( self, *params, **model ):
        """A nonstationary observer that recalibrates its decision process according to CST

        Criterion setting theory (CST) is an extension to signal detection theory that was originally
        proposed by Treisman & Williams (1984): each isolated psychophysical decision is described in
        terms of signal detection theory. In addition, the criterion used by the observer is updated
        from trial to trial. There are two processes describing the updating of the criterion:

        1. *stabilization* : The criterion is moved into the direction of new incoming stimuli. This
            results in the criterion being placed roughly in the middle of the stimulus stream. Thus,
            stabilization can be said to maximize information transmission.
        2. *tracking* : The criterion is placed in a direction that increases repetitions of the same
            response. This strategy can be considered as a prior for an unchanging outside world.

        Both, stabilization and tracking are described in the form of linearly decreasing traces.
        If Ec(j) ist the effective criterion on trial j and E(j) is the sensory input on trial j, then
        this will result in a stabilization trace of the form Ds*( E(j) - Ec(j) ). After one trial,
        this trace will have decreased on the next trial j+1 to Ds*( E(j) - Ec(j) ) - ds. Thus,
        Ec(j+1) = Ec(j) + Ds*( E(j) - Ec(j))-ds. Similarly for tracking, there are two parameters Dr
        and dr describing the size of an indicator trace directly after a response has been made and
        describing the decrease of an indicator trace from one trial to another.

        :Parameters:
            *params* :
                a,b,lapse parameters of the psychometric function describing the interrelation between
                        stimulus intensity and sensitivity (where sensitivity is given in fraction of
                        correct responses in 2AFC).
                E0      the reference criterion, i.e. Ec(0) := E0
                Ds,ds   parameters of the stabilization mechanism
                Dr,dr   parameters of the tracking mechanism

        :Example:
        >>> O = CriterionSettingObserver ( 4., 1., 0.02, 0.5, 0.1, 0.005, 0.1, 0.01 )
        >>> O.seed ( 0 )
        >>> O.DoATrial ( 1 )
        0
        >>> O.DoABlock ( 3, 30 )
        24
        >>> O.DoABlock ( 3, 30, 0 )
        (5, 11, 8, 19)
        """
        # Make sure, the parameters of the psychometric function are interpreted as 2AFC
        model["nafc"] = 2
        # Initialize the Observer
        Observer.__init__ ( self, *(params[:3]), **model )

        # Store CST parameters
        self.E0, self.Ds, self.ds, self.Dr, self.dr = params[3:]

        self.Straces = []
        self.Ttraces = []

    def DoATrial ( self, stimulus_intensity=1, nAFC=2 ):
        """Simulate a single trial with criterion setting

        :Parameters:
            *stimulus_intensity* :
                intensity of the stimulus the is presented
            *nAFC* :
                number of alternatives presented / task to be performed.
                - nAFC=-1 Yes/No task, stimulus presence is random.
                - nAFC=0  Yes/No task, noise only trial
                - nAFC=1  Yes/No task, signal+noise trial
                - nAFC>1  nAFC task, number of alternatives

        :Output:
            A number 1 or 0, depending on the task, this means correct/incorrect (nAFC)
            or signal present, signal absent (Yes/No)
        """
        prob = float( interface.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=2 ) )
        d = stats.norm.ppf(prob) - stats.norm.ppf(1-prob)

        if nAFC==-1:
            if N.random.rand() < 0.5:
                stim = [0]
            else:
                stim = [d]
        elif nAFC==0:
            stim = [0]
        elif nAFC==1:
            stim = [d]
        else:
            stim = N.random.multinomial ( 1, [self.guess]*nAFC ) * d
        maxx,maxi = -100,-1

        dat = []

        # The response traces are updated only once per trial
        Ec = self.E0 + self.__updateTrace ( self.Ttraces, self.dr )
        for i,s in enumerate(stim):
            Ej = s + N.random.randn()
            # The remaining traces are updated after every stimulus
            Ecs = Ec + self.__updateTrace ( self.Straces, self.ds )
            x = Ej - Ecs

            # Store
            self.Straces.append( self.Ds*x )
            dat.append ( (Ej, Ecs) )

            if x>maxx:
                maxx = x
                maxi = i

        # self.data.append ( dat )
        self.data.append ( Ecs )

        # Decision rule depends on task
        if nAFC < 2:
            if x > 0:
                self.Ttraces.append ( self.Dr )
                return 1
            else:
                self.Ttraces.append( -self.Dr )
                return 0
        else:
            if stim[maxi] == 0:
                return 0
            else:
                return 1

    def DoABlock ( self, stimulus_intensity=1, ntrials=50, nAFC=2 ):
        """Simulate a whole block of trials with criterion setting

        :Parameters:
            *stimulus_intensity* :
                intensity of the stimulus presented
            *ntrials* :
                number of trials in the block
            *nAFC* :
                number of alternatives presented / task to be performed.
                - nAFC in {0,1,-1} Yes/No task, stimulus presence is random.
                - nAFC>1  nAFC task, number of alternatives

        :Output:
            Number of Correct responses (nAFC) or a tuple with
            (number of hits, number of signal trials, number of false alarms, number of noise trials)
        """
        if nAFC<2:
            return self.__DoABlock_YesNo ( stimulus_intensity, ntrials )
        else:
            return self.__DoABlock_nAFC ( stimulus_intensity, ntrials, nAFC )

    def __DoABlock_YesNo ( self, stimulus_intensity=1, ntrials=50 ):
        """Yes/No block"""
        nhits,nsignals,nfalsealarms,nnoise = 0,0,0,0
        for j in xrange ( ntrials ):
            s_or_n = N.random.rand () < 0.5
            if s_or_n == 0:
                nnoise += 1
                nfalsealarms += self.DoATrial ( stimulus_intensity, s_or_n )
            else:
                nsignals += 1
                nhits += self.DoATrial ( stimulus_intensity, s_or_n )
        return nhits,nsignals,nfalsealarms,nnoise

    def __DoABlock_nAFC ( self, stimulus_intensity=1, ntrials=50, nAFC=2 ):
        """nAFC block"""

        assert nAFC >= 2

        ncorrect = 0

        for k in xrange ( ntrials ):
            ncorrect += self.DoATrial ( stimulus_intensity, nAFC )

        return ncorrect

    def __updateTrace ( self, traces, d ):
        """sum up the traces and update"""
        Ec = 0
        for k,t in enumerate(traces):
            newtrace = t-N.sign(t)*d
            if newtrace*t < 0: # Sign change of traces ~> the trace died here
                traces.pop(k)
            else:
                traces[k] = newtrace
                Ec += newtrace
        return Ec

    def __str__ ( self ):
        return "< CriterionSettingObserver E0=%g,Ds=%g,ds=%g,Dr=%g,dr=%g >" % (self.E0,self.Ds,self.ds,self.Dr,self.dr)

class BetaBinomialObserver ( Observer ):
    def __init__ ( self, *params, **model ):
        """An overdispersed observer that otherwise follows binomial assumptions

        For every block, this observer draws a p value from a Beta distribution
        with mean given by the psychometric function. This p-value is then used
        to generate the responses of that block.

        :Parameters:
            *params* :
                a,b,lapse(,guess) are parameters for the classical psychometric function,
                m gives the dispersion of the Beta distribution. For a Beta distribution
                with parameters alpha and beta, m=alpha+beta
            *model* :
                a list of keywords to further specify the model

        :Example:
        >>> O = BetaBinomialObserver ( 4., 1., .02, 30 )
        >>> O.seed ( 0 )
        >>> O.DoABlock ( 4, 30 )
        24
        >>> O.DoABlock ( 4, 30 )
        18
        >>> O.DoAnExperiment ( [4,2,8,10,6], 50 )
        [[4, 37, 50], [2, 24, 50], [8, 47, 50], [10, 49, 50], [6, 49, 50]]
        """
        self.dispersion = params[-1]
        Observer.__init__ ( self, *(params[:-1]), **(model) )

    def DoABlock ( self, stimulus_intensity=1, ntrials=50 ):
        """Perform a block of trials

        :Parameters:
            *stimulus_intensity* :
                intensity of the presented stimulus
            *ntrials* :
                number of trials in the block

        :Output:
            number of correct responses (nAFC) or number of yes responses (Yes/No-Task)
        """
        psi = float( interface.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

        alpha = psi*self.dispersion
        beta  = (1-psi)*self.dispersion
        prob = stats.beta.rvs ( alpha, beta )

        resp = N.random.binomial ( ntrials, prob )

        self.data.append ( [stimulus_intensity, resp, ntrials] )

        return resp

    def __str__ ( self ):
        if self.model["nafc"] < 2:
            return "< BetaBinomialObserver a=%g,b=%g,lapse=%g,guess=%g,M=%g >" % (self.a,self.b,self.lapse,self.guess,self.dispersion)
        else:
            return "< BetaBinomialObserver a=%g,b=%g,lapse=%g,M=%g >" % (self.a,self.b,self.lapse,self.dispersion)


############################################################
# Utilities

def yesno2dprime ( hits, signals, falsealarms, noises ):
    """Determine dprime from hits and false alarms

    :Parameters:
        *hits* :
            number of hits
        *signals* :
            number of signal trials
        *falsealarms* :
            number of false alarm trials
        *noises* :
            number of noise only trials

    :Output:
        dprime, variance (based on bootstrap)
    """
    hrate = N.clip(float(hits)/signals,       .001, .999 )
    frate = N.clip(float(falsealarms)/noises, .001, .999 )

    hsample = N.clip ( N.random.binomial ( signals, hrate, size=5000 ).astype('d') / signals, .001, .999 )
    fsample = N.clip ( N.random.binomial ( noises,  frate, size=5000 ).astype('d') / noises , .001, .999 )

    dprime = stats.norm.ppf ( hrate ) - stats.norm.ppf ( frate )
    dprimes = stats.norm.ppf ( hsample ) - stats.norm.ppf ( fsample )

    return dprime, N.var ( dprimes )

def dprime2Pcorrect ( dprime ):
    """Convert dprime to 2AFC fraction of correct responses

    From a given d' the probability of a correct response in a 2AFC trial can be inferred as

    P(correct) = P(s>n) = P(n<s) = E[ Phi(s) ].

    The expectation is evaluated by solving the integral using trapezoid integration.
    """
    x = N.mgrid [ -20:20:1000j ]
    Phiphi = stats.norm.cdf(x)*stats.norm.pdf(x-dprime)
    return N.trapz ( Phiphi, x )

if __name__ == "__main__":
    import doctest
    doctest.testmod()
