==============
psignifit TODO
==============

Python interface (pypsignifit)
------------------------------

- enable fixing the guessing rate to an arbitrary value (? which?)        <- This is actually possible by putting a strong prior on the guessing rate...
- what is going on with type checking in psigniplot.GoodnessOfFit?
- add a dictionary to the data objects to make them know what they are (stimulus intensity, ...)
- ROC curves for 1AFC

- Blob size reduction for PF
- Deviance Bootstrap only one sided test
- Error bars on PF. Simply streched bars? (Not sure what is meant)
- 68 Confidence Interval
- In 3d/4d, we could also integrate the posterior without monte carlo methods. Advantage: We get a somehow analytical formula for the posterior. Disadvantage: Might still take quite long and complicates things like CI estimation. In addition, it would be difficult to estimate the errors originating from the fact that we approximate R by a compact interval.
- Plot MCMC diagnostics  for all parameters (Not sure what is meant)
- Plot the blue line on top of all the white ones, in the diagnostic plot.
- Alternative view: Shaded region of the posterior instead of 20 sample PFs
- In diagnostic plots, evtl. ad a message on top (what does this mean)
- Replace "model correction" with psi(x) or psi(stimulus intensity) (for people with less elaborated statistics background)
- Influence for Slope and Influence for Threshold seperately instead of using an
  aggregation in the Plot for the Influential Observers.
- directly use swignifit

Python interface (swignifit/psypy)
----------------------------------

- speed!

C++ engine
----------

- noninformative Jeffreys Priors
- seed the random number generator
- Goodness of fit
- Unit test with gamma as a free parameter ( works in principle, but only very rough correspondence to psignifit -- should be stricter)
- Bugs with different cores -> no convergence, exceptionally high deviance with ok looking fit, ... This is only due to the logCore-Kernel that returns -inf for contrasts of 0
- codes in sigmoids should be constants with ' const unsigned int $VARIBALE'
- add 'void' to the clone method calls
- In PsiMCList, getdeviance and setdeviance should be getDeviance and
  setDeviance
- Review the documentation for BootstrapList
- in PsiSigmoid getcode() should be getCode()
- Use pure virtual functions in PsiSigmoid, PsiCore, and PsiPrior, instead of
  NotImplemetedError
- Make proper copy constructors for rng.h instead of relying on implicit ones

- When using [] to access elements from the vector class, we need to manually
  check that the index value is within range. Alternatively we could use the
  'at()' method, which is index safe, and throws an 'std::out_of_range'
  exception. Consider replacing parts of the code that do th manual checking
  with calls to 'at()'. In this case we would also benefit from the swig
  exception handling for all STL exceptions.

MayBe:
- implement fullmodel, nullmodel
- Implement alternative SimplexAlgorithm (Numerical Recipes?, gsl-Wrapper-Object?) The current one depends heavily on the initial simplex!

    * Perform some gradient based steps after simplex optimization (doesn't work good)
    * alternative: Perform some gradient based steps as a special case of simplex optimization (e.g. particularly good points are moved in the direction of the gradient?)

Documentation
-------------

- add a delete/uninstall possibility
- Makefile/setup.py call for doctests
- List _all_ dependencies in terms of packages for debian (python/numpy headers
  sphinx etc...)

- Figure for sigmoid/core philosophy

R statistical computing language
--------------------------------

- continue R interface -- maybe swig?
- Diagnostic plots: Influential, Parameter inspection
- MCMC convergence diagnostics and Raftery & Lewis stuff

Recently done
-------------

+ GoodnessOfFit looks better now
+ math symbols in documentation according to
  http://matplotlib.sourceforge.net/sampledoc/extensions.html
+ allow for the gamma=lambda prior
+ optimize automatic proposal distribution generation
+ swignifit solves psipy problems

+ negative Gamma prior
+ Influential observations and outliers for Bayes                          [OK]
+ improved search for starting values                                      [OK]
+ influential observations marked graded                                   [OK]
+ posterior predictive Rkd, Rpd                                            [OK]
+ more meaningfull errors if sample based plot functions are used before sampling [OK]
+ Inference objects take relative probabilities, too                       [OK]
+ nonparametric bootstrap                                                  [OK]
+ Sensitivity analysis                                                     [OK]
+ Add ThresholdPlot to Tutorial                                            [OK]
+ resampling of chains in BayesInference objects                           [OK]
+ Like ParameterPlot but for thresholds                                    [OK]
+ move numbers further away from the axes.                                 [OK]
+ warning message for Rpd: "Try other sigmoid!"                            [OK]
+ unit tests                                                               [OK]
+ write a number of simulated observers                                    [OK]
+ complete tutorial                                                        [OK]
+ setup.py                                                                 [OK]
+ More Sigmoids (gumbel, weibull, gauss, ...)                              [OK] at least for now
+ log-core to allow fitting data on log contrast (i.e. gumbel to weibull)  [OK]
+ unit tests for logCore and linearCore                                    [OK]
+ linear core ax+b                                                         [OK]
+ Unit test for mwCore                                                     [OK]
+ Outliers and Influential observations                                    [OK]
+ MCMC
     implement dlposteri und dnegllikeli                                   [OK]
     check hybrid MCMC versus MH-MCMC                                      [OK]
     can we put both MCMC strategies together to have the same base class? [OK]
+ Documentation                                                            [OK]
+ pointer arithmetic for datasets                                          [OK]
+ low level Python interface
+ refactor the python toolbox to have "strict" data objects and plot functions working on top of these  [OK]
+ Convergence diagnostics for MCMC                                         [OK]
+ posterior intervals and posterior histograms for model parameters        [OK]
+ Using linalg matrix routines in leastfavourable                          [OK]
+ Don't use asymptotic values for the correlations.                        [OK] only for Rkd, Rpd seemed be be based on all blocks (Why?)
+ copy Core, Sigmoid, ... in psychometric                                  [OK] done for priors too
+ migrate to boost-python?                                                 [OK] decided to use SWIG instead
