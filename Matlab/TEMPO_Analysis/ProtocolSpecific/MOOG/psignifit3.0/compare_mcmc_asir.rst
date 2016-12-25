Comparison: MCMC vs ASIR -- Are they equal?
===========================================

Analyze 3 datasets per condition, repeat analysis of each dataset 3 times.
Restrict sampler to 100.000 > number of mcmc samples (psignifit code!)

Questions
---------
* Numerically the same results
* Computation time
* If MCMC does not converge, do we get "something useful" from ASIR?

Numerical comparison of the following percentiles: 2.5%, 31.6%, 50%, 68.4%, 97.5%
Also compare mean and variance.

Simulation parameters
---------------------

sigmoid="logistic"
core="mw0.1"

* 20, 40, 60 trials/block
* 4, 6, 8, 12 blocks
* w_gen = 1, 2, 4
* Sampling: Best Wichmann & Hill Sampling auf der Basis m=4, w=2, for yes/no and 2afc
* yes/no, 2afc

Should give a total of 648 runs!

Sampling schemes
----------------

2afc-S7 (best 2afc): 0.34, 0.44, 0.54, 0.80, 0.90, 0.98
yes/no-y1: 0.01,0.1,0.2,0.8,0.9,0.99
