=================================
Changes between psignifit 2 and 3
=================================

This page contains some of the changes in conversion between psignifit 2 and 3. As of now, this list is nots complete, while we are in the process of updating this list, please let us know if you find any differences that we have not added so far.

Order of Parameters
--------------------

In psignifit 2 the different parameters were used in the following oder:

alpha beta gamma lambda

Psignifit 3 changes the order of the parameters slightly when you are using the classical notation 

alpha beta lambda gamma

This is because you do not have to use gamma in all situations. When you are analysing a 2AFC task you will only have to specify 3 priors. If you are analysing a yes-no task, for example, you will have to specify your 4th prior (for gamma) as well.


Specifying the psychometric function
-------------------------------------

In psignifit 2 you specify a function from the list
    - Logistic
    - cumulative Gaussian
    - Linear
    - Gumbel
    - Weibull

This has changed. In psignifit 3 this is split into the sigmoid and the core. Below we have summarized how the old functions map onto the new sigmoid-core framework:

======================  ==================== =====================
 psignifit 2             psignifit 3 core     psignifit 3 sigmoid
======================  ==================== =====================
 Logistic                ab                    logistic
 cumulative Gaussian     ab                    gauss
 Linear                  linear                id
 Gumbel                  ab                    gumbel_l
 Weibull                 poly                  exp
======================  ==================== =====================

You should be aware though, that by having the simoid-core framework you have several other combinations that you can choose from. For a description of the different combinations that are possible in psignifit 3 have a look at the section :doc:`PSYCHOMETRICFUNCTIONS`.


Specifying priors
------------------

In psignifit 2 you could access the priors independently. For example, if you only wanted to change your lamda prior you only had to specify the lambda prior.
In psignifit 3 you have two options, either you leave all priors at their default setting or (even if you just want to change one of them) you make all priors explicit.
More information can be found in :doc:`BAYESINTRO`.
