=============
psignifit API
=============

Data Types
==========
.. automodule:: pypsignifit.psignidata

.. currentmodule:: pypsignifit.psignidata

.. autosummary::

    BootstrapInference
    BayesInference

.. autoclass:: BootstrapInference
   :members:

.. autoclass:: BayesInference
   :members:

Diagnostic Plots
================

Default functions
-----------------

The following functions will be imported by default:

.. currentmodule:: pypsignifit.psigniplot

.. autosummary::

    GoodnessOfFit
    ConvergenceMCMC
    ParameterPlot
    ThresholdPlot
    plotSensitivity
    plotInfluential
    plotMultiplePMFs

.. automodule:: pypsignifit.psigniplot
   :members:

Subfunctions
------------

To get access to all plot functions in isolation, they can also be imported separately. Here is the documentation

.. currentmodule:: pypsignifit.psigniplot

.. autosummary::

    drawaxes
    plotRd
    plotHistogram
    plotPMF
    plotThres
    plotGeweke
    plotChains
    plotParameterDist

.. autofunction:: drawaxes
.. autofunction:: plotRd
.. autofunction:: plotHistogram
.. autofunction:: plotPMF
.. autofunction:: plotThres
.. autofunction:: plotGeweke
.. autofunction:: plotChains
.. autofunction:: plotParameterDist

Simulated Observers
===================

psignifit allows to simulate a number of observers to access stability of psychometric functions.

.. currentmodule:: pypsignifit.psigobservers

.. autosummary::

    Observer
    LinearSystemLearner
    CriterionSettingObserver
    BetaBinomialObserver

.. autoclass:: Observer
   :members:

.. autoclass:: LinearSystemLearner
   :members:

.. autoclass:: CriterionSettingObserver
   :members:

.. autoclass:: BetaBinomialObserver
   :members:

Errors
======

.. automodule:: pypsignifit.psignierrors
    :members:
