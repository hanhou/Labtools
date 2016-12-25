============================
A quick start to pypsignifit
============================

This document presents two example analyses of psychometric data using pypsignifit.
The :doc:`Bootstrap Tutorial </TUTORIAL_BOOTSTRAP>` explains how to fit a psychometric function using constrained maximum
likelihood as described in the papers by [Wichmann_and_Hill_2001a]_, [Wichmann_and_Hill_2001b]_. 
The :doc:`Bayes Tutorial </TUTORIAL_BAYES>` explains how to fit a psychometric function using a Bayesian approach. Parts of 
the ideas that are implemented here can be found in the paper by [Kuss_et_al_2005]_, the rest was new at the time of this writing.


Getting started
===============
To get you started with pypsignifit, open a Python interpreter and type the following:

>>> import pypsignifit as psi
>>> dir(psi)
['ASIRInference',
 'BayesInference',
 'BootstrapInference',
 'ConvergenceMCMC',
 'GoodnessOfFit',
 'ParameterPlot',
 'ThresholdPlot',
 '__builtins__',
 '__doc__',
 '__docformat__',
 '__file__',
 '__name__',
 '__package__',
 '__path__',
 '__test__',
 '__version__',
 'dump_info',
 'interface',
 'plotInfluential',
 'plotMultiplePMFs',
 'plotSensitivity',
 'psignidata',
 'psignierrors',
 'psigniplot',
 'psignipriors',
 'pygibbsit',
 'set_seed',
 'show',
 'subprocess',
 'sys',
 'version']

With the first command you import the complete functionality of the Python
module ``pypsignifit`` to your current workspace. ``dir( <module_name> )``
provides you with a list of functions and data types that come with pypsignifit.
To get help and documentation about one of these functions, you can use the
online Python help by typing ``help( <object_name> )``. For instance,

>>> help ( psi.BayesInference )

will show you the documentation of the ``BayesInference`` object.

Hint: if you would like to copy and paste the examples from this website we
recommend using the `IPython <http://ipython.scipy.org/moin/>`_ interpreter.
This has a special magic command ``%cpaste`` which ignores prefixing ``>>>``
from its input.

If you want to obtain the version identifier (for inclusion in support requests
and bug reports), type:

>>> psi.version
'snap-2011-05-17'

Experimental scenario and data format
=====================================
The data [1]_ that will be used in the following tutorials have been gathered in a 2-alternative forced-choice discrimination experiment. Observers had to discriminate between two simultaneously presented stimuli. One of them  was the original (standard) and the other one was a comparison of five different stimulus intensities which were all larger than the standard. Different comparison intensities were presented in different experimental blocks (num_of_block = 5). One block contained 50 trials (num_of_trials = 50), 25 of which contained the original and the other 25 contained one of the five different stimulus intensities. Data for all stimulus intensities were repeatedly gathered in three sessions (num_of_sess = 3). Different experimental designs are described in detail in the section `specifying your experimental design <http://psignifit.sourceforge.net/MODELSPECIFICATION.html#specifiing-the-experimental-design>`_.

We will now create our example data set for which we want to estimate a psychometric function. The data format should be a numpy array consisting of the following three columns: stimulus intensities, relative/absolute frequencies of correct (or 'yes') responses, number of observations per stimulus intensity:

    >>> import numpy as np # numpy module required
    >>> num_of_sess   = 3  # experimental parameters
    >>> num_of_block  = 5
    >>> num_of_trials = 50
    >>> stimulus_intensities = [0.021, 0.079, 0.154, 0.255,  0.30] # stimulus levels
    >>> percent_correct_1    = [0.5 ,  0.84,  0.96,  1.,   1.]     # percent correct sessions 1-3
    >>> percent_correct_2    = [0.64,  0.92,  1.  ,  0.96, 1.]
    >>> percent_correct_3    = [0.58,  0.76,  0.98,  1.,   1.]
    >>> num_observations     = [num_of_trials] * num_of_block      # observations per block
    >>> data_1 = np.c_[stimulus_intensities, percent_correct_1, num_observations]
    >>> data_2 = np.c_[stimulus_intensities, percent_correct_2, num_observations]
    >>> data_3 = np.c_[stimulus_intensities, percent_correct_3, num_observations]
    >>> data_single_sessions = np.r_[ data_1, data_2, data_3 ]       # concatenate data from all sessions

Numpy arrays data_1, data_2, data_3 summarize data from each session with each line representing a single experimental block. It is assumed that data are entered in the same sequence in which they have been acquired (often in ascending stimulus intensity as in classical signal detection tasks [Blackwell_1952]_). The last line of the code concatenates data from single sessions into a single numpy array. Again, the information about the sequence of acquisition is coded by the ordering of blocks (rows) and it will be used for the assessment of stability of performance in the :ref:`goodness of fit diagnostics <goodness_of_fit>`.


Now as you generated your data, it is time to choose whether you want to fit your psychometric function using the Bootstrap approach based on Maximum Likelihood estimation
:doc:`Maximum Likelihood Bootstrap </TUTORIAL_BOOTSTRAP>` or to chose the  :doc:`Bayesian Inference Approach </TUTORIAL_BAYES>`. 
Large scale simulations show, that especially for small datasets (n < 750) confidence intervals estimated via the Bootstrap procedure are often too small, a problem which does not occur in the Bayesian Inference approach. 

.. [1] Data courtesty of M. Maertens.
