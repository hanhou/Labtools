Specifying the experimental design
==================================

In order to estimate a guess rate, psignifit needs to know how many alternatives you are presenting participants with. 
However, there are vast differences in naming conventions for experimental designs. The way psignifit approaches the number of alternatives is stimulus based. 
Note that a large number of labs use a response based nomenclature.

If you don't set 'nafc' to a specific value, psignifit will assume that the data come from a 2 alternative-forced-choice
experiment (AFC). That means, that on each trial two stimuli were presented and the
observer knows that one and only one of these stimuli was the target stimulus. 
For example, in a contrast detection task for luminance only the target stimulus would have contained a different luminance.
Psignifit follows a classical signal-detection approach to experimental design. For all other situations you should create your Inference object with the keyword 'nafc' set to the number of stimulus alternatives that were presented. So, if your standard is always smaller or always larger than all the other alternatives 
then you use the number of alternatives to determine the value that you set in 'nafc'.


The easiest way of deciding which 'nafc' you want to set in your code is to consider the range of your psychometric function (PF) on the y-axis:

	- PF from 0 to 1 -> nafc = 1
	- PF from 0.5 to 1 -> nafc = 2
	- PF from 0.33 to 1 -> nafc = 3
	- PF from 0.25 to 1 -> nafc = 4
	- etc.

Setting the keyword 'nafc' to a value of 2 or larger results in a fixed guessing rate of 1/nafc.


If you are running an experiment that doesn't quite fit the classical AFC design, the way you think about your experiment might not correspond to the names you will use when analyzing your data in psignifit.

Two aspects that you will want to think about are:
	- How many alternatives are you presenting?
	- Is the standard stimulus the smallest of the comparisons?
 
Again, these aspects crucially influence your guess rate, but also the way you interpret your diagnostic plots. 

In many experiments we will record which response an observer chose and we will then
fit the number of "stimulus present", "stimulus left", "stimulus longer" responses (or
whatever is suitable in the present context). In those cases our guess rate does not necessarily equal 1/nafc.
For instance, in a detection task, the observer might be very conservative and virtually
never report the presence of a target for low stimulus intensities. Or, the observer might
always respond "stimulus left" if the stimulus is presented sufficiently far to the left
of a mark.

Another common approach is to present only one stimulus per trial. The observers might then have to indicate whether the target stimulus
was presented or not (typically called detection task, or yes-no-task). 
If your standard is in the middle of the stimulus intensities you will also have to change your approach. For example, in the discrimination example (ADD LINK) 
the stimulus intensities of the test stimulus vary symmetrically around the standard intensity. 

We will summarize all these designs as "yes-no designs" although the term yes-no is typically restricted to detection like tasks. 
The crucial difference between yes-no designs and forced choice designs for fitting
psychometric functions is that yes-no designs allow for arbitrarily set "guessing" rates.
In all these situations, the lower asymptote of the psychometric function will
be a free parameter by setting the keyword 'nafc' to 1. 
Note that in this case you also need to specify priors for four parameters instead of the three parameters in an nAFC experiment. 

If you would like to read more about how you specify the shape of the psychometric functions, have a look at the section on :doc:`PSYCHOMETRICFUNCTIONS`

A note on yes-no designs
-------------------------

Psignifit follows a classical signal-detection approach to experimental design. This means that the goodness-of-fit plots are intended for such designs as well. 

While you are able to use psignifit with both forced choice and yes-no designs, keep in mind that you will have to approach your results differently (this is especially important for the interpretation of non-stationarities which we talk about in the ADD LINK & SECTION!) 

In a detection experiment, we typically have two types of trials. In one case a target
stimulus is presented, in the other case, no target stimulus is presented. In terms of signal
detection theory, these two cases are called "signal+noise" and "noise only". The observer
responds to these two cases with either "yes, a signal was present" or "no, only noise
was presented". This results in four different outcomes on each trial: hits (the observer
correctly reports the presence of a signal), false alarms (the observer reports a signal
although "noise only" was presented), misses (the observer reports no signal
although a signal was presented), and correct rejections (the observer correctly reports
the absence of a signal). These different experimental outcomes are discussed in large detail
in standard books on signal detection theory [Green_and_Swets_1966]_. What matters with respect
to fitting psychometric functions is that depending on the strategy, both, hit rate as well
as the correct response rate might change. So which of these two should be fitted with a
psychometric function? There are two general objectives that an optimal observer could
follow in a yes-no task.

1. maximize the hit rate while keeping a fixed false alarm rate. In this case, we would
   like to fit the hit rate with a psychometric function (the false alarm rate is
   constant anyhow). Thus, if we want to fit the hit rate with a psychometric function,
   we should check that the observers maintained a more of less fixed false alarm
   rate. A future release of psignifit will contain more formal tools for this check.
2. maximize the number of correct responses. In this case, we would like to fit
   the correct response rate with a psychometric function. To check that an observer
   really uses this strategy, we should check that the false alarm rate decreases with
   the hit rate.

Furthermore, signal detection theory also offers a number of criterion free discriminability parameters,
like area under the ROC curve and the famous d' index. However, these indices can not generally
be assumed to have binomial variance (or anything similar to that). Therefore, psignifit
does not attempt to fit such data.


References
----------

.. [Green_and_Swets_1966] Green, DM and Swets, JA (1966): Signal Detection Theory and
    Psychophysics. New York: Wiley.
