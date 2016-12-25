#/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################
import numpy as np
import swignifit.swignifit_raw as sfr
import swignifit.utility as sfu
import operator as op

def bootstrap(data, start=None, nsamples=2000, nafc=2, sigmoid="logistic",
        core="ab", priors=None, cuts=None, parametric=True, gammaislambda=False ):
    """ Parametric bootstrap of a psychometric function.

    Parameters
    ----------

    data : A list of lists or an array of data.
        The first column should be stimulus intensity, the second column should
        be number of correct responses (in 2AFC) or number of yes- responses (in
        Yes/No), the third column should be number of trials. See also: the examples
        section below.

    start : sequence of floats of length number of model parameters
        Generating values for the bootstrap samples. If this is None, the
        generating value will be the MAP estimate. Length should be 4 for Yes/No
        and 3 for nAFC.

    nsamples : number
        Number of bootstrap samples to be drawn.

    nafc : int
        Number of alternatives for nAFC tasks. If nafc==1 a Yes/No task is
        assumed.

    sigmoid : string
        Name of the sigmoid to be fitted. Valid sigmoids include:
                logistic
                gauss
                gumbel_l
                gumbel_r
        See `swignifit.utility.available_sigmoids()` for all available sigmoids.

    core : string
        \"core\"-type of the psychometric function. Valid choices include:
                ab       (x-a)/b
                mw%g     midpoint and width
                linear   a+bx
                log      a+b log(x)
        See `swignifit.utility.available_cores()` for all available sigmoids.

    priors : sequence of strings length number of parameters
        Constraints on the likelihood estimation. These are expressed in the form of a list of
        prior names. Valid prior choices include:
                Uniform(%g,%g)
                Gauss(%g,%g)
                Beta(%g,%g)
                Gamma(%g,%g)
                nGamma(%g,%g)
                if an invalid prior or `None` is selected, no constraints are imposed at all.
        See `swignifit.utility.available_priors()` for all available sigmoids.

    cuts : a single number or a sequence of numbers.
        Cuts indicating the performances that should be considered 'threshold'
        performances. This means that in a 2AFC task, cuts==0.5 the 'threshold'
        is somewhere around 75%% correct performance, depending on the lapse
        rate parametric boolean to indicate whether or not the bootstrap
        procedure should be parametric or not.

    parametric : boolean
        If `True` do parametric, otherwise do a non-parametric bootstrap.

    gammaislambda : boolean
        Set the gamma == lambda prior.

    Returns
    -------

    (samples,estimates,deviance,
    threshold, th_bias, th_acceleration,
    slope, slope_bias, slope_accelerateion
    Rkd,Rpd,outliers,influential)

    samples : numpy array, shape: (nsamples, nblocks)
        the bootstrap sampled data

    estimates : numpy array, shape: (nsamples, nblocks)
        estimated parameters associated with the data sets

    deviance : numpy array, length: nsamples
        deviances for the bootstraped datasets

    threshold : numpy array, shape: (nsamples, ncuts)
        thresholds/cuts for each bootstraped datasets

    th_bias : numpy array, shape: (ncuts)
        the bias term associated with the threshold

    th_acc : numpy array, shape: (ncuts)
        the acceleration constant associated with the threshold

    slope : numpy array, shape: (nsamples, ncuts)
        slope at each cuts for each bootstraped datasets

    sl_bias : numpy array, shape: (ncuts)
        bias term associated with the slope

    sl_acc : numpy array, shape: (ncuts)
        acceleration term associated with the slope

    Rkd : numpy array, length: nsamples
        correlations between block index and deviance residuals

    Rpd : numpy array, length: nsamples
        correlations between model prediction and deviance residuals

    outliers : numpy array of booleans, length nblocks
        points that are outliers

    influential : numpy array of booleans, length nblocks
        points that are influential observations

    Example
    -------
    >>> x = [float(2*k) for k in xrange(6)]
    >>> k = [34,32,40,48,50,48]
    >>> n = [50]*6
    >>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    >>> priors = ('flat','flat','Uniform(0,0.1)')
    >>> samples,est,D,thres,thbias,thacc,slope,slbias,slacc,Rkd,Rpd,out,influ \
            = bootstrap(d,nsamples=2000,priors=priors)
    >>> np.mean(est[:,0])
    2.7547034408466811
    >>> mean(est[:,1])
    1.4057297989923003

    """
    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors, gammaislambda=gammaislambda)

    cuts = sfu.get_cuts(cuts)
    ncuts = len(cuts)
    if start is not None:
        start = sfu.get_start(start, nparams)

    bs_list = sfr.bootstrap(nsamples, dataset, pmf, cuts, start, True, parametric)
    jk_list = sfr.jackknifedata(dataset, pmf)

    nblocks = dataset.getNblocks()

    # construct the massive tuple of return values
    samples = np.zeros((nsamples, nblocks), dtype=np.int32)
    estimates = np.zeros((nsamples, nparams))
    deviance = np.zeros((nsamples))
    thres = np.zeros((nsamples, ncuts))
    slope = np.zeros((nsamples, ncuts))
    Rpd = np.zeros((nsamples))
    Rkd = np.zeros((nsamples))
    for row_index in xrange(nsamples):
        samples[row_index] = bs_list.getData(row_index)
        estimates[row_index] = bs_list.getEst(row_index)
        deviance[row_index] = bs_list.getdeviance(row_index)
        thres[row_index] = [bs_list.getThres_byPos(row_index, j) for j in xrange(ncuts)]
        slope[row_index] = [bs_list.getSlope_byPos(row_index, j) for j in xrange(ncuts)]
        Rpd[row_index] = bs_list.getRpd(row_index)
        Rkd[row_index] = bs_list.getRkd(row_index)

    thacc = np.zeros((ncuts))
    thbias = np.zeros((ncuts))
    slacc = np.zeros((ncuts))
    slbias = np.zeros((ncuts))
    for cut in xrange(ncuts):
        thacc[cut] = bs_list.getAcc_t(cut)
        thbias[cut] = bs_list.getBias_t(cut)
        slacc[cut] = bs_list.getAcc_t(cut)
        slbias[cut] = bs_list.getBias_t(cut)

    ci_lower = sfr.vector_double(nparams)
    ci_upper = sfr.vector_double(nparams)

    for param in xrange(nparams):
        ci_lower[param] = bs_list.getPercentile(0.025, param)
        ci_upper[param] = bs_list.getPercentile(0.975, param)

    outliers = np.zeros((nblocks), dtype=np.bool)
    influential = np.zeros((nblocks))

    for block in xrange(nblocks):
        outliers[block] = jk_list.outlier(block)
        influential[block] = jk_list.influential(block, ci_lower, ci_upper)

    return samples, estimates, deviance, thres, thbias, thacc, slope, slbias, slacc, Rpd, Rkd, outliers, influential

def mcmc( data, start=None, nsamples=10000, nafc=2, sigmoid='logistic',
        core='mw0.1', priors=None, stepwidths=None, sampler="MetropolisHastings", gammaislambda=False):
    """ Markov Chain Monte Carlo sampling for a psychometric function.

    Parameters
    ----------

    data : A list of lists or an array of data.
        The first column should be stimulus intensity, the second column should
        be number of correct responses (in 2AFC) or number of yes- responses (in
        Yes/No), the third column should be number of trials. See also: the examples
        section below.

    start : sequence of floats of length number of model parameters
        Starting values for the markov chain. If this is None, the MAP estimate
        will be used.

    nsamples : int
        Number of samples to be taken from the posterior (note that due to
        suboptimal sampling, this number may be much lower than the effective
        number of samples.

    nafc : int
        Number of responses alternatives for nAFC tasks. If nafc==1 a Yes/No task is
        assumed.

    sigmoid : string
        Name of the sigmoid to be fitted. Valid sigmoids include:
                logistic
                gauss
                gumbel_l
                gumbel_r
        See `swignifit.utility.available_sigmoids()` for all available sigmoids.

    core : string
        \"core\"-type of the psychometric function. Valid choices include:
                ab       (x-a)/b
                mw%g     midpoint and width
                linear   a+bx
                log      a+b log(x)
        See `swignifit.utility.available_cores()` for all available sigmoids.

    priors : sequence of strings length number of parameters
        Prior distributions on the parameters of the psychometric function.
        These are expressed in the form of a list of prior names.
        Valid prior choices include:
                Uniform(%g,%g)
                Gauss(%g,%g)
                Beta(%g,%g)
                Gamma(%g,%g)
                nGamma(%g,%g)
                if an invalid prior or `None` is selected, no constraints are imposed at all.
        See `swignifit.utility.available_priors()` for all available sigmoids.

        if an invalid prior is selected, no constraints are imposed on that
        parameter resulting in an improper prior distribution.

    stepwidths : sequence of floats of length number of model parameters
        Standard deviations of the proposal distribution. The best choice is
        sometimes a bit tricky here. However, as a rule of thumb we can
        state: if the stepwidths are too small, the samples might not cover
        the whole posterior, if the stepwidths are too large, most steps
        will leave the area of high posterior density and will therefore be
        rejected.  Thus, in general stepwidths should be somewhere in the
        middle.

    sampler : string
        The type of MCMC sampler to use.
        See: `sw.utility.available_samplers()` for a list of available samplers.

    gammaislambda : boolean
        Set the gamma == lambda prior.


    Output
    ------

    (estimates, deviance, posterior_predictive_data,
    posterior_predictive_deviances, posterior_predictive_Rpd,
    posterior_predictive_Rkd, logposterior_ratios, accept_rate)

    estimates : numpy array, shape: (nsamples, nparameters)
        Parameters sampled from the posterior density of parameters given the data.

    deviances : numpy array, length: nsamples
        Associated deviances for each estimate

    posterior_predictive_data : numpy array, shape: (nsamples, nblocks)
        Data that are simulated by sampling from the joint posterior of data and
        parameters. They are important for model checking.

    posterior_predictive_deviances : numpy array, length: nsamples
        The deviances that are associated with the posterior predictive data. A
        particular way of model checking could be to compare the deviances and the
        posterior predicitive deviances. For a good model these should be relatively
        similar.

    posterior_predictive_Rpd : numpy array, length: nsamples
        Correlations between psychometric function and deviance residuals
        associated with posterior predictive data

    posterior_predictive_Rkd : numpy array, length: nsamples
        Correlations between block index and deviance residuals associated with
        posterior predictive data.

    logposterior_ratios :  numpy array, shape: (nsamples, nblocks)
        Ratios between the full posetrior and the posterior for a single block
        for all samples. Used for calculating the KL-Divergence to detrmine
        influential observations in the Bayesian paradigm.

    accept_rate : float
        The number of proposed MCMC samples that were accepted.

    Example
    -------
    >>> x = [float(2*k) for k in xrange(6)]
    >>> k = [34,32,40,48,50,48]
    >>> n = [50]*6
    >>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    >>> priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')
    >>> stepwidths = (1.,1.,0.01)
    >>> (estimates, deviance, posterior_predictive_data,
         posterior_predictive_deviances, posterior_predictive_Rpd,
         posterior_predictive_Rkd, logposterior_ratios, accept_rate) \
         = mcmc(d,nsamples=10000,priors=priors,stepwidths=stepwidths)
    >>> mean(estimates[:,0])
    2.4811791550665272
    >>> mean(estimates[:,1])
    7.4935217545849184
    """

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors, gammaislambda=gammaislambda)

    if start is not None:
        start = sfu.get_start(start, nparams)
    else:
        # use mapestimate
        opt = sfr.PsiOptimizer(pmf, dataset)
        start = opt.optimize(pmf, dataset)

    proposal = sfr.GaussRandom()
    if sampler not in sfu.sampler_dict.keys():
        raise sfu.PsignifitException("The sampler: " + sampler + " is not available.")
    else:
        sampler  = sfu.sampler_dict[sampler](pmf, dataset, proposal)
    sampler.setTheta(start)

    if stepwidths != None:
        stepwidths = np.array(stepwidths)
        if len(stepwidths.shape)==2:
            if isinstance ( sampler, sfr.GenericMetropolis ):
                sampler.findOptimalStepwidth ( sfu.make_pilotsample ( stepwidths ) )
            elif isinstance ( sampler, sfr.MetropolisHastings ):
                sampler.setStepSize ( sfr.vector_double( stepwidths.std(0) ) )
            else:
                raise sfu.PsignifitException("You provided a pilot sample but the selected sampler does not support pilot samples")
        elif len(stepwidths) != nparams:
            raise sfu.PsignifitException("You specified \'"+str(len(stepwidths))+\
                    "\' stepwidth(s), but there are \'"+str(nparams)+ "\' parameters.")
        else:
            if isinstance ( sampler, sfr.DefaultMCMC ):
                for i,p in enumerate(stepwidths):
                    p = sfu.get_prior(p)
                    sampler.set_proposal(i, p)
            else:
                sampler.setStepSize(sfr.vector_double(stepwidths))

    post = sampler.sample(nsamples)

    nblocks = dataset.getNblocks()

    estimates = np.zeros((nsamples, nparams))
    deviance = np.zeros(nsamples)
    posterior_predictive_data = np.zeros((nsamples, nblocks))
    posterior_predictive_deviances = np.zeros(nsamples)
    posterior_predictive_Rpd = np.zeros(nsamples)
    posterior_predictive_Rkd = np.zeros(nsamples)
    logposterior_ratios = np.zeros((nsamples, nblocks))

    for i in xrange(nsamples):
        for j in xrange(nparams):
            estimates[i, j] = post.getEst(i, j)
        deviance[i] = post.getdeviance(i)
        for j in xrange(nblocks):
            posterior_predictive_data[i, j] = post.getppData(i, j)
            logposterior_ratios[i,j] = post.getlogratio(i,j)
        posterior_predictive_deviances[i] = post.getppDeviance(i)
        posterior_predictive_Rpd[i] = post.getppRpd(i)
        posterior_predictive_Rkd[i] = post.getppRkd(i)

    accept_rate = post.get_accept_rate()

    return (estimates, deviance, posterior_predictive_data,
        posterior_predictive_deviances, posterior_predictive_Rpd,
        posterior_predictive_Rkd, logposterior_ratios, accept_rate)

def mapestimate ( data, nafc=2, sigmoid='logistic', core='ab', priors=None,
        cuts = None, start=None, gammaislambda=False):
    """ MAP or constrained maximum likelihood estimation for a psychometric function.

    Parameters
    ----------

    data : A list of lists or an array of data.
        The first column should be stimulus intensity, the second column should
        be number of correct responses (in 2AFC) or number of yes- responses (in
        Yes/No), the third column should be number of trials. See also: the examples
        section below.

    nafc : int
        Number of responses alternatives for nAFC tasks. If nafc==1 a Yes/No task is
        assumed.

    sigmoid : string
        Name of the sigmoid to be fitted. Valid sigmoids include:
                logistic
                gauss
                gumbel_l
                gumbel_r
        See `swignifit.utility.available_sigmoids()` for all available sigmoids.

    core : string
        \"core\"-type of the psychometric function. Valid choices include:
                ab       (x-a)/b
                mw%g     midpoint and width
                linear   a+bx
                log      a+b log(x)
        See `swignifit.utility.available_cores()` for all available sigmoids.

    priors : sequence of strings length number of parameters
        Prior distributions on the parameters of the psychometric function.
        These are expressed in the form of a list of prior names.
        Valid prior choices include:
                Uniform(%g,%g)
                Gauss(%g,%g)
                Beta(%g,%g)
                Gamma(%g,%g)
                nGamma(%g,%g)
                if an invalid prior or `None` is selected, no constraints are imposed at all.
        See `swignifit.utility.available_priors()` for all available sigmoids.

        if an invalid prior is selected, no constraints are imposed on that
        parameter resulting in an improper prior distribution.

    cuts : sequence of floats
        Cuts at which thresholds should be determined.  That is if cuts =
        (.25,.5,.75), thresholds (F^{-1} ( 0.25 ), F^{-1} ( 0.5 ), F^{-1} ( 0.75
        )) are returned.  Here F^{-1} denotes the inverse of the function
        specified by sigmoid. If cuts==None, this is modified to cuts=[0.5].

    start : sequence of floats of length number of model parameters
        Values at which to start the optimization, if None the starting value is
        determined using a coarse grid search.

    Output
    ------

    estimate, fisher, thres, slope, deviance

    estimate : numpy array length nparams
        the map/cml estimate

    fisher : numpy array shape (nparams, nparams)
        the fisher matrix

    thres : numpy array length ncuts
        the model prediction at the cuts

    slope : numpy array length ncuts
        the gradient of the psychometric function at the cuts

    deviance : numpy array length 1
        the deviance for the estimate

    Example
    -------
    >>> x = [float(2*k) for k in xrange(6)]
    >>> k = [34,32,40,48,50,48]
    >>> n = [50]*6
    >>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    >>> priors = ('flat','flat','Uniform(0,0.1)')
    >>> estimate, fisher, thres, slope, deviance = mapestimate ( d, priors=priors )
    >>> estimate
    array([ 2.75180624,  1.45717745,  0.01555658])
    >>> deviance
    array(8.0713313642328242)

    """

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data,
            nafc, sigmoid, core, priors, gammaislambda=gammaislambda)

    cuts = sfu.get_cuts(cuts)

    opt = sfr.PsiOptimizer(pmf, dataset)
    estimate = opt.optimize(pmf, dataset, sfu.get_start(start, nparams) if start is not
            None else None)
    H = pmf.ddnegllikeli(estimate, dataset)
    thres = [pmf.getThres(estimate, c) for c in cuts]
    slope = [pmf.getSlope(estimate, th) for th in thres]
    deviance = pmf.deviance(estimate, dataset)

    # convert to numpy stuff
    estimate = np.array(estimate)
    fisher = np.zeros((nparams,nparams))
    for (i,j) in ((i,j) for i in xrange(nparams) for j in xrange(nparams)):
        fisher[i,j] = sfr.doublep_value(H(i,j))
    thres = np.array(thres)
    slope = np.array(slope)
    deviance = np.array(deviance)

    return estimate, fisher, thres, slope, deviance

def diagnostics(data, params, nafc=2, sigmoid='logistic', core='ab', cuts=None, gammaislambda=False):
    """ Some diagnostic statistics for a psychometric function fit.

    This function is a bit messy since it has three functions depending on the
    type of the `data` argument.

    Parameters
    ----------

    data : variable
        real data : A list of lists or an array of data.
            The first column should be stimulus intensity, the second column should
            be number of correct responses (in 2AFC) or number of yes- responses (in
            Yes/No), the third column should be number of trials. See also: the examples
            section below.
        intensities : sequence of floats
            The x-values of the psychometric function, then we obtain only the
            predicted values.
        no data : empty sequence
            In this case we evaluate the psychometric function at the cuts. All
            other return values are then irrelevant.

    params : sequence of len nparams
        parameter vector at which the diagnostic information should be evaluated

    nafc : int
        Number of responses alternatives for nAFC tasks. If nafc==1 a Yes/No task is
        assumed.

    sigmoid : string
        Name of the sigmoid to be fitted. Valid sigmoids include:
                logistic
                gauss
                gumbel_l
                gumbel_r
        See `swignifit.utility.available_sigmoids()` for all available sigmoids.

    core : string
        \"core\"-type of the psychometric function. Valid choices include:
                ab       (x-a)/b
                mw%g     midpoint and width
                linear   a+bx
                log      a+b log(x)
        See `swignifit.utility.available_cores()` for all available sigmoids.

    cuts : sequence of floats
        Cuts at which thresholds should be determined.  That is if cuts =
        (.25,.5,.75), thresholds (F^{-1} ( 0.25 ), F^{-1} ( 0.5 ), F^{-1} ( 0.75
        )) are returned.  Here F^{-1} denotes the inverse of the function
        specified by sigmoid. If cuts==None, this is modified to cuts=[0.5].

    Output
    ------

    (predicted, deviance_residuals, deviance, thres, Rpd, Rkd)

    predicted : numpy array of length nblocks
        predicted values associated with the respective stimulus intensities

    deviance_residuals : numpy array of length nblocks
        deviance residuals of the data

    deviance float
        deviance of the data

    thres : numpy array length ncuts
        the model prediction at the cuts

    Rpd : float
        correlation between predicted performance and deviance residuals

    Rkd : float
        correlation between block index and deviance residuals

    Example
    -------
    >>> x = [float(2*k) for k in xrange(6)]
    >>> k = [34,32,40,48,50,48]
    >>> n = [50]*6
    >>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    >>> prm = [2.75, 1.45, 0.015]
    >>> pred,di,D,thres,slope,Rpd,Rkd = diagnostics(d,prm)
    >>> D
    8.0748485860836254
    >>> di[0]
    1.6893279652591433
    >>> Rpd
    -0.19344675783032755

    """

    # here we need to hack stuff, since data can be either 'real' data, or just
    # a list of intensities, or just an empty sequence

    # in order to remain compatible with psipy we must check for an empty
    # sequence here, and return a specially crafted return value in that case.
    # sorry..
    # TODO after removal of psipy we can probably change this.
    if op.isSequenceType(data) and len(data) == 0:
        pmf, nparams =  sfu.make_pmf(sfr.PsiData([0],[0],[0],1), nafc, sigmoid, core, None, gammaislambda=gammaislambda )
        thres = np.array([pmf.getThres(params, cut) for cut in sfu.get_cuts(cuts)])
        slope = np.array([pmf.getSlope(params, th ) for th in thres])
        return np.array([]), np.array([]), 0.0, thres, np.nan, np.nan

    shape = np.shape(np.array(data))
    intensities_only = False
    if len(shape) == 1:
        # just intensities, make a dataset with k and n all zero
        k = n = [0] * shape[0]
        data  = [[xx,kk,nn] for xx,kk,nn in zip(data,k,n)]
        intensities_only = True
    else:
        # data is 'real', just do nothing
        pass

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, None, gammaislambda=gammaislambda)
    cuts = sfu.get_cuts(cuts)
    params = sfu.get_params(params, nparams)
    predicted = np.array([pmf.evaluate(intensity, params) for intensity in
            dataset.getIntensities()])

    if intensities_only:
        return predicted
    else:
        deviance_residuals = pmf.getDevianceResiduals(params, dataset)
        deviance = pmf.deviance(params, dataset)
        thres = np.array([pmf.getThres(params, cut) for cut in cuts])
        slope = np.array([pmf.getSlope(params, th ) for th in thres])
        rpd = pmf.getRpd(deviance_residuals, params, dataset)
        rkd = pmf.getRkd(deviance_residuals, dataset)
        return predicted, deviance_residuals, deviance, thres, slope, rpd, rkd

def asir ( data, nsamples=2000, nafc=2, sigmoid="logistic",
        core="mw0.1", priors=None, gammaislambda=False, propose=25 ):
    dataset, pmf, nparams = sfu.make_dataset_and_pmf ( data, nafc, sigmoid, core, priors, gammaislambda=gammaislambda )

    posterior = sfr.independent_marginals ( pmf, dataset )
    if nsamples > 0:
        samples   = sfr.sample_posterior ( pmf, dataset, posterior, nsamples, propose )
        sfr.sample_diagnostics ( pmf, dataset, samples )

        out = {'mcestimates': np.array( [ [samples.getEst ( i, par ) for par in xrange ( nparams ) ] for i in xrange ( nsamples )]),
            'mcdeviance': np.array( [ samples.getdeviance ( i ) for i in xrange ( nsamples ) ] ),
            'mcRpd':                    np.array ( [ samples.getRpd ( i ) for i in xrange ( nsamples ) ] ),
            'mcRkd':                    np.array ( [ samples.getRkd ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_data': np.array ( [ samples.getppData ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_deviance': np.array ( [ samples.getppDeviance ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_Rpd': np.array ( [ samples.getppRpd ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_Rkd': np.array ( [ samples.getppRkd ( i ) for i in xrange ( nsamples ) ] ),
            'logposterior_ratios':      np.array ( [ [samples.getlogratio ( i,j ) for j in xrange(len(data)) ] for i in xrange ( nsamples ) ] ),
            'duplicates':               samples.get_accept_rate (),
            'posterior_approximations_py': [posterior.get_posterior(i) for i in xrange ( nparams ) ],
            'posterior_approximations_str': [r"$\mathcal{N}(%.2f,%.2f)$" % (posterior.get_posterior(0).getprm(0),posterior.get_posterior(0).getprm(1)),
                r"$\mathrm{Gamma}(%.2f,%.2f)$" % (posterior.get_posterior(1).getprm(0),posterior.get_posterior(1).getprm(1)),
                r"$\mathrm{Beta}(%.2f,%.2f)$" % (posterior.get_posterior(2).getprm(0),posterior.get_posterior(2).getprm(1))],
            'posterior_grids':          [ posterior.get_grid ( i ) for i in xrange ( nparams ) ],
            'posterior_margin':         [ posterior.get_margin ( i ) for i in xrange ( nparams ) ],
            'resampling-entropy':       samples.get_entropy ()
            }

    else:

        out = {'posterior_approximations_py': [posterior.get_posterior(i) for i in xrange ( nparams ) ],
            'posterior_approximations_str': [r"$\mathcal{N}(%.2f,%.2f)$" % (posterior.get_posterior(0).getprm(0),posterior.get_posterior(0).getprm(1)),
            r"$\mathrm{Gamma}(%.2f,%.2f)$" % (posterior.get_posterior(1).getprm(0),posterior.get_posterior(1).getprm(1)),
            r"$\mathrm{Beta}(%.2f,%.2f)$" % (posterior.get_posterior(2).getprm(0),posterior.get_posterior(2).getprm(1))],
            'posterior_grids':          [ posterior.get_grid ( i ) for i in xrange ( nparams ) ],
            'posterior_margin':         [ posterior.get_margin ( i ) for i in xrange ( nparams ) ] }

    if nparams==4:
        out['posterior_approximations_str'].append ( r"$\mathrm{Beta}(%.2f,%.2f)$" % (posterior.get_posterior(3).getprm(0),posterior.get_posterior(3).getprm(1)) )

    return out
