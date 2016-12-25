#/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

"""Variety of utilities for working with swig generated code."""

__docformat__ = "restructuredtext"

import operator as op
import re
import numpy as np
import swignifit_raw as sfr

class PsignifitException(Exception):
    pass

def extract_subclasses(base, sub_func):
    """Recursively extract subclasses, given a swig base class.

    Parameters
    ----------
    base : swig class
        The base class from which to start.
    sub_func : string
        The function or attribute to use as name for subclass.

    Returns
    -------
    subs : dict
        A dictionary mapping subclass names to constructors.

    """
    to_visit = base.__subclasses__()
    subclasses = dict()
    for cl in to_visit:
        descriptor = eval("cl."+sub_func)
        if descriptor not in subclasses.keys():
            subclasses[descriptor] = cl
            to_visit.extend(cl.__subclasses__())
    return subclasses

def extract_subclasses_descriptor(base):
    """Recursively extract subclasses, using the `getDescriptor()` method."""
    return extract_subclasses(base, "getDescriptor()")

def extract_subclasses_reflection(base):
    """Recursively extract subclasses, using the `__name__` attribute."""
    return extract_subclasses(base, "__name__")

sig_dict = extract_subclasses_descriptor(sfr.PsiSigmoid)
core_dict = extract_subclasses_descriptor(sfr.PsiCore)
prior_dict = extract_subclasses_reflection(sfr.PsiPrior)
sampler_dict = extract_subclasses_reflection(sfr.PsiSampler)

def available_cores():
    print "The following cores are availabe:"
    print core_dict.keys()

def available_sigmoids():
    print "The following sigmoids are available:"
    print sig_dict.keys()

def available_priors():
    print "The following priors are available:"
    print prior_dict.keys()

def available_samplers():
    print "The following mcmc samplers are available:"
    print sampler_dict.keys()

def make_dataset(data, nafc):
    """Create a PsiData object from column based input.

    Parameters
    ----------
    data : sequence on length 3 sequences
        Psychometric data in colum based input,
        e.g.[[1, 1, 5], [2, 3, 5] [3, 5, 5]].
    nafc : int
        Number of alternative choices in forced choice procedure.

    Returns
    -------
    data: PsiData
        Dataset object.

    """
    data = np.array(data).T
    x = sfr.vector_double(map(float, data[0]))
    k = sfr.vector_int(map(int, data[1]))
    N = sfr.vector_int(map(int, data[2]))
    return sfr.PsiData(x,N,k,nafc)

def make_pmf(dataset, nafc, sigmoid, core, priors, gammaislambda=False):
    """Assemble PsiPsychometric object from model parameters.

    Parameters
    ----------
    dataset: sequence of length 3 sequences
        Psychometric data in colum based input,
        e.g.[[1, 1, 5], [2, 3, 5] [3, 5, 5]].
    nafc : int
        Number of alternative choices in forced choice procedure.
    sigmoid : string
        Description of model sigmoid.
    core : string
        Description of model core.
    priors : sequence of strings
        The model priors.
    gammaislambda : bool
        Constrain guessing rate and lapsing rate to be equal

    Returns
    -------
    pmf : PsiPsychometric
        Model object.
    nparams : int
        Number of free parameters in model..

    """
    sigmoid = get_sigmoid(sigmoid)
    core = get_core(core, dataset, sigmoid)
    pmf = sfr.PsiPsychometric(nafc, core, sigmoid)
    if gammaislambda:
        pmf.setgammatolambda()
    nparams = pmf.getNparams()
    set_priors(pmf,priors)
    return pmf, nparams

def make_dataset_and_pmf(data, nafc, sigmoid, core, priors, gammaislambda=False ):
    """Assemble PsiData and PsiPsychometric objects simultaneously.

    Parameters
    ----------
    see: make_dataset and make_pmf

    Returns
    -------
    data : PsiData
        Dataset object.
    pmf : PsiPsychometric
        Model object.
    nparams : int
        Number of free parameters.

    """
    dataset = make_dataset(data, nafc)
    pmf, nparams = make_pmf(dataset, nafc, sigmoid, core, priors, gammaislambda)
    return dataset, pmf, nparams

def get_sigmoid(descriptor):
    """Convert string representation of sigmoid to PsiSigmoid object.

    Parameters
    ----------
    descriptor : string
        Description of sigmoid.

    Returns
    -------
    sigmoid : subclass of PsiSigmoid
        An instantiated sigmoid of the requested type.

    Raises
    ------
    ValueError
        If the requested sigmoid is not available.

    See Also
    --------
    `available_sigmoids()`

    """
    if not sig_dict.has_key(descriptor):
        raise ValueError("The sigmoid '%s' you requested, is not available." %
                descriptor)
    return sig_dict[descriptor]()

def get_core(descriptor, data, sigmoid):
    """Convert string representation of core to PsiCore object.

    Parameters
    ----------
    descriptor : string
        Description of core.
    data : PsiData
        Instantiated dataset.
    sigmoid : PsiSigmoid
        Instantiated sigmoid.

    Returns
    -------
    prior : subclass of PsiCore
        An instantiated core of the requested type.

    Raises
    ------
    ValueError
        If the requested core is not available.

    See Also
    --------
    `available_cores()`

    Notes
    -----
    The core objects may require a dataset and a sigmoid type identifier to be
    instantiated. See the Psi++ code for details.

    """
    descriptor, parameter = re.match('([a-z]+)([\d\.]*)', descriptor).groups()
    sigmoid_type = sigmoid.getcode()
    if descriptor not in core_dict.keys():
        raise ValueError("The core '%s' you requested, is not available." %
                descriptor)
    if len(parameter) > 0:
        return core_dict[descriptor](data, sigmoid_type, float(parameter))
    else:
        return core_dict[descriptor](data, sigmoid_type)

def get_prior(prior):
    """Convert string based representation of prior to PsiPrior object.

    Parameters
    ----------
    prior : string
        Description of prior, with paramters.

    Returns
    -------
    prior : PsiPrior
        An instantiated prior of the requested type.

    See Also
    --------
    `available_priors()`

    Notes
    -----
    This function does not raise any error and silently returns `None` if the
    prior does not exist.

    """
    try:
        prior = "sfr."+"Prior(".join(prior.split('('))
        return eval(prior)
    except Exception, e:
        return None

def set_priors(pmf, priors):
    """Set the priors to be used in the model object.

    Parameters
    ----------
    pmf : PsiPsychometric object
        Instantiated model.
    priors : Sequence of strings of length of free parameters of `pmf`
         Model priors.

    Raises
    ------
    ValueError
        If the number of priors is not equal to the number of free parameters in
        the model.

    """
    if priors is not None:
        nparams = pmf.getNparams()
        if len(priors) != nparams:
            raise ValueError("You specified %d priors, " % len(priors) +
                    "but there are %d free parameters in the model." % nparams)
        for (i,p) in enumerate((get_prior(p) for p in priors)):
            if p is not None:
                pmf.setPrior(i, p)

def get_start(start, nparams):
    """Convert sequence of starting values to vector_double type.

    Parameters
    ----------
    start : sequence of numbers
        Starting values.
    nparams : int
        Number of free parameters of the model.

    Returns
    -------
    start: vector_double
        Starting values.

    Raises
    ------
    ValueError
        If the length of the sequence is not equal to the number of free
        parameters.

    """
    if len(start) != nparams:
            raise ValueError("You specified %d starting value(s), " % len(start)
                    +"but there are %d parameters." % nparams)
    else:
        return sfr.vector_double(start)

def get_params(params, nparams):
    """Convert sequence of parameter values to vector_double type.

    Parameters
    ----------
    params : sequence of numbers
        Parameter values.
    nparams : int
        Number of free parameters in the model.

    Returns
    -------
    params : vector_double
        Parameter values.

    Raises
    ------
    ValueError
        If the length of the sequence is not equal to the number of free
        parameters.

    """
    if len(params) != nparams:
                raise ValueError("You specified %d parameters, " % len(params) +
                        "but the model has parameters." % nparams)
    else:
        return sfr.vector_double(params)

def get_cuts(cuts):
    """ Convert `cuts` argument to vector_double type.

    Argument can be None, a number or a sequence of numbers. If None, there is
    only one cut at 0.5. If `cuts` is a number, function returns a vector_double
    with that number as a single element. If its a sequence, that sequence will
    be converted to vector_double type.

    Parameters
    ----------
    cuts : None, number or sequence of numbers
        Cut values

    Returns
    -------
    cuts : vector_double
        Cut values.

    Raises
    ------
    TypeError
        If `cuts` is not None, a number or a sequence of numbers.

    """
    if cuts is None:
        return sfr.vector_double([0.5])
    elif op.isSequenceType(cuts) and np.array([op.isNumberType(a) for a in cuts]).all():
        return sfr.vector_double(cuts)
    elif op.isNumberType(cuts):
        return sfr.vector_double([cuts])
    else:
        raise TypeError("'cuts' must be either None, a number or a "+
                "sequence of numbers.")

def make_pilotsample ( mcsamples ):
    """create an MCList from a set of pilot samples

    Deviances of the pilot list will be meaningless!

    Parameters
    ----------
    mcsamples : array of Nsamples x Nparams
        pilot samples

    Returns
    -------
    pilot : PsiMClist
        MClist with the pilot samples
    """
    N,nprm = mcsamples.shape
    pilot = sfr.PsiMClist ( N, nprm )
    for i in xrange ( N ):
        pilot.setEst ( i, mcsamples[i,:], -1 )
    return pilot
