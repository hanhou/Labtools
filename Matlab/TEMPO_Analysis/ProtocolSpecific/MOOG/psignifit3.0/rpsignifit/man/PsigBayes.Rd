\name{PsigBayes}
\alias{PsigBayes}
\title{Bayesian inference for psychometric functions}
\description{This function performs markov-chain monte carlo (MCMC) to generate samples form the posterior distribution
    of parameters (and also data) of the psychometric function. A number of different parameters is calculated on the fly.
    A difficult point about MCMC is the fact that although the samples will eventually be from the target distribution, this
    does not necessarily hold for the finite number of samples drawn during a simulation run. There are at least two packages
    that provide a mass of convergence diagnostics in R. These are the coda and the boa packages. Use these to avoid wrong and/or
    misleading results! Previous MCMC approaches to fitting psychometric functions (Kuss et al, 2005) used Hybrid MCMC for sampling
    the posterior. However, in our experience, a classical Metropolis-Hastings sampler performed best. That is why, the implementation
    is based on a Metropolis-Hastings sampler.}
\usage{PsigBayes ( psignidata, number.of.samples=2000, start=NULL, proposal=NULL )}
\arguments{
\item{psignidata}{A PsigniSetup that contains the data to be fitted. Depending on the type of task
    there will be different numbers of parameters that are fitted. In nAFC tasks, there are three parameters.
    The first two determine the shape of the sigmoid and the third determines the upper asymptote. In Yes/No tasks,
    there is a fourth parameter that also determines the lower asymptote.}
\item{number.of.samples}{Total number of MCMC samples to be generated. Of these, the first half is discarded and inference is based
    on the second half only.}
\item{start}{Starting values for the sampler. By default, the MAP estimate is taken. If more than one chain is to be drawn (e.g.
    for convergence diagnostics) it might be useful to start these at different points in parameter space.}
\item{proposal}{A vector of standard deviances for the (gaussian) proposal distributions. By default, these are set to c(.4,.4,.01) for
    nAFC and to c(.4,.4,.01,.01) for Yes/No. This will in many cases be a bad choice. It might be helpful to base these on standard deviations
    from pilot sample runs.}
}
\references{
    Kuss, M, Jäkel, F & Wichmann, FA ( 2005 ): Bayesian inference for psychometric functions. Journal of Vision, 5, 478-492.
}
\examples{
    x <- seq(0,10,2)            # Stimulus intensities
    k <- c(34,32,40,48,50,48)   # Numbers of correct responses
    n <- rep(50,6)              # Numbers of stimuli presented
    D <- PsigniSetup ( x, k, n, list("Gauss(0,100)","Gauss(0,100)","Beta(2,50)"), cuts=seq(.25,.75,.25) )
    mcmc <- PsigBayes ( D )
    \dontrun{print ( mcmc )}
}
