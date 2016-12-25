\name{PsigniSetup}
\alias{PsigniSetup}
\title{Data setup for fitting psychometric functions}
\description{
    In order to fit psychometric functions, there are a number of that need to be fixed.
    This function returns a list that contains everything that needs to be known to fit a psychometric
    function.
}
\usage{
    PsigniSetup ( x, k, n, priors=list("","","Uniform(0,.1)"), sigmoid="logistic", core="mw0.1", number.of.alternatives=2, cuts=0.5 )
}
\arguments{
\item{x}{Stimulus intensities at which stimuli were presented. This should be a sequence in the order in which the stimuli were presented
    during the experiment.}
\item{k}{Number of correct responses (in nAFC trials) or number of YES responses (in Yes/No trials). These numbers should correspond to the
    entries in x.}
\item{n}{Number of repetitions of the respective stimuli.}
\item{priors}{a list of priors for the respective parameters. Psychometric functions for nAFC trials are governed by three parameters, the first two
    determine the shape of the sigmoid function and the last determines the upper asympote. Psychometric functions for Yes/No trials are described by
    four parameters. These are the same as for the nAFC case and an additional parameter for the lower asympote.
    Valid priors are "Gauss(mu,sigma)", "Gamma(shape,rate)", "Uniform(low,high)", "Beta(alpha,beta)". Invalid priors are considered flat, i.e. no
    restriction on the parameter. That means, the following priors are equivalent: "unconstrained", "flat", "", "improper". Note however, that also
    "gauss(mu,sigma)" and "gamma shape rate" are invalid and thus pose no restriction on the respective parameter!}
\item{sigmoid}{Shape of the sigmoidal nonlinearity. Valid choices are "logistic", "exp", "gumbel_l", "gumbel_r", "gauss", "cauchy". These
    correspond to the cumulative distribution functions of the respective distributions.}
\item{core}{By combining a sigmoid and a core different parameterizations of the psychometric function can be created. Valid cores are
    "ab", "mw%g" with "%g" a number larger than 0 and less than 0.5, "linear", "log", "poly", and "weibull". "ab" corresponds to a core of the form
    (x-alpha)/beta. If this is combined with logistic or gauss sigmoids it corresponds to the respective functions in psignifit 1.
    "mw%g" corresponds to a parameterization in terms of midpoint and width of the rising part of the sigmoid. This width is defined as the length of the
    interval on which the sigmoidal part reaches from "%g" to 1-"%g". The "linear" core corresponds to a+bx. The "log" core is a+b*log(x).
    If this is combined with a "gumbel_l" sigmoid, the resulting cdf is a weibull. The "weibull" core also fits parameters on a log-scale, thus resulting
    in a weibull cdf if combined with the "gumbel_l" sigmoid. In this case, the two parameters are midpoint and slope at the midpoint. Finally, the
    "poly" core corresponds to a linear core of the form (x/a)**b. In combination with an "exp" sigmoid, this results in the classical parameterization
    of the weibull cdf. Typically this also gives the best convergence of the weibull cdf for MCMC and for MAP estimation.}
\item{number.of.alternatives}{Is the number of stimuli presented in an nAFC task. If number.of.alternatives<2, a Yes/No task is assumed and both lower
    and upper asymptote are fitted.}
\item{cuts}{A sequence of cuts at which the thresholds should be determined. These thresholds are given as the inverse of the sigmoid applied to the
    inverse of the core. In a 2AFC task with no lapsing, the cuts correspond to the stimulus levels that yield (0.5+0.5*cut)*100\% performance.}
}
\references{
    Wichmann, FA & Hill, NJ ( 2001 ): The psychometric function: I. Fitting, sampling, and goodness of fit. Perception & Psychophysics, 63(8), 1293-1313.
    Wichmann, FA & Hill, NJ ( 2001 ): The psychometric function: II. Bootstrap-based confidence intervals and sampling. Perception & Psychophysics, 63(8), 1314-1329.
    Kuss, M, Jäkel, F & Wichmann, FA ( 2005 ): Bayesian inference for psychometric functions. Journal of Vision, 5, 478-492.
}
\examples{
    x <- seq(0,10,2)            # Stimulus intensities
    k <- c(34,32,40,48,50,48)   # Numbers of correct responses
    n <- rep(50,6)              # Numbers of stimuli presented
    D <- PsigniSetup ( x, k, n, list("Gauss(0,100)","Gauss(0,100)","Beta(2,50)"), cuts=seq(.25,.75,.25) )
    \dontrun{print (D)}
}
