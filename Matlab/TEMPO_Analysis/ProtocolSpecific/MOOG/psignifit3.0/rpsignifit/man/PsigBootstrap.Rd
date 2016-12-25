\name{PsigBootstrap}
\alias{PsigBootstrap}
\title{Bootstrapping psychometric functions}
\description{This function performs bootstrap inference of a psychometric function. The function also determines the bias correction and
    acceleration constants that are needed to adjust the bootstrap percentiles.}
\usage{PsigBootstrap ( psignidata, number.of.samples=2000, generating=-999 )}
\arguments{
\item{psignidata}{A PsigniSetup that contains the data to be fitted. Depending on the type of task
    there will be different numbers of parameters that are fitted. In nAFC tasks, there are three parameters.
    The first two determine the shape of the sigmoid and the third determines the upper asymptote. In Yes/No tasks,
    there is a fourth parameter that also determines the lower asymptote.}
\item{number.of.samples}{Number of monte carlo samples to be generated}
\item{generating}{The generating parameters for a parametric bootstrap. Typically, the map estimate will be used here
    for parametric bootstrap. If generating is set to -999 a nonparametric bootstrap will be performed. Parametric bootstrap
    around the MAP estimate is performed by Wichmann & Hill (2001). If a nonparametric bootstrap is requested, there will
    be an additional parametric bootstrap run around the map estimate to determine the values that are needed for the
    goodness of fit assessment.}
}
\references{
    Wichmann, FA & Hill, NJ ( 2001 ): The psychometric function: I. Fitting, sampling, and goodness of fit. Perception & Psychophysics, 63(8), 1293-1313.
    Wichmann, FA & Hill, NJ ( 2001 ): The psychometric function: II. Bootstrap-based confidence intervals and sampling. Perception & Psychophysics, 63(8), 1314-1329.
}
\examples{
    x <- seq(0,10,2)
    k <- c(34,32,40,48,50,48)
    n <- rep(50,6)
    D <- PsigniSetup ( x, k, n )
    boots <- PsigBootstrap ( D )
    \dontrun{print ( boots )}
}


