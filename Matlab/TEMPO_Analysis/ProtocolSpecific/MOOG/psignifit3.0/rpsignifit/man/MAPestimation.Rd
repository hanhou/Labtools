\name{MAPestimation}
\alias{MAPestimation}
\title{Maximum a Posteriori estimation}
\description{This function performs maximum a Posteriori point estimation of a psychometric function.
    Maximum a Posteriori estimation of psychometric functions can be considered a generalization of
    constrained maximum likelihood as used for the bootstrap procedures.
}
\usage{MAPestimation ( psignidata )}
\arguments{
\item{psignidata}{A PsigniSetup that contains the data to be fitted. Depending on the type of task
    there will be different numbers of parameters that are fitted. In nAFC tasks, there are three parameters.
    The first two determine the shape of the sigmoid and the third determines the upper asymptote. In Yes/No tasks,
    there is a fourth parameter that also determines the lower asymptote.}
}
\references{
    Wichmann, FA & Hill, NJ ( 2001 ): The psychometric function: I. Fitting, sampling, and goodness of fit. Perception & Psychophysics, 63(8), 1293-1313.
}
\examples{
    x <- seq(0,10,2)
    k <- c(34,32,40,48,50,48)
    n <- rep(50,6)
    D <- PsigniSetup ( x, k, n )
    fit <- MAPestimation ( D )
    \dontrun{print ( fit )}
}


