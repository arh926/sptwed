% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saddle_likelihood.R
\name{saddle_likelihood}
\alias{saddle_likelihood}
\title{Saddlepoint approximation to the Tweedie Compound Poisson-gamma likelihood}
\usage{
saddle_likelihood(y = NULL, p = NULL, mu = NULL, phi = NULL, epsilon = 1/6)
}
\arguments{
\item{y}{response}

\item{p}{index parameter}

\item{mu}{mean vector}

\item{phi}{dispersion vector}
}
\description{
Saddlepoint approximation to the Tweedie Compound Poisson-gamma likelihood
}
\examples{
\dontrun{
require(tweedie)
N <- 10
mu <- rnorm(N, 10, 2)
phi <- rgamma(N, 2, 2)
y <- rtweedie(N, 1.5, mu, phi)
saddle_likelihood(y=y,p=1.5,mu=mu,phi=phi)
}
}
\keyword{saddle_likelihood}
