% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_gaussian.R
\name{boot_gaussian}
\alias{boot_gaussian}
\title{Parametric Bootstrap of time-to-event data following a gaussian distribution}
\usage{
boot_gaussian(t0, B = 1000, theta, C, N)
}
\arguments{
\item{t0}{time point of interest}

\item{B}{number of bootstrap repetitions. The default is B=1000}

\item{theta}{parameter of the gaussian distribution, theta=(mean,sd)}

\item{C}{rate of the exponential distribution specifiying the censoring}

\item{N}{size of the dataset = number of observations}
}
\value{
A vector of length B containing the estimated survival at t0
}
\description{
Function generating bootstrap data according to a gaussian distribution (specified by a model parameter \eqn{\theta}),
assuming exponentially distributed right-censoring (specified by a rate C). After data generation again a model is fitted and evaluated
at a pre-specified time point \eqn{t_0} yielding the response vector.
}
\examples{
t0<-2
N<-30
C<-1
boot_gaussian(t0=t0,theta=c(1.7,1),C=C,N=N)
}
