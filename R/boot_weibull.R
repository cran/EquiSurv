################################################################################
#' @title Parametric Bootstrap of time-to-event data following a Weibull distribution
#'
#' @description Function generating bootstrap data according to a Weibull distribution (specified by a model parameter \eqn{\theta}),
#' assuming exponentially distributed right-censoring (specified by a rate C). After data generation again a model is fitted and evaluated
#' at a pre-specified time point \eqn{t_0} yielding the response vector.
#'
#' @name boot_weibull
#' @export
#' @import survival
#' @param t0 time point of interest
#' @param B number of bootstrap repetitions. The default is B=1000
#' @param theta parameter of the Weibull distribution, theta=(shape,scale)
#' @param C rate of the exponential distribution specifiying the censoring
#' @param N size of the dataset = number of observations
#' @return A vector of length B containing the estimated survival at t0
#' @examples
#' t0<-2
#' N<-30
#' C<-1
#' boot_weibull(t0=t0,theta=c(1,3),C=C,N=N)
################################################################################

boot_weibull <- function(t0,B=1000,theta,C,N){

  boot <- numeric()

  #data simulating function
  simulWeib <- function(N, shape=shape, scale=scale, rateC)
  {
    # survival time assuming a Weibull distribution
    Tlat <- rweibull(N, shape = shape, scale = scale)

    # censoring times (exponential distributed)
    C <- rexp(n=N, rate=rateC)

    # follow-up times and event indicators
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C) #uncensored (experiences an event)

    # data set
    data.frame(id=1:N,
               time=time,
               status=status
    )
  }

  #bootstrap to obtain the standard error
  for(l in 1:B){
    data_rb <- simulWeib(N=N, shape=theta[1], scale=theta[2], rateC=C)
    mod1b <- survreg(Surv(time,status)~1,data=data_rb,dist="weibull")
    theta_rb <- c(1/mod1b$scale,exp(mod1b$coefficients))
    F <- function(theta,x) {
      pweibull(x,shape=theta[1],scale=theta[2])
    }
    surv_rb <- function(t){
      1-F(theta_rb,t)
    }

    boot[l] <- surv_rb(t0)
  }
  return(boot)
}
