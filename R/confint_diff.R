################################################################################
#' @title Lower and upper confidence bounds for the difference of two parametric survival curves
#'
#' @description Function fitting parametric survival curves \eqn{S_1}, \eqn{S_2} to two groups and
#' yielding lower and upper (1-\eqn{\alpha})-confidence bounds for the difference \eqn{S_1-S_2} of these
#' two curves at a specific time point, based on approximating the variance via bootstrap.
#' For the bootstrap exponentially distributed random censoring is assumed and the parameters estimated from the datasets.
#' \eqn{m_1} and \eqn{m_2} are parametric survival models following a Weibull, exponential, gaussian, logistic, log-normal or log-logistic distribution.
#' For the generation of the bootstrap data exponentially distributed right-censoring is assumed and the rates estimated from the datasets.
#' See Moellenhoff and Tresch <arXiv:2009.06699> for details.
#'
#' @name confint_diff
#' @export
#' @import survival eha
#' @importFrom graphics curve legend grid par
#' @importFrom stats optimize dexp dlnorm dlogis dnorm dweibull pexp plnorm plogis pnorm pweibull qnorm rexp rlnorm rlogis rnorm rweibull var
#' @param alpha confidence level
#' @param t0 time point of interest
#' @param m1,m2 type of parametric model. Possible model types are "weibull", "exponential", "gaussian", "logistic", "lognormal" and "loglogistic"
#' @param B number of bootstrap repetitions. The default is B=1000
#' @param data_r,data_t datasets containing time and status for each individual (have to be referenced as this)
#' @param plot if TRUE, a plot of the two survival curves will be given
#' @return A list containing the difference \eqn{S_1(t_0)-S_2(t_0)}, the lower and upper (1-\eqn{\alpha})-confidence bounds and a summary of the two model fits. Further a plot of the curves is given.
#' @examples
#' data(veteran)
#' veteran_r <- veteran[veteran$trt==1,]
#' veteran_t <- veteran[veteran$trt==2,]
#' alpha<-0.05
#' t0<-80
#' confint_diff(alpha=alpha,t0=t0,m1="weibull",m2="weibull",data_r=veteran_r,data_t=veteran_t)
#' @references K.Moellenhoff and A.Tresch: Survival analysis under non-proportional hazards: investigating non-inferiority or equivalence in time-to-event data <arXiv:2009.06699>

################################################################################

confint_diff <- function(alpha,t0,m1,m2,B=1000,data_r,data_t,plot=TRUE){

  #small B
  if (B<200) {
    warning("Warning: A larger B should be choosen for a higher accuracy.")
  }

  mod1 <- survreg(Surv(time,status)~1,data=data_r,dist=m1)
  mod2 <- survreg(Surv(time,status)~1,data=data_t,dist=m2)
  N1 <- summary(mod1)$n
  N2 <- summary(mod2)$n

  g <- function(psi,v) {
    dexp(v,rate=psi)
  }
  G <- function(psi,v) {
    pexp(v,rate=psi)
  }

  likelihood <- function(data){
    f<- function(psi){-sum(log(
      g(psi,data$time)^(1-data$status)*(1-G(psi,data$time))^(data$status)
    ))
    }
    return(f)
  }

  C_r <- optimize(f=likelihood(data_r),interval=c(0,5/max(data_r$time)))$minimum
  C_t <- optimize(f=likelihood(data_t),interval=c(0,5/max(data_t$time)))$minimum

  if(m1=="weibull"){
    theta_r <- c(1/mod1$scale,exp(mod1$coefficients))
    f1 <- function(theta,x) {
      dweibull(x,shape=theta[1],scale=theta[2])
    }
    F1 <- function(theta,x) {
      pweibull(x,shape=theta[1],scale=theta[2])
    }
    v1 <- boot_weibull(B=B,t0=t0,theta=theta_r,C=C_r,N=N1)
  } else if (m1=="exponential"){
    theta_r <- 1/exp(mod1$coefficients)
    f1 <- function(theta,x) {
      dexp(x,rate=theta)
    }
    F1 <- function(theta,x) {
      pexp(x,rate=theta)
    }
    v1 <- boot_exponential(B=B,t0=t0,theta=theta_r,C=C_r,N=N1)
  } else if (m1=="gaussian"){
    theta_r <- c(mod1$coefficients,mod1$scale)
    f1 <- function(theta,x) {
      dnorm(x,mean=theta[1],sd=theta[2])
    }
    F1 <- function(theta,x) {
      pnorm(x,mean=theta[1],sd=theta[2])
    }
    v1 <- boot_gaussian(B=B,t0=t0,theta=theta_r,C=C_r,N=N1)
  } else if (m1=="logistic"){
    theta_r <- c(mod1$coefficients,mod1$scale)
    f1 <- function(theta,x) {
      dlogis(x,location=theta[1],scale=theta[2])
    }
    F1 <- function(theta,x) {
      plogis(x,location=theta[1],scale=theta[2])
    }
    v1 <- boot_logistic(B=B,t0=t0,theta=theta_r,C=C_r,N=N1)
  } else if (m1=="lognormal"){
    theta_r <- c(mod1$coefficients,mod1$scale)
    f1 <- function(theta,x) {
     dlnorm(x,meanlog=theta[1],sdlog=theta[2])
    }
    F1 <- function(theta,x) {
      plnorm(x,meanlog=theta[1],sdlog=theta[2])
    }
    v1 <- boot_lognormal(B=B,t0=t0,theta=theta_r,C=C_r,N=N1)
      } else if (m1=="loglogistic"){
        theta_r <- c(1/mod1$scale,exp(mod1$coefficients))
        f1 <- function(theta,x) {
          dllogis(x,shape=theta[1],scale=theta[2])
        }
        F1 <- function(theta,x) {
          pllogis(x,shape=theta[1],scale=theta[2])
        }
        v1 <- boot_loglogistic(B=B,t0=t0,theta=theta_r,C=C_r,N=N1)
        } else {return("m1: Invalid type of model")}


         if(m2=="weibull"){
          theta_t <- c(1/mod2$scale,exp(mod2$coefficients))
          f2 <- function(theta,x) {
            dweibull(x,shape=theta[1],scale=theta[2])
          }
          F2 <- function(theta,x) {
            pweibull(x,shape=theta[1],scale=theta[2])
          }
          v2 <- boot_weibull(B=B,t0=t0,theta=theta_t,C=C_t,N=N2)
        } else if (m2=="exponential"){
          theta_t <- 1/exp(mod2$coefficients)
          f2 <- function(theta,x) {
            dexp(x,rate=theta)
          }
          F2 <- function(theta,x) {
            pexp(x,rate=theta)
          }
          v2 <- boot_exponential(B=B,t0=t0,theta=theta_t,C=C_t,N=N2)
        } else if (m2=="gaussian"){
          theta_t <- c(mod2$coefficients,mod2$scale)
          f2 <- function(theta,x) {
            dnorm(x,mean=theta[1],sd=theta[2])
          }
          F2 <- function(theta,x) {
            pnorm(x,mean=theta[1],sd=theta[2])
          }
          v2 <- boot_gaussian(B=B,t0=t0,theta=theta_t,C=C_t,N=N2)
        } else if (m2=="logistic"){
          theta_t <- c(mod2$coefficients,mod2$scale)
          f2 <- function(theta,x) {
            dlogis(x,location=theta[1],scale=theta[2])
          }
          F2 <- function(theta,x) {
            plogis(x,location=theta[1],scale=theta[2])
          }
          v2 <- boot_logistic(B=B,t0=t0,theta=theta_t,C=C_t,N=N2)
        } else if (m2=="lognormal"){
          theta_t <- c(mod2$coefficients,mod2$scale)
          f2 <- function(theta,x) {
            dlnorm(x,meanlog=theta[1],sdlog=theta[2])
          }
          F2 <- function(theta,x) {
            plnorm(x,meanlog=theta[1],sdlog=theta[2])
          }
          v2 <- boot_lognormal(B=B,t0=t0,theta=theta_t,C=C_t,N=N2)
        } else if (m2=="loglogistic"){
          theta_t <- c(1/mod2$scale,exp(mod2$coefficients))
          f2 <- function(theta,x) {
            dllogis(x,shape=theta[1],scale=theta[2])
          }
          F2 <- function(theta,x) {
            pllogis(x,shape=theta[1],scale=theta[2])
          }
          v2 <- boot_loglogistic(B=B,t0=t0,theta=theta_t,C=C_t,N=N2)
          } else {return("m2: Invalid type of model")}

surv_r <- function(t){
            1-F1(theta_r,t)
          }
surv_t <- function(t){
            1-F2(theta_t,t)
          }

diff <- function(t){surv_r(t)-surv_t(t)}

t1 <- diff(t0)

up_df <- diff(t0)+qnorm(1-alpha)*sqrt(var(v1-v2,na.rm=TRUE))
lo_df <- diff(t0)-qnorm(1-alpha)*sqrt(var(v1-v2,na.rm=TRUE))

if(plot==TRUE){
  x <- NULL; rm(x);
  curve(1-F1(theta_r,x),xlim=c(0,max(data_r$time,data_t$time)),xlab="Time (days)",ylab="Survival",lwd=1.5)
  curve(1-F2(theta_t,x),xlim=c(0,max(data_r$time,data_t$time)),lwd=1.5,lty=2,add=TRUE)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topright", lty=1:2, lwd=1.5, legend=c("Reference treatment", "Test treatment"), cex=0.9)
}

return(list(diff=t1,lower.bound=lo_df,upper.bound=up_df,modell.1=mod1,modell.2=mod2))
}
