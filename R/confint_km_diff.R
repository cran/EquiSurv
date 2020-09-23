################################################################################
#' @title Lower and upper confidence bounds for the difference of two Kaplan-Meier curves
#'
#' @description Function fitting Kaplan-Meier curves \eqn{S_1}, \eqn{S_2} to two groups and yielding
#' lower and upper (1-\eqn{\alpha})-confidence bounds for the difference \eqn{S_1-S_2} of these
#' two curves at a specific time point by using Greenwood's formula.
#'
#' @name confint_km_diff
#' @export
#' @import survival
#' @param alpha confidence level
#' @param t0 time point of interest
#' @param data_r,data_t datasets containing time and status for each individual
#' @param plot if TRUE, a plot of the two Kaplan Meier curves will be given
#' @return A list containing the difference \eqn{S_1(t_0)-S_2(t_0)} and the lower and upper (1-\eqn{\alpha})-confidence bounds. Further a plot of the curves is given.
#' @examples
#' data(veteran)
#' veteran_r <- veteran[veteran$trt==1,]
#' veteran_t <- veteran[veteran$trt==2,]
#' alpha<-0.05
#' t0<-80
#' confint_km_diff(alpha=alpha,t0=t0,data_r=veteran_r,data_t=veteran_t)
################################################################################

confint_km_diff <- function(alpha,t0,data_r,data_t,plot=TRUE){

    km_fit1 <- survfit(Surv(time,status)~ 1, type="kaplan-meier", data = data_r)
    km_fit2 <- survfit(Surv(time,status)~ 1, type="kaplan-meier", data = data_t)

      var1 <- summary(km_fit1, times = t0)[7]$std.err^2
      var2 <- summary(km_fit2, times = t0)[7]$std.err^2

      t1 <- summary(km_fit1, times = t0)[6]$surv-summary(km_fit2, times = t0)[6]$surv
      conf_int_l <- t1 - qnorm(1-alpha)*sqrt(var1+var2)
      conf_int_u <- t1 + qnorm(1-alpha)*sqrt(var1+var2)

      if(plot==TRUE){
        plot(km_fit1, xlab="Time (days)", ylab="Survival", xlim=c(0,1000),conf.int=FALSE)
        par(new=TRUE)
        plot(km_fit2, axes=FALSE,lty=2,conf.int=FALSE)
        grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
        legend(x="topright", lty=1:2, lwd=1.5, legend=c("Reference treatment", "Test treatment"), cex=0.9)
      }
      return(list(diff=t1,lower.bound=conf_int_l,upper.bound=conf_int_u))
}
