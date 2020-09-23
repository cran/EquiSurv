################################################################################
#' @title Non-inferiority and equivalence test for the difference of two parametric survival curves
#'
#' @description Function for fitting and testing two parametric survival curves \eqn{S_1}, \eqn{S_2} at \eqn{t_0} concerning the
#'hypotheses of non-inferiority \deqn{H_0:S_1(t_0)-S_2(t_0)\geq \epsilon\ vs.\ H_1: S_1(t_0)-S_2(t_0)< \epsilon}
#'or equivalence \deqn{H_0:|S_1(t_0)-S_2(t_0)|\geq \epsilon\ vs.\ H_1: |S_1(t_0)-S_2(t_0)|< \epsilon.}
#'\eqn{m_1} and \eqn{m_2} are parametric survival models following a Weibull, exponential, gaussian, logistic, log-normal or log-logistic distribution.
#'The test procedure is based on confidence intervals obtained via bootstrap.
#'For the generation of the bootstrap data exponentially distributed random censoring is assumed and the rates estimated from the datasets.
#'See Moellenhoff and Tresch <arXiv:2009.06699> for details.
#'
#' @name test_diff
#' @export
#' @import survival eha
#' @importFrom stats dexp dlnorm dlogis dnorm dweibull pexp plnorm plogis pnorm pweibull qnorm rexp rlnorm rlogis rnorm rweibull var
#' @param epsilon non-inferiority/equivalence margin
#' @param alpha significance level
#' @param t0 time point of interest
#' @param type type of the test. "ni" for non-inferiority, "eq" for equivalence test
#' @param m1,m2 type of parametric model. Possible model types are "weibull", "exponential", "gaussian", "logistic", "lognormal" and "loglogistic"
#' @param B number of bootstrap repetitions. The default is B=1000
#' @param data_r,data_t datasets containing time and status for each individual (have to be referenced as this)
#' @param plot if TRUE, a plot of the two survival curves will be given
#' @return A list containing the difference \eqn{S_1(t_0)-S_2(t_0)}, the lower and upper (1-\eqn{\alpha})-confidence bounds, the summary of the two model fits, the chosen margin and significance level and the test decision. Further a plot of the curves is given.
#' @examples
#' data(veteran)
#' veteran_r <- veteran[veteran$trt==1,]
#' veteran_t <- veteran[veteran$trt==2,]
#' alpha<-0.05
#' t0<-80
#' epsilon<-0.15
#' test_diff(epsilon=epsilon,alpha=alpha,t0=t0,type="eq",m1="weibull",m2="weibull",
#' data_r=veteran_r,data_t=veteran_t)
#' @references K.Moellenhoff and A.Tresch: Survival analysis under non-proportional hazards: investigating non-inferiority or equivalence in time-to-event data <arXiv:2009.06699>
################################################################################

test_diff <- function(epsilon,alpha,t0,type,m1,m2,B=1000,plot=TRUE,data_r,data_t){

  conf_int <- confint_diff(alpha=alpha,t0=t0,m1=m1,m2=m2,B,plot,data_r=data_r,data_t=data_t)

  if (type=="ni") {
    if(conf_int$upper.bound <= epsilon){decision="reject H0 (non-inferiority)"}else{decision="don't reject H0"}
  }
  else if (type=="eq") {
    if(-epsilon <= conf_int$lower.bound & conf_int$upper.bound <= epsilon){decision="reject H0 (equivalent)"}else{decision="don't reject H0 (not equivalent)"}
  } else {return("Type is not valid, choose ni or eq")}

  return(list(conf_int=conf_int,margin=epsilon,level=alpha,decision=decision))
}
