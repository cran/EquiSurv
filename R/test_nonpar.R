################################################################################
#' @title Non-inferiority and equivalence test for the difference of two Kaplan-Meier curves
#'
#' @description Function for fitting and testing two Kaplan Meier curves \eqn{S_1}, \eqn{S_2} at \eqn{t_0} concerning the
#'hypotheses of non-inferiority \deqn{H_0:S_1(t_0)-S_2(t_0)\geq \epsilon\ vs.\ H_1: S_1(t_0)-S_2(t_0)< \epsilon}
#'or equivalence \deqn{H_0:|S_1(t_0)-S_2(t_0)|\geq \epsilon\ vs.\ H_1: |S_1(t_0)-S_2(t_0)|< \epsilon.}
#'
#' @name test_nonpar
#' @export
#' @import survival
#' @param epsilon non-inferiority/equivalence margin
#' @param alpha significance level
#' @param t0 time point of interest
#' @param type type of the test. "ni" for non-inferiority, "eq" for equivalence test
#' @param data_r,data_t datasets containing time and status for each individual
#' @param plot if TRUE, a plot of the two Kaplan Meier curves will be given
#' @return A list containing the difference \eqn{S_1(t_0)-S_2(t_0)}, the lower and upper (1-\eqn{\alpha})-confidence bounds, the chosen margin and significance level and the test decision. Further a plot of the curves is given.
#' @examples
#' data(veteran)
#' veteran_r <- veteran[veteran$trt==1,]
#' veteran_t <- veteran[veteran$trt==2,]
#' alpha<-0.05
#' t0<-80
#' epsilon<-0.15
#' test_nonpar(epsilon=epsilon,alpha=alpha,t0=t0,type="eq",data_r=veteran_r,data_t=veteran_t)
################################################################################

test_nonpar <- function(epsilon,alpha,t0,type,data_r,data_t,plot=TRUE){

  conf_int <- confint_km_diff(alpha=alpha,t0=t0,data_r=data_r,data_t=data_t)

  if (type=="ni") {
    if(conf_int$upper.bound <= epsilon){decision="reject H0 (non-inferiority)"}else{decision="don't reject H0"}
  }
  else if (type=="eq") {
    if(-epsilon <= conf_int$lower.bound & conf_int$upper.bound <= epsilon){decision="reject H0 (equivalent)"}else{decision="don't reject H0 (not equivalent)"}
  } else {return("Type is not valid, choose ni or eq")}

  return(list(conf_int=conf_int,margin=epsilon,level=alpha,decision=decision))
}
