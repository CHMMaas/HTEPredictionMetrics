#' @title val_surv_mi
#' @description This function calculates intercept, slope, and C-index for risk predictions of multiple imputed data set(s).
#'
#' @importFrom survival coxph
#' @importFrom rms Predict
#' @importFrom Hmisc rcorr.cens
#'
#' @param est vector of estimates
#' @param se vector of standard errors
#'
#' @return The output of the val_surv_mi function is a "list" with the following components.
#'
#' est
#'
#' Combined estimate
#'
#'
#' se
#'
#' Combined standard error
#' @export
#'
#' @examples
#' library(HTEPredictionMetrics)
#' m <- 5 # number of imputations
#' est <- runif(m)
#' se <- runif(m)
#' Rubin.combine(est, se)
Rubin.combine<-function(est,se){
  m<-length(est)
  est.mi<-mean(est)
  var.w<-mean(se^2)
  var.b<-0
  if (m>1) var.b<-sum((est - est.mi)^2)/(m - 1)
  se.mi<-sqrt(var.w+(1+1/m)*var.b)
  return(list(est=est.mi,se=se.mi))
}