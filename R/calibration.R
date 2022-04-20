#' @title E-for-benefit
#' @description This function calculates E-for-benefit statistics.
#'
#' @importFrom dplyr slice
#' @importFrom stats quantile
#' @importFrom stats loess
#' @importFrom stats predict
#'
#' @param Y a vector of outcomes
#' @param W a vector of treatment assignment; 1 for active treatment; 0 for control
#' @param X a matrix of patient characteristics or individualized treatment effect predictions
#' @param p.0 a vector of outcome probabilities under control
#' @param p.1 a vector of outcome probabilities under active treatment
#' @param tau.hat a vector of individualized treatment effect predictions
#' @param CI boolean; TRUE compute confidence interval; default=FALSE do not compute confidence interval (default=FALSE)
#' @param nr.bootstraps boolean; number of bootstraps to use for confidence interval computation (default=1)
#' @param message boolean; TRUE display computation time message; FALSE do not display message (default=TRUE)
#' @param matched.patients dataframe; optional if you want to provide your own dataframe of matched patients, otherwise patients will be matched (default=NULL)
#' @param measure measure option of matchit function from MatchIt package (default="nearest")
#' @param distance distance option of matchit function from MatchIt package (default="mahalanobis)
#' @param estimand default ATC meaning treated units are selected to be matched with control units
#' @param replace boolean; TRUE if matching with replacement, FALSE if matching without replacement
#' @param ... additional arguments for matchit function from MatchIt package
#'
#' @return The output of the E.for.benefit function is a "list" with the following components.
#'
#' matched.patients
#'
#' a dataframe containing the matched patients.
#'
#'
#' Eavg.for.benefit
#'
#' the resulting Eavg-for-Benefit value.
#'
#'
#' Eavg.lower.CI
#'
#' the lower bound of the confidence interval of the Eavg-for-Benefit (if CI = TRUE).
#'
#'
#' Eavg.upper.CI
#'
#' the upper bound of the confidence interval of the Eavg-for-Benefit (if CI = TRUE).
#'
#'
#' E50.for.benefit
#'
#' the resulting E50-for-Benefit value.
#'
#'
#' E50.lower.CI
#'
#' the lower bound of the confidence interval of the E50-for-Benefit (if CI = TRUE).
#'
#'
#' E50.upper.CI
#'
#' the upper bound of the confidence interval of the E50-for-Benefit (if CI = TRUE).
#'
#'
#' E90.for.benefit
#'
#' the resulting E90-for-Benefit value.
#'
#'
#' E90.lower.CI
#'
#' the lower bound of the confidence interval of the E90-for-Benefit (if CI = TRUE).
#'
#'
#' E90.upper.CI
#'
#' the upper bound of the confidence interval of the E90-for-Benefit (if CI = TRUE).
#' @export
#'
#' @examples
#' library(HTEPredictionMetrics)
#' n <- 100
#' Y <- sample(0:1, n, replace=TRUE)
#' W <- sample(0:1, n, replace=TRUE)
#' X <- matrix(rnorm(n), n, 4)
#' p.0 <- runif(n)
#' p.1 <- runif(n)
#' tau.hat <- runif(n)
#' EB.out <- E.for.Benefit(Y=Y, W=W, X=X, p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE,
#'                         matched.patients=NULL,
#'                         measure="nearest", distance="mahalanobis",
#'                         estimand=NULL, replace=FALSE)
#' EB.out
E.for.Benefit <- function(Y, W, X,
                          p.0, p.1, tau.hat,
                          CI=FALSE, nr.bootstraps=50, message=TRUE,
                          matched.patients=NULL,
                          measure="nearest", distance="mahalanobis",
                          estimand=NULL, replace=FALSE, ...){
  # ensure correct data types
  stopifnot("W must be a vector" = is.vector(W))
  stopifnot("X must be a vector or matrix" = is.matrix(X) | is.vector(X))
  stopifnot("p.0 must be a vector" = is.vector(p.0))
  stopifnot("p.1 must be a vector" = is.vector(p.1))
  stopifnot("tau.hat must be a vector" = is.vector(tau.hat))
  stopifnot("nr.bootstraps must be a scalar" = length(nr.bootstraps)==1)

  stopifnot("Y must be numeric" = is.numeric(Y))
  stopifnot("W must be numeric" = is.numeric(W))
  stopifnot("X must be numeric" = is.numeric(X))
  stopifnot("p.0 must be numeric" = is.numeric(p.0))
  stopifnot("p.1 must be numeric" = is.numeric(p.1))
  stopifnot("tau.hat must be numeric" = is.numeric(tau.hat))
  stopifnot("nr.bootstraps must be numeric" = is.numeric(nr.bootstraps))

  stopifnot("W must only consists of zeros and ones" = !sum(sort(unique(W))-c(0, 1)))
  stopifnot("CI must be a boolean (TRUE or FALSE)" = isTRUE(CI)|isFALSE(CI))
  stopifnot("message must be a boolean (TRUE or FALSE)" = isTRUE(message)|isFALSE(message))

  if (is.null(matched.patients)){
    # compute matched pairs
    matched.patients <- match.patients(Y, W, X,
                                       p.0, p.1, tau.hat,
                                       measure="nearest", distance="mahalanobis",
                                       estimand=NULL, replace=FALSE, ...)$matched.patients
  }
  else{
    # use the dataframe provided by the user
    stopifnot("matched.patients must be a dataframe" = is.data.frame(matched.patients))
  }

  # perform smoothing on matched patient pairs
  loess.calibrate <- stats::loess(matched.tau.obs ~ matched.tau.hat,
                                  data=matched.patients)

  # calculate calibration metrics
  Eavg.for.benefit <- mean(abs(matched.patients$tau.smoothed - matched.patients$matched.tau.hat))
  E50.for.benefit <- stats::median(abs(matched.patients$tau.smoothed - matched.patients$matched.tau.hat))
  E90.for.benefit <- as.numeric(stats::quantile(abs(matched.patients$tau.smoothed - matched.patients$matched.tau.hat), probs=0.9))

  if (CI){
    if (message){
      cat('Calculating confidence interval... Taking too long? Lower the number of bootstraps. \n')
    }
    Eavg.for.CI <- c()
    E50.for.CI <- c()
    E90.for.CI <- c()
    for (B in 1:nr.bootstraps){
      # bootstrap matched patient pairs
      subclass.IDs <- unique(matched.patients$subclass)
      sample.subclass <- sample(subclass.IDs, length(subclass.IDs), replace=TRUE)
      dup.subclass.IDs <- c()
      for (i in sample.subclass){
        dup.subclass.IDs <- c(dup.subclass.IDs, matched.patients[matched.patients$subclass==i, 'match.id'])
      }
      duplicated.matched.patients <- dplyr::slice(matched.patients, dup.subclass.IDs)

      # perform smoothing on matched patient pairs
      loess.calibrate.B <- stats::loess(matched.tau.obs ~ matched.tau.hat,
                                        data=duplicated.matched.patients)
      tau.smoothed.B <- stats::predict(loess.calibrate, newdata=duplicated.matched.patients)

      # calculate calibration metrics
      Eavg.for.CI <- c(Eavg.for.CI, mean(abs(tau.smoothed.B - duplicated.matched.patients$matched.tau.hat)))
      E50.for.CI <- c(E50.for.CI, stats::median(abs(tau.smoothed.B - duplicated.matched.patients$matched.tau.hat)))
      E90.for.CI <- c(E90.for.CI, as.numeric(stats::quantile(abs(tau.smoothed.B - duplicated.matched.patients$matched.tau.hat), probs=0.9)))
    }
    Eavg.lower.CI <- as.numeric(stats::quantile(Eavg.for.CI, 0.025))
    Eavg.upper.CI <- as.numeric(stats::quantile(Eavg.for.CI, 0.975))
    E50.lower.CI <- as.numeric(stats::quantile(E50.for.CI, 0.025))
    E50.upper.CI <- as.numeric(stats::quantile(E50.for.CI, 0.975))
    E90.lower.CI <- as.numeric(stats::quantile(E90.for.CI, 0.025))
    E90.upper.CI <- as.numeric(stats::quantile(E90.for.CI, 0.975))
  }
  else{
    Eavg.lower.CI <- NA
    Eavg.upper.CI <- NA
    E50.lower.CI <- NA
    E50.upper.CI <- NA
    E90.lower.CI <- NA
    E90.upper.CI <- NA
  }

  return(list(matched.patients=matched.patients,
              Eavg.for.benefit=Eavg.for.benefit, Eavg.lower.CI=Eavg.lower.CI, Eavg.upper.CI=Eavg.upper.CI,
              E50.for.benefit=E50.for.benefit, E50.lower.CI=E50.lower.CI, E50.upper.CI=E50.upper.CI,
              E90.for.benefit=E90.for.benefit, E90.lower.CI=E90.lower.CI, E90.upper.CI=E90.upper.CI))
}
