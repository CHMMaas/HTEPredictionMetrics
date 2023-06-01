#' @title E-for-benefit
#' @description This function calculates E-for-benefit statistics. Please note,
#' this function is only applicable for binary outcomes.
#'
#' @importFrom dplyr slice
#' @importFrom stats quantile
#' @importFrom stats loess
#' @importFrom stats predict
#'
#' @param Y a vector of binary outcomes; 1 if an (unfavourable) event; 0 if not
#' @param W a vector of treatment assignment; 1 for active treatment; 0 for control
#' @param X a matrix or data.frame of patient characteristics or individualized treatment effect predictions, categorical variables may be coded as.factor() to create dummy variables when matching, do not include Y or W in this matrix
#' @param p.0 a vector of outcome probabilities under control
#' @param p.1 a vector of outcome probabilities under active treatment
#' @param tau.hat a vector of individualized treatment effect predictions
#' @param matched.patients dataframe; optional if you want to provide your own dataframe (including matched.tau.hat, matched.tau.obs, and also include subclass if confidence interval needs to be computed) of matched patients, otherwise patients will be matched (default=NULL)
#' @param CI boolean; TRUE compute confidence interval; default=FALSE do not compute confidence interval (default=FALSE)
#' @param nr.bootstraps boolean; number of bootstraps to use for confidence interval computation (default=1)
#' @param message boolean; TRUE display computation time message; FALSE do not display message (default=TRUE)
#' @param measure measure option of matchit function from MatchIt package (default="nearest")
#' @param distance distance option of matchit function from MatchIt package (default="mahalanobis)
#' @param estimand default ATC meaning treated units are selected to be matched with control units
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
#' set.seed(1)
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
#'                         estimand=NULL)
#' EB.out
#'
#' # alternatively, use a dataframe of matched patients and calculate the C-for-Benefit
#' out.matched <- match.patients(Y=Y, W=W, X=X,
#'                               p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                               print=TRUE, measure="nearest",
#'                               distance="mahalanobis", estimand=NULL)
#' EB.out <- E.for.Benefit(matched.patients=out.matched$df.matched.pairs,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE)
#' EB.out
E.for.Benefit <- function(Y=NULL, W=NULL, X=NULL,
                          p.0=NULL, p.1=NULL, tau.hat=NULL,
                          matched.patients=NULL,
                          CI=FALSE, nr.bootstraps=50, message=TRUE,
                          measure="nearest", distance="mahalanobis",
                          estimand=NULL, ...){
  # check user input
  stopifnot("nr.bootstraps must be a scalar" = length(nr.bootstraps)==1)
  stopifnot("nr.bootstraps must be numeric" = is.numeric(nr.bootstraps))
  stopifnot("CI must be a boolean (TRUE or FALSE)" = isTRUE(CI)|isFALSE(CI))
  stopifnot("message must be a boolean (TRUE or FALSE)" = isTRUE(message)|isFALSE(message))

  if (is.null(matched.patients)){
    # check user input
    stopifnot("Y must be numeric" = is.numeric(Y))
    stopifnot("W must be numeric" = is.numeric(W))
    stopifnot("p.0 must be numeric" = is.numeric(p.0))
    stopifnot("p.1 must be numeric" = is.numeric(p.1))
    stopifnot("tau.hat must be numeric" = is.numeric(tau.hat))

    stopifnot("W must be a vector" = is.vector(W))
    stopifnot("W must only consists of zeros and ones" = !sum(sort(unique(W))-c(0, 1)))
    stopifnot("p.0 must be a vector" = is.vector(p.0))
    stopifnot("p.1 must be a vector" = is.vector(p.1))
    stopifnot("tau.hat must be a vector" = is.vector(tau.hat))

    stopifnot("Y and W must be the same length" = length(Y)==length(W))
    stopifnot("Y and X must have the same number of observations" = length(Y)==nrow(X))
    stopifnot("Y and p.0 must be the same length" = length(Y)==length(p.0))
    stopifnot("Y and p.1 must be the same length" = length(Y)==length(p.1))
    stopifnot("Y and tau.hat must be the same length" = length(Y)==length(tau.hat))

    stopifnot("W and X must have the same number of observations" = length(W)==nrow(X))
    stopifnot("W and p.0 must be the same length" = length(W)==length(p.0))
    stopifnot("W and p.1 must be the same length" = length(W)==length(p.1))
    stopifnot("W and tau.hat must be the same length" = length(W)==length(tau.hat))

    stopifnot("X and p.0 must have the same number of observations" = length(p.0)==nrow(X))
    stopifnot("X and p.1 must have the same number of observations" = length(p.1)==nrow(X))
    stopifnot("X and tau.hat must have the same number of observations" = length(tau.hat)==nrow(X))

    stopifnot("p.0 and p.1 must be the same length" = length(p.1)==length(p.0))
    stopifnot("p.0 and tau.hat must be the same length" = length(p.1)==length(p.0))

    stopifnot("p.1 and tau.hat must be the same length" = length(p.1)==length(p.0))

    # compute matched pairs
    matched.patients <- match.patients(Y=Y, W=W, X=X,
                                       p.0=p.0, p.1=p.1, tau.hat=tau.hat,
                                       measure=measure, distance=distance,
                                       estimand=estimand, ...)$df.matched.pairs
  }
  else{
    # use the dataframe provided by the user
    stopifnot("matched.patients must be a dataframe" = is.data.frame(matched.patients))
  }

  # perform smoothing on matched patient pairs
  loess.calibrate <- stats::loess(matched.tau.obs ~ matched.tau.hat,
                                  data=matched.patients)
  matched.patients$tau.smoothed <- predict(loess.calibrate)

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
    # bootstrap matched patient pairs
    for (B in 1:nr.bootstraps){
      # obtain all subclass IDs
      subclass.IDs <- sort(unique(matched.patients$subclass))
      # sample randomly from subclass IDs
      sample.subclass <- sample(subclass.IDs, length(subclass.IDs), replace=TRUE)
      # data frame of duplicated patients
      matched.patients.dup <- dplyr::slice(matched.patients, sample.subclass)

      # perform smoothing on matched patient pairs
      loess.calibrate.B <- stats::loess(matched.tau.obs ~ matched.tau.hat,
                                        data=matched.patients.dup)
      tau.smoothed.B <- stats::predict(loess.calibrate, newdata=matched.patients.dup)

      # calculate calibration metrics
      Eavg.for.CI <- c(Eavg.for.CI, mean(abs(tau.smoothed.B - matched.patients.dup$matched.tau.hat)))
      E50.for.CI <- c(E50.for.CI, stats::median(abs(tau.smoothed.B - matched.patients.dup$matched.tau.hat)))
      E90.for.CI <- c(E90.for.CI, as.numeric(stats::quantile(abs(tau.smoothed.B - matched.patients.dup$matched.tau.hat), probs=0.9)))
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
