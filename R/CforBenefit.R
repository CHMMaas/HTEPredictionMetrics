#' @title C-for-benefit
#' @description This function calculates the C-for-benefit, as proposed by
#' D. van Klaveren et al. (2018), which measures the discriminative ability of
#' models predicting individualized treatment effect. The C-for-benefit
#' corresponds to the probability that from two randomly chosen matched patient
#' pairs with unequal observed treatment effect, the pair with greater observed
#' treatment effect also has a higher predicted treatment effect. Please note,
#' this function is only applicable for binary outcomes.
#'
#' @importFrom Hmisc rcorr.cens
#' @importFrom dplyr slice
#' @importFrom stats quantile
#'
#' @param Y a vector of binary outcomes; 1 if an (unfavourable) event; 0 if not
#' @param W a vector of treatment assignment; 1 for active treatment; 0 for control (or alternative treatment), note: make sure that W consists only of values of 0 and 1
#' @param X a matrix of patient characteristics or individualized treatment effect predictions, do not include Y or W in this matrix
#' @param p.0 a vector of outcome probabilities under control (or alternative treatment)
#' @param p.1 a vector of outcome probabilities under active treatment
#' @param tau.hat a vector of individualized treatment effect predictions
#' @param matched.patients dataframe; optional if you want to provide your own dataframe (including matched.tau.hat, matched.tau.obs, and also include subclass if confidence interval needs to be computed) of matched patients, otherwise patients will be matched (default=NULL)
#' @param CI boolean; TRUE compute confidence interval; default=FALSE do not compute confidence interval (default=FALSE)
#' @param nr.bootstraps boolean; number of bootstraps to use for confidence interval computation (default=1)
#' @param message boolean; TRUE display computation time message; FALSE do not display message (default=TRUE)
#' @param measure measure option of matchit function from MatchIt package (default="nearest")
#' @param distance distance option of matchit function from MatchIt package (default="mahalanobis)
#' @param estimand default ATC meaning treated units are selected to be matched with control units
#' @param replace boolean; TRUE if matching with replacement, FALSE if matching without replacement
#' @param ... additional arguments for matchit function from MatchIt package
#'
#' @return The output of the C.for.Benefit function is a "list" with the following components.
#'
#' matched.patients
#'
#' a dataframe containing the matched patients.
#'
#'
#' C.for.benefit
#'
#' the resulting C-for-benefit value.
#'
#'
#' lower.CI
#'
#' the lower bound of the confidence interval (if CI = TRUE).
#'
#'
#'  upper.CI
#'
#' the upper bound of the confidence interval (if CI = TRUE).
#' @export
#'
#' @examples
#' library(HTEPredictionMetrics)
#' set.seed(1)
#' n <- 100
#' Y <- sample(0:1, n, replace=TRUE)
#' W <- sample(0:1, n, replace=TRUE)
#' X <- matrix(rnorm(n), n, 3)
#' p.0 <- runif(n)
#' p.1 <- runif(n)
#' tau.hat <- runif(n)
#' CB.out <- C.for.Benefit(Y=Y, W=W, X=X, p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE,
#'                         matched.patients=NULL,
#'                         measure="nearest", distance="mahalanobis",
#'                         estimand=NULL, replace=FALSE)
#' CB.out
#'
#' # alternatively, use a dataframe of matched patients and calculate the C-for-Benefit
#' out.matched <- match.patients(Y=Y, W=W, X=X,
#'                               p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                               CI=FALSE, message=TRUE,
#'                               measure="nearest", distance="mahalanobis",
#'                               estimand=NULL, replace=FALSE)
#' CB.out <- C.for.Benefit(matched.patients=out.matched$df.matched.pairs,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE, replace=FALSE)
#' CB.out
C.for.Benefit <- function(Y=NULL, W=NULL, X=NULL,
                          p.0=NULL, p.1=NULL, tau.hat=NULL,
                          matched.patients=NULL,
                          CI=FALSE, nr.bootstraps=50, message=TRUE,
                          measure="nearest", distance="mahalanobis",
                          estimand=NULL, replace=FALSE, ...){
  # check user input
  stopifnot("nr.bootstraps must be a scalar" = length(nr.bootstraps)==1)
  stopifnot("nr.bootstraps must be numeric" = is.numeric(nr.bootstraps))
  stopifnot("CI must be a boolean (TRUE or FALSE)" = isTRUE(CI)|isFALSE(CI))
  stopifnot("message must be a boolean (TRUE or FALSE)" = isTRUE(message)|isFALSE(message))

  if (is.null(matched.patients)){
    # check user input
    stopifnot("Y must be numeric" = is.numeric(Y))
    stopifnot("W must be numeric" = is.numeric(W))
    stopifnot("X must be numeric" = is.numeric(as.matrix(X)))
    stopifnot("p.0 must be numeric" = is.numeric(p.0))
    stopifnot("p.1 must be numeric" = is.numeric(p.1))
    stopifnot("tau.hat must be numeric" = is.numeric(tau.hat))

    stopifnot("W must be a vector" = is.vector(W))
    stopifnot("W must only consists of zeros and ones" = !sum(sort(unique(W))-c(0, 1)))
    stopifnot("X must be a vector or matrix, use as.matrix(X) instead" = is.vector(X) | is.matrix(X))
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
                                       estimand=estimand, replace=replace, ...)$df.matched.pairs
  }
  else{
    # use the dataframe provided by the user
    stopifnot("matched.patients must be a dataframe of matched pairs" = is.data.frame(matched.patients))
  }

  # calculate C-for-benefit
  cindex <- Hmisc::rcorr.cens(matched.patients$matched.tau.hat, matched.patients$matched.tau.obs)
  C.for.benefit <- cindex["C Index"][[1]]
  if (CI){
    if (message){
      cat('Calculating confidence interval... Taking too long? Lower the number of bootstraps. \n')
    }
    CB.for.CI <- c()
    # bootstrap matched patient pairs
    for (B in 1:nr.bootstraps){
      # obtain all subclass IDs
      subclass.IDs <- sort(unique(matched.patients$subclass))
      # sample randomly from subclass IDs
      sample.subclass <- sample(subclass.IDs, length(subclass.IDs), replace=TRUE)
      # data frame of duplicated patients
      matched.patients.dup <- dplyr::slice(matched.patients, sample.subclass)

      # calculate C-for-benefit for duplicated matched pairs
      duplicated.cindex <- Hmisc::rcorr.cens(matched.patients.dup$matched.tau.hat, matched.patients.dup$matched.tau.obs)
      CB.for.CI <- c(CB.for.CI, duplicated.cindex["C Index"][[1]])
    }
    lower.CI <- as.numeric(stats::quantile(CB.for.CI, 0.025))
    upper.CI <- as.numeric(stats::quantile(CB.for.CI, 0.975))
  }
  else{
    lower.CI <- NA
    upper.CI <- NA
  }

  return(list(df.matched.patients=matched.patients, C.for.benefit=C.for.benefit, lower.CI=lower.CI, upper.CI=upper.CI))
}
