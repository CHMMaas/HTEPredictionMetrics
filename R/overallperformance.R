#' @title Overall performance metrics for benefit
#' @description This function calculates logistic-loss-for-benefit and Brier-for-benefit.
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
#' @param matched.patients dataframe; optional if you want to provide your own dataframe (including p.0, p.1, matched.tau.hat, matched.tau.obs, and also include subclass if confidence interval needs to be computed) of matched patients, otherwise patients will be matched (default=NULL)
#' @param CI boolean; TRUE compute confidence interval; default=FALSE do not compute confidence interval (default=FALSE)
#' @param nr.bootstraps boolean; number of bootstraps to use for confidence interval computation (default=1)
#' @param message boolean; TRUE display computation time message; FALSE do not display message (default=TRUE)
#' @param measure measure option of matchit function from MatchIt package (default="nearest")
#' @param distance distance option of matchit function from MatchIt package (default="mahalanobis)
#' @param estimand default ATC meaning treated units are selected to be matched with control units
#' @param replace boolean; TRUE if matching with replacement, FALSE if matching without replacement
#' @param ... additional arguments for matchit function from MatchIt package
#'
#' @return The output of the OP.for.Benefit function is a "list" with the following components.
#'
#' matched.patients
#'
#' a dataframe containing the matched patients.
#'
#'
#' Log.Loss.for.Benefit
#'
#' the resulting logistic-loss-for-Benefit value.
#'
#'
#' Log.Loss.lower.CI
#'
#' the lower bound of the confidence interval of the logistic-loss-for-Benefit (if CI = TRUE).
#'
#'
#' Log.Loss.upper.CI
#'
#' the upper bound of the confidence interval of the logistic-loss-for-Benefit (if CI = TRUE).
#'
#'
#' Brier.for.Benefit
#'
#' the resulting Brier-for-Benefit value.
#'
#'
#' Brier.lower.CI
#'
#' the lower bound of the confidence interval of the Brier-for-Benefit (if CI = TRUE).
#'
#'
#' Brier.upper.CI
#'
#' the upper bound of the confidence interval of the Brier-for-Benefit (if CI = TRUE).
#'
#' @export
#'
#' @examples
#' library(HTEPredictionMetrics)
#' n <- 100
#' Y <- sample(0:1, n, replace=TRUE)
#' W <- sample(0:1, n, replace=TRUE)
#' X <- matrix(rnorm(n), n, 3)
#' p.0 <- runif(n)
#' p.1 <- runif(n)
#' tau.hat <- runif(n)
#' OP.out <- OP.for.Benefit(Y=Y, W=W, X=X, p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE,
#'                         matched.patients=NULL,
#'                         measure="nearest", distance="mahalanobis",
#'                         estimand=NULL, replace=FALSE)
#' OP.out
#'
#' # alternatively, use a dataframe of matched patients and calculate the overall performance metrics
#' out.matched <- match.patients(Y=Y, W=W, X=X,
#'                               p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                               CI=FALSE, nr.bootstraps=50, message=TRUE,
#'                               measure="nearest", distance="mahalanobis",
#'                               estimand=NULL, replace=FALSE)
#' OP.out <- OP.for.Benefit(matched.patients=out.matched$df.matched.patients,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE, replace=FALSE)
#' OP.out
OP.for.Benefit <- function(Y=NULL, W=NULL, X=NULL,
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
    stopifnot("W must be a vector" = is.vector(W))
    stopifnot("X must be a vector or matrix" = is.matrix(X) | is.vector(X))
    stopifnot("p.0 must be a vector" = is.vector(p.0))
    stopifnot("p.1 must be a vector" = is.vector(p.1))
    stopifnot("tau.hat must be a vector" = is.vector(tau.hat))

    stopifnot("Y must be numeric" = is.numeric(Y))
    stopifnot("W must be numeric" = is.numeric(W))
    stopifnot("X must be numeric" = is.numeric(X))
    stopifnot("p.0 must be numeric" = is.numeric(p.0))
    stopifnot("p.1 must be numeric" = is.numeric(p.1))
    stopifnot("tau.hat must be numeric" = is.numeric(tau.hat))

    # compute matched pairs
    matched.patients <- match.patients(Y=Y, W=W, X=X,
                                        p.0=p.0, p.1=p.1, tau.hat=tau.hat,
                                        measure=measure, distance=distance,
                                        estimand=estimand, replace=replace, ...)$df.matched.patients
  }
  else{
    # use the dataframe provided by the user
    stopifnot("matched.patients must be a dataframe" = is.data.frame(matched.patients))
  }

  # set up quantities
  t.1 <- (1-matched.patients$p.1)*matched.patients$p.0
  t.0 <- (1-matched.patients$p.1)*(1-matched.patients$p.0) + matched.patients$p.1*matched.patients$p.0
  t.min1 <- matched.patients$p.1*(1-matched.patients$p.0)

  I.1 <- matched.patients$matched.tau.obs==1
  I.0 <- matched.patients$matched.tau.obs==0
  I.min1 <- matched.patients$matched.tau.obs==-1

  # Brier score for benefit
  Brier.for.Benefit <- (sum((t.1-I.1)^2)
                        +sum((t.0-I.0)^2)
                        +sum((t.min1-I.min1)^2))/(2*length(matched.patients$matched.tau.obs))

  # Logistic loss for benefit
  omit <- which(t.1<0 | t.0<0 | t.min1<0)
  if (length(omit)>0 & message){
    cat('nr. omitted observations for log loss:', length(omit), '\n')
    Log.Loss.for.Benefit <- -(sum(I.1[-omit]*log(t.1[-omit]))+sum(I.0[-omit]*log(t.0[-omit]))+sum(I.min1[-omit]*log(t.min1[-omit])))/length(matched.patients$matched.tau.obs)
  }
  else{
    Log.Loss.for.Benefit <- -(sum(I.1*log(t.1))+sum(I.0*log(t.0))+sum(I.min1*log(t.min1)))/length(matched.patients$matched.tau.obs)
  }

  if (CI){
    if (message){
      cat('Calculating confidence interval... Taking too long? Lower the number of bootstraps. \n')
    }
    Log.Loss.for.CI <- c()
    Brier.for.CI <- c()
    for (B in 1:nr.bootstraps){
      # bootstrap matched patient pairs
      subclass.IDs <- unique(matched.patients$subclass)
      sample.subclass <- sample(subclass.IDs, length(subclass.IDs), replace=TRUE)
      dup.subclass.IDs <- c()
      for (i in sample.subclass){
        dup.subclass.IDs <- c(dup.subclass.IDs, matched.patients[matched.patients$subclass==i, 'match.id'])
      }
      duplicated.matched.patients <- dplyr::slice(matched.patients, dup.subclass.IDs)

      # set up quantities
      t.1.B <- (1-duplicated.matched.patients$p.1)*duplicated.matched.patients$p.0
      t.0.B <- (1-duplicated.matched.patients$p.1)*(1-duplicated.matched.patients$p.0) + duplicated.matched.patients$p.1*duplicated.matched.patients$p.0
      t.min1.B <- duplicated.matched.patients$p.1*(1-duplicated.matched.patients$p.0)

      I.1.B <- duplicated.matched.patients$matched.tau.obs==1
      I.0.B <- duplicated.matched.patients$matched.tau.obs==0
      I.min1.B <- duplicated.matched.patients$matched.tau.obs==-1

      # Brier score for benefit
      Brier.for.Benefit.B <- (sum((t.1.B-I.1.B)^2)
                            +sum((t.0.B-I.0.B)^2)
                            +sum((t.min1.B-I.min1.B)^2))/(2*length(duplicated.matched.patients$matched.tau.obs))

      # Logistic loss for benefit
      omit.B <- which(t.1.B<0 | t.0.B<0 | t.min1.B<0)
      if (length(omit.B)>0 & message){
        cat('nr. omitted observations for log loss:', length(omit.B), '\n')
        Log.Loss.for.Benefit.B <- -(sum(I.1.B[-omit.B]*log(t.1.B[-omit.B]))
                                    +sum(I.0.B[-omit.B]*log(t.0.B[-omit.B]))
                                    +sum(I.min1.B[-omit.B]*log(t.min1[-omit.B])))/length(duplicated.matched.patients$matched.tau.obs)
      }
      else{
        Log.Loss.for.Benefit.B <- -(sum(I.1.B*log(t.1.B))
                                  +sum(I.0.B*log(t.0.B))
                                  +sum(I.min1.B*log(t.min1.B)))/length(duplicated.matched.patients$matched.tau.obs)
      }

      # calculate calibration metrics
      Log.Loss.for.CI <- c(Log.Loss.for.CI, Log.Loss.for.Benefit.B)
      Brier.for.CI <- c(Brier.for.CI, Brier.for.Benefit.B)
    }
    Log.Loss.lower.CI <- as.numeric(stats::quantile(Log.Loss.for.CI, 0.025))
    Log.Loss.upper.CI <- as.numeric(stats::quantile(Log.Loss.for.CI, 0.975))
    Brier.lower.CI <- as.numeric(stats::quantile(Brier.for.CI, 0.025))
    Brier.upper.CI <- as.numeric(stats::quantile(Brier.for.CI, 0.975))
  }
  else{
    Log.Loss.lower.CI <- NA
    Log.Loss.upper.CI <- NA
    Brier.lower.CI <- NA
    Brier.upper.CI <- NA
  }

  return(list(matched.patients=matched.patients,
              Log.Loss.for.Benefit=Log.Loss.for.Benefit,
              Log.Loss.lower.CI=Log.Loss.lower.CI,
              Log.Loss.upper.CI=Log.Loss.upper.CI,
              Brier.for.Benefit=Brier.for.Benefit,
              Brier.lower.CI=Brier.lower.CI,
              Brier.upper.CI=Brier.upper.CI))
}
