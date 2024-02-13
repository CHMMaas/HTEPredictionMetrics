#' @title Overall performance metrics for benefit
#' @description This function calculates logistic-loss-for-benefit and
#' Brier-for-benefit.  Please note, this function is only applicable for binary
#' outcomes.
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
#' @param matched.patients dataframe; optional if you want to provide your own dataframe (including p.0, p.1, matched.tau.hat, matched.tau.obs, and also include subclass if confidence interval needs to be computed) of matched patients, otherwise patients will be matched (default=NULL)
#' @param CI boolean; TRUE compute confidence interval; default=FALSE do not compute confidence interval (default=FALSE)
#' @param nr.bootstraps boolean; number of bootstraps to use for confidence interval computation (default=1)
#' @param message boolean; TRUE display computation time message; FALSE do not display message (default=TRUE)
#' @param measure measure option of matchit function from MatchIt package (default="nearest")
#' @param distance distance option of matchit function from MatchIt package (default="mahalanobis)
#' @param estimand default ATC meaning treated units are selected to be matched with control units
#' @param ... additional arguments for matchit function from MatchIt package
#'
#' @return The output of the OP.for.Benefit function is a "list" with the following components.
#'
#' matched.patients
#'
#' a dataframe containing the matched patients.
#'
#'
#' Cross.entropy.for.benefit
#'
#' the resulting logistic-loss-for-Benefit value.
#'
#'
#' Cross.entropy.lower.CI
#'
#' the lower bound of the confidence interval of the logistic-loss-for-Benefit (if CI = TRUE).
#'
#'
#' Cross.entropy.upper.CI
#'
#' the upper bound of the confidence interval of the logistic-loss-for-Benefit (if CI = TRUE).
#'
#'
#' Brier.for.benefit
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
#' set.seed(1)
#' n <- 100
#' Y <- sample(0:1, n, replace=TRUE)
#' W <- sample(0:1, n, replace=TRUE)
#' X <- matrix(rnorm(n*3), n, 3)
#' p.0 <- runif(n)
#' p.1 <- runif(n)
#' tau.hat <- runif(n)
#' OP.out <- OP.for.Benefit(Y=Y, W=W, X=X, p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE,
#'                         matched.patients=NULL,
#'                         measure="nearest", distance="mahalanobis",
#'                         estimand=NULL)
#' OP.out
#'
#' # alternatively, use a dataframe of matched patients and calculate the overall performance metrics
#' out.matched <- match.patients(Y=Y, W=W, X=X,
#'                               p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                               print=TRUE, measure="nearest",
#'                               distance="mahalanobis", estimand=NULL)
#' OP.out <- OP.for.Benefit(matched.patients=out.matched$df.matched.pairs,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE)
#' OP.out
OP.for.Benefit <- function(Y=NULL, W=NULL, X=NULL,
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

  # prepare tau and indicator functions
  t.1 <- (1-matched.patients$matched.p.1)*matched.patients$matched.p.0
  t.0 <- (1-matched.patients$matched.p.1)*(1-matched.patients$matched.p.0) + matched.patients$matched.p.1*matched.patients$matched.p.0
  t.min1 <- matched.patients$matched.p.1*(1-matched.patients$matched.p.0)

  I.1 <- matched.patients$matched.tau.obs==1
  I.0 <- matched.patients$matched.tau.obs==0
  I.min1 <- matched.patients$matched.tau.obs==-1

  # Brier score for benefit
  n.p <- nrow(matched.patients)
  Brier.for.benefit <- (sum((t.1-I.1)^2)
                        +sum((t.0-I.0)^2)
                        +sum((t.min1-I.min1)^2))/(2*n.p)

  # Logistic loss for benefit
  omit <- which(t.1<0 | t.0<0 | t.min1<0)
  if (length(omit)>0){
    if (message){
      cat('nr. omitted observations for log loss:', length(omit), '\n')
    }
    Cross.entropy.for.benefit <- -(sum(I.1[-omit]*log(t.1[-omit]))+sum(I.0[-omit]*log(t.0[-omit]))+sum(I.min1[-omit]*log(t.min1[-omit])))/n.p
  }
  else{
    Cross.entropy.for.benefit <- -(sum(I.1*log(t.1))+sum(I.0*log(t.0))+sum(I.min1*log(t.min1)))/n.p
  }

  if (CI){
    if (message){
      cat('Calculating confidence interval... Taking too long? Lower the number of bootstraps. \n')
    }
    Cross.Entropy.for.CI <- c()
    Brier.for.CI <- c()
    # bootstrap matched patient pairs
    for (B in 1:nr.bootstraps){
      # obtain all subclass IDs
      subclass.IDs <- sort(unique(matched.patients$subclass))
      # sample randomly from subclass IDs
      sample.subclass <- sample(subclass.IDs, length(subclass.IDs), replace=TRUE)
      # data frame of duplicated patients
      matched.patients.dup <- dplyr::slice(matched.patients, sample.subclass)

      # set up quantities
      t.1.B <- (1-matched.patients.dup$matched.p.1)*matched.patients.dup$matched.p.0
      t.0.B <- (1-matched.patients.dup$matched.p.1)*(1-matched.patients.dup$matched.p.0) + matched.patients.dup$matched.p.1*matched.patients.dup$matched.p.0
      t.min1.B <- matched.patients.dup$matched.p.1*(1-matched.patients.dup$matched.p.0)

      I.1.B <- matched.patients.dup$matched.tau.obs==1
      I.0.B <- matched.patients.dup$matched.tau.obs==0
      I.min1.B <- matched.patients.dup$matched.tau.obs==-1

      # Brier score for benefit
      Brier.for.benefit.B <- (sum((t.1.B-I.1.B)^2)
                              +sum((t.0.B-I.0.B)^2)
                              +sum((t.min1.B-I.min1.B)^2))/(2*n.p)

      # Logistic loss for benefit
      omit.B <- which(t.1.B<0 | t.0.B<0 | t.min1.B<0)
      if (length(omit.B)>0){
        if (message){
          cat('nr. omitted observations for log loss:', length(omit.B), '\n')
        }
        Cross.entropy.for.benefit.B <- -(sum(I.1.B[-omit.B]*log(t.1.B[-omit.B]))
                                    +sum(I.0.B[-omit.B]*log(t.0.B[-omit.B]))
                                    +sum(I.min1.B[-omit.B]*log(t.min1[-omit.B])))/n.p
      }
      else{
        Cross.entropy.for.benefit.B <- -(sum(I.1.B*log(t.1.B))
                                  +sum(I.0.B*log(t.0.B))
                                  +sum(I.min1.B*log(t.min1.B)))/n.p
      }

      # calculate calibration metrics
      Cross.Entropy.for.CI <- c(Cross.Entropy.for.CI, Cross.entropy.for.benefit.B)
      Brier.for.CI <- c(Brier.for.CI, Brier.for.benefit.B)
    }
    Cross.entropy.lower.CI <- as.numeric(stats::quantile(Cross.Entropy.for.CI, 0.025))
    Cross.entropy.upper.CI <- as.numeric(stats::quantile(Cross.Entropy.for.CI, 0.975))
    Brier.lower.CI <- as.numeric(stats::quantile(Brier.for.CI, 0.025))
    Brier.upper.CI <- as.numeric(stats::quantile(Brier.for.CI, 0.975))
  }
  else{
    Cross.entropy.lower.CI <- NA
    Cross.entropy.upper.CI <- NA
    Brier.lower.CI <- NA
    Brier.upper.CI <- NA
  }

  return(list(matched.patients=matched.patients,
              Cross.entropy.for.benefit=Cross.entropy.for.benefit,
              Cross.entropy.lower.CI=Cross.entropy.lower.CI,
              Cross.entropy.upper.CI=Cross.entropy.upper.CI,
              Brier.for.benefit=Brier.for.benefit,
              Brier.lower.CI=Brier.lower.CI,
              Brier.upper.CI=Brier.upper.CI))
}
