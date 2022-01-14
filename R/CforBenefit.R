#' @title C-for-benefit
#' @description This function calculates the C-for-benefit, as proposed by
#' D. van Klaveren et al. (2018), which measures the discriminative ability of
#' models predicting individualized treatment effect. The C-for-benefit
#' corresponds to the probability that from two randomly chosen matched patient
#' pairs with unequal observed treatment effect, the pair with greater observed
#' treatment effect also has a higher predicted treatment effect. Observed
#' treatment effect was defined as the difference between outcomes in pairs of
#' patients matched on patient characteristics (or on individualized treatment
#' effect predictions). Predicted treatment effect of a matched pair was defined
#' as the difference between the predicted outcome probability of the untreated
#' patient minus the predicted outcome probability of the treated patient.
#'
#' @importFrom MatchIt matchit
#' @importFrom dplyr slice
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
#' @param measure measure option of matchit function from MatchIt package (default="nearest")
#' @param distance distance option of matchit function from MatchIt package (default="mahalanobis)
#' @param ... additional arguments for matchit function from MatchIt package
#'
#' @return The output of the C.for.benefit function is a "list" with the following components.
#'
#' matched.patients
#'
#' a dataframe containing the matched patients.
#'
#'
#' c.for.benefit
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
#' n <- 100
#' Y <- sample(0:1, n, replace=TRUE)
#' W <- sample(0:1, n, replace=TRUE)
#' X <- matrix(rnorm(n), n, 3)
#' p.0 <- runif(n)
#' p.1 <- runif(n)
#' tau.hat <- runif(n)
#' CB.out <- C.for.Benefit(Y=Y, W=W, X=X, p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                         CI=TRUE, nr.bootstraps=100, message=TRUE,
#'                         measure="nearest", distance="mahalanobis")
#' CB.out
C.for.Benefit <- function(Y, W, X,
                          p.0, p.1, tau.hat,
                          CI=FALSE, nr.bootstraps=50, message=TRUE,
                          measure="nearest", distance="mahalanobis", ...){
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

  stopifnot("CI must be a boolean (TRUE or FALSE)" = isTRUE(CI)|isFALSE(CI))
  stopifnot("message must be a boolean (TRUE or FALSE)" = isTRUE(message)|isFALSE(message))

  # patient can only be matched to one other patient from other treatment arm
  if (sum(W==1) <= sum(W==0)){
    # ATT: all treated patients get matched with control patient
    estimand.meth <- "ATT"
  } else if (sum(W==1) > sum(W==0)){
    # ATC: all control patients get matched with treated patient
    estimand.meth <- "ATC"
  }

  # combine all data in one dataframe
  data.df <- data.frame(match.id=1:length(Y),
                        W=W, X=X, Y=Y,
                        p.0=p.0, p.1=p.1, tau.hat=tau.hat)

  # match on covariates
  matched <- MatchIt::matchit(W ~ X, data=data.df,
                              method=measure, distance=distance,
                              estimand=estimand.meth, ...) # TODO: add to documentations that the ... are for matchit
  matched.patients <- MatchIt::match.data(matched)
  matched.patients$subclass <- as.numeric(matched.patients$subclass)

  # sort on subclass and W
  matched.patients <- matched.patients[with(matched.patients, order(subclass, 1-W)), ]

  # matched observed treatment effect
  observed.TE <- stats::aggregate(matched.patients, list(matched.patients$subclass), diff)$Y
  matched.patients$matched.tau.obs <- rep(observed.TE, each=2)

  # matched p.0 = P[Y = 1| W = 0] so the probability of an outcome given no treatment of the untreated patient
  matched.p.0 <- (1-matched.patients$W)*matched.patients$p.0
  matched.patients$matched.p.0 <- rep(matched.p.0[matched.p.0!=0], each=2)

  # matched p.1 = P[Y = 1| W = 1] so the probability of an outcome given no treatment of the treated patient
  matched.p.1 <- matched.patients$W*matched.patients$p.1
  matched.patients$matched.p.1 <- rep(matched.p.1[matched.p.1!=0], each=2)

  # matched treatment effect
  matched.patients$matched.tau.hat <- matched.patients$matched.p.0 - matched.patients$matched.p.1

  # C-for-benefit
  cindex <- Hmisc::rcorr.cens(matched.patients$matched.tau.hat, matched.patients$matched.tau.obs)
  c.for.benefit <- cindex["C Index"][[1]]

  if (CI){
    if (message){
      cat('Calculating confidence interval... Taking too long? Lower the number of bootstraps. \n')
    }
    CB.for.CI <- c()
    for (B in 1:nr.bootstraps){
      # bootstrap matched patient pairs
      subclass.IDs <- unique(matched.patients$subclass)
      sample.subclass <- sample(subclass.IDs, length(subclass.IDs), replace=TRUE)
      dup.subclass.IDs <- c()
      for (i in sample.subclass){
        dup.subclass.IDs <- c(dup.subclass.IDs, matched.patients[matched.patients$subclass==i, 'match.id'])
      }
      duplicated.matched.patients <- slice(matched.patients, dup.subclass.IDs)

      # calculate C-for-benefit for duplicated matched pairs
      duplicated.cindex <- Hmisc::rcorr.cens(duplicated.matched.patients$matched.tau.hat, duplicated.matched.patients$matched.tau.obs)
      CB.for.CI <- c(CB.for.CI, duplicated.cindex["C Index"][[1]])
    }
    lower.CI <- as.numeric(stats::quantile(CB.for.CI, 0.025))
    upper.CI <- as.numeric(stats::quantile(CB.for.CI, 0.975))
  }
  else{
    lower.CI <- NA
    upper.CI <- NA
  }

  return(list(matched.patients=matched.patients, c.for.benefit=c.for.benefit, lower.CI=lower.CI, upper.CI=upper.CI))
}
