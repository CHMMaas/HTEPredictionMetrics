#' @title Match patients
#' @description This function matches patient pairs using the MatchIt package.
#' Observed treatment effect was defined as the difference between outcomes in
#' pairs of patients matched on patient characteristics (or on individualized
#' treatment effect predictions). Predicted treatment effect of a matched pair
#' was defined as the difference between the predicted outcome probability of
#' the untreated patient minus the predicted outcome probability of the treated
#' patient.  Please note, this function is only applicable for binary outcomes.
#'
#' @importFrom MatchIt matchit
#' @importFrom stats aggregate
#' @importFrom dplyr setdiff
#'
#' @param Y a vector of binary outcomes; 1 if an (unfavourable) event; 0 if not
#' @param W a vector of treatment assignment; 1 for active treatment; 0 for control
#' @param X a matrix of patient characteristics or individualized treatment effect predictions, do not include Y or W in this matrix
#' @param p.0 a vector of outcome probabilities under control
#' @param p.1 a vector of outcome probabilities under active treatment
#' @param tau.hat a vector of individualized treatment effect predictions
#' @param measure measure option of matchit function from MatchIt package (default="nearest")
#' @param distance distance option of matchit function from MatchIt package (default="mahalanobis)
#' @param estimand default NULL meaning that all patients in the smallest treatment arm get matched with one patient in other treatment arm; estimand can also be "ATT", "ATC" or "ATE", see Details of matchit function of MatchIt for more information
#' @param replace boolean; TRUE if matching with replacement, FALSE if matching without replacement
#' @param ... additional arguments for matchit function from MatchIt package
#'
#' @return The output of the match.patients function is
#'
#' df.matched.patients
#'
#' a dataframe containing the matched patients, thus each individual is included in this dataframe
#'
#'
#' df.matched.pairs
#'
#' a dataframe containing only the matched pairs, and not the individuals
#'
#'
#' discarded
#'
#' a vector of ID's of the patients omitted that are not matched
#'
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
#' matched.patients <- match.patients(Y=Y, W=W, X=X,
#'                                    p.0=p.0, p.1=p.1, tau.hat=tau.hat,
#'                                    CI=FALSE, message=TRUE,
#'                                    measure="nearest", distance="mahalanobis",
#'                                    estimand=NULL, replace=FALSE)
#' matched.patients
match.patients <- function(Y, W, X,
                          p.0, p.1, tau.hat,
                          measure="nearest", distance="mahalanobis",
                          estimand=NULL, replace=FALSE, ...){
  # ensure correct data types
  stopifnot("Y must be numeric" = is.numeric(Y))
  stopifnot("W must be numeric" = is.numeric(W))
  stopifnot("X must be numeric" = is.numeric(as.matrix(X)))
  stopifnot("p.0 must be numeric" = is.numeric(p.0))
  stopifnot("p.1 must be numeric" = is.numeric(p.1))
  stopifnot("tau.hat must be numeric" = is.numeric(tau.hat))
  stopifnot("estimand must be ATT, ATC, ATE, or NULL" = is.null(estimand)|estimand=="ATT"|estimand=="ATC"|estimand=="ATE")

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

  # combine all data in one dataframe
  data.df <- data.frame(match.id=1:length(Y),
                        W=W, X=X, Y=Y,
                        p.0=p.0, p.1=p.1, tau.hat=tau.hat)

  # match on covariates
  if (is.null(estimand)){
    if (sum(W==1) <= sum(W==0)){
      # ATT: all control patients get matched with treated patient
      estimand <- "ATT"
    } else if (sum(W==1) > sum(W==0)){
      # ATC: all treated patients get matched with control patient
      estimand <- "ATC"
    }
  }
  matched <- MatchIt::matchit(W ~ X, data=data.df,
                              method=measure, distance=distance,
                              estimand=estimand, replace=replace, ...)
  matched.patients <- MatchIt::match.data(matched)
  matched.patients$subclass <- as.numeric(matched.patients$subclass)

  # patient IDs of those who weren't matched
  discarded <- dplyr::setdiff(data.df$match.id, MatchIt::get_matches(matched, data=data.df)$id)

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

  # set up df to calculate OP
  matched.patients.undup <- matched.patients[, c("subclass", "matched.tau.obs", "matched.p.0", "matched.p.1", "matched.tau.hat")]
  # remove duplicates from matched.patients data frame
  matched.patients.undup <- matched.patients.undup[rep(c(TRUE, FALSE), nrow(matched.patients.undup)/2),]

  return(list(matched.out=matched,
              df.matched.patients=matched.patients,
              df.matched.pairs=matched.patients.undup,
              discarded=discarded))
}
