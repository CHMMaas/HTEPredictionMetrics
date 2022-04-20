#' @title E-for-benefit
#' @description This function calculates E-for-benefit statistics.
#'
#' @importFrom dplyr slice
#' @importFrom stats quantile
#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 .data
#'
#' @param matched.patients dataframe; dataframe of matched patients, which contains a vector of predicted treatment effect (tau.hat), predicted treatment effect of matched patients (matched.tau.hat), and observed treatment effect (matched.tau.obs) of matched patients
#' @param limits list; indicating the x-axis and y-axis limits, e.g. list(ymin=-1, ymax=1, xmin=-1, xmax=1)
#' @param plot.CI boolean; TRUE if you want to plot the confidence interval of the calibration plot of predicted versus observed treatment effect of matched patients
#' @param ... additional arguments for loess function from loess package
#'
#' @return The output of the E.for.benefit function is a "list" with the following components.
#'
#' matched.patients
#'
#' a dataframe containing the matched patients.
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
#' limits <- list(ymin=-1, ymax=1, xmin=-1, xmax=1)
#' calibration.plot(matched.patients=EB.out$matched.patients, limits=limits, plot.CI=TRUE)
calibration.plot <- function(matched.patients=NULL,
                             limits=list(ymin=-1, ymax=1, xmin=-1, xmax=1),
                             plot.CI=FALSE, ...){
  # ensure correct data types
  stopifnot("matched.patients must be a dataframe" = is.data.frame(matched.patients))
  stopifnot("CI must be a boolean (TRUE or FALSE)" = isTRUE(plot.CI)|isFALSE(plot.CI))

  # compute the smoothed calibration curve if it is not in the matched.patient dataframe
  if (is.null(matched.patients$tau.smoothed) | plot.CI){
    # perform smoothing on matched patient pairs
    loess.calibrate <- stats::loess(matched.tau.obs ~ matched.tau.hat,
                                    data=matched.patients, ...)

    if (plot.CI){
      # compute standard error if plot around LOESS needs to be computed
      loess.result <- predict(loess.calibrate,
                              newdata=matched.patients,
                              se=TRUE)
      matched.patients$tau.smoothed <- loess.result$fit
    } else {
      matched.patients$tau.smoothed <- predict(loess.calibrate, newdata=matched.patients)
    }
  }

  # omit 2.5% and 97.5% quantiles
  quantiles <- as.numeric(quantile(matched.patients$matched.tau.hat, c(0.025, 0.975)))
  included.rows <- which(matched.patients$matched.tau.hat > quantiles[1] & matched.patients$matched.tau.hat < quantiles[2])
  matched.patients <- matched.patients[included.rows, ]

  # create plot
  build.plot <- ggplot2::ggplot(data=matched.patients, ggplot2::aes(x= .data$matched.tau.hat),
                          show.legend=TRUE)                     # set data
  build.plot <- build.plot+ggplot2::theme_light(base_size=22)               # increase font size
  build.plot <- build.plot+ggplot2::geom_line(ggplot2::aes(y= .data$tau.smoothed), # plot LOESS line
                                  color="blue", size=1)
  build.plot <- build.plot+ggplot2::geom_abline(intercept=0, linetype="dashed")# 45-degree line
  build.plot <- build.plot+ggplot2::labs(x="Predicted treatment effect",
                             y="Observed treatment effect", color=" ")   # axis names

  # edit limits
  build.plot <- build.plot+ggplot2::ylim(limits$ymin, limits$ymax)
  build.plot <- build.plot+ggplot2::xlim(limits$xmin, limits$xmax)

  # plot confidence interval
  if (plot.CI){
    y.min <- matched.patients$tau.smoothed-stats::qt(0.975, loess.result$df)*loess.result$se.fit[included.rows]
    y.max <- matched.patients$tau.smoothed+stats::qt(0.975, loess.result$df)*loess.result$se.fit[included.rows]
    build.plot <- build.plot+ggplot2::geom_ribbon(ggplot2::aes(ymin=y.min, ymax=y.max), alpha=0.2)
  }

  # show plot
  methods::show(build.plot)

  return(list(matched.patients=matched.patients, build.plot=build.plot))
}
