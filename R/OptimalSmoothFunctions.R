## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to select best smoothing parameter
## ------------------------------------------------------------------------


# CalculateLOOCV() --------------------------------------------------------
#' Calculate Leave-One-Out Cross Validation
#'
#' Returns the value of the leave-one-out cross validation score for a given
#'   data set and smoothing parameter when a P-Spline fit is used.
#'
#' @inheritParams PredictPSpline
#'
#' @return Single numeric value representing the mean sqaure CV error.
#'
#' @export
#'
CalculateLOOCV <- function(x, y, lambda) {
  n <- length(x)

  # Get CV error for each observation
  CVerror <- vector(mode="numeric", length = length(x))
  for (i in 1:n) {
    TestYEst <-   PredictPSpline(x[-i], y[-i], lambda, x[i])
    CVerror[i] <- y[i] - TestYEst
  }

  # CV(lambda) is mean of CV error sqaured
  return(mean(CVerror^2))
}


# MinimiseLOOCV() ---------------------------------------------------------
#' Minimised Leave-One-Out Cross Validation
#'
#' Finds the value from a given range that minimised the mean sqaure cross
#'   validation error.
#'
#' @inheritParams CalculateLOOCV
#' @param Lambdas Vector of smoothing paramaeters to test.
#' @param toPlot Logical. Should the plot of CV error against
#'   \code{Lambdas} appear?
#'
#' @export
#'
MinimiseLOOCV <- function(x, y, Lambdas, toPlot = T) {
  if(missing(Lambdas)) stop("You must provide a range of smoothing values")

  CVs <-  vector(mode = "numeric", length = length(Lambdas))

  for (i in 1:length(Lambdas)) {
    CVs[i] <- CalculateLOOCV(x, y, Lambdas[i])
  }

  if (toPlot) {
    graphics::plot(Lambdas, CVs, type="l")
  }

  OptimalLambda = Lambdas[which(CVs == min(CVs))]
  return(OptimalLambda)
}
