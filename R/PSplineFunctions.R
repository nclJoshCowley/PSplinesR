## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to fit P-Spline (Natural cubic B-Spline with n-2 interior knots)
## ------------------------------------------------------------------------


# FitPSpline() ------------------------------------------------------------
#' Fit Natural cubic spline with data points as knots (P-Spline)
#'
#' This will utilise the B-Spline functions in this package to fit a
#'   P-Spline. That is, a cubic B-spline with interior and exterior knots
#'   defined by each data point.
#'
#' @param x Vector of data points.
#' @param y Vector of observations
#' @param lambda Smoothing parameter.
#'
#' @return Named list with the following components:
#'
#' \code{AlphaEst} is a vector of estimated coefficients of the
#'   regression model.
#'
#' \code{YEst} is a vector of fitted values.
#'
#' \code{VarYEst} is the dispersion matrix of the fitted value estimates.
#'
#' @export
#'
FitPSpline <- function(x, y, lambda) {
  # Get natural cubic B-spline
  n <- length(x)
  IK <- x[-c(1, n)]
  EK <- x[c(1, n)]
  B <- PSplinesR::GetBSpline(x, deg = 3, IK, EK)

  # Get difference matrix (see DiffMatrix.R)
  D <- GetDiffMatrix(dim(B)[1], 2)

  # Calculate \hat{\bm{\alpha} and Hat matrix
  temp <-  solve((t(B) %*% B) + lambda * (t(D) %*% D), t(B))
  AlphaEst <-  temp %*% y
  H <- B %*% temp

  # Using theory from Molinair (2014) and LaTeX doc.
  YEst <- H %*% y
  Sig2Est <- sum((y - YEst)^2) / (n - sum(diag(H)))
  YVarEst <- Sig2Est * (H %*% t(H))

  return(list(AlphaEst = AlphaEst,
              YEst = as.vector(t(YEst)),
              YVarEst = YVarEst))
}


# PredictPSpline() --------------------------------------------------------
#' Prediction via P-Spline Model
#'
#' Predict a value for \eqn{y'}, given new value \eqn{x'} based on regression
#' model and parameter estimates found via P-Spline method.
#'
#' @inheritParams FitPSpline
#' @param XNew A single data point, for which to predict a \code{y} value.
#'   (Need to update code to accept vectors of XNew)
#'
#' @return Single predicted \code{y} value. (See note above)
#'
#' @export
#'
PredictPSpline <- function(x, y, lambda, XNew) {
  if(length(XNew) != 1) stop("XNew must be a single value")

  # Get B-spline basis functions at XNew
  IK <- x[-c(1, length(x))]
  EK <- x[c(1, length(x))]
  BNew <- PSplinesR::GetBSpline(XNew, deg = 3, IK, EK)

  # y* = sum_{j=1,..,q} b_q(x*) alpha_q
  AlphaEst <- PSplinesR::FitPSpline(x, y, lambda)$AlphaEst
  YNew = sum(BNew %*% AlphaEst)

  return(as.numeric(YNew))
}


# PlotPSplineFit() --------------------------------------------------------
#' Plot P-Spline Fit to Data
#'
#' Plot given data with curve showing fitted values as well as overlaying
#'   confidence intervals for given probability.
#'
#' @inheritParams FitPSpline
#' @param CI Given confidence level.
#' @param ... Optional graphical parameters.
#'
#' @export
#'
PlotPSplineFit <- function(x, y, lambda, CI = 0.95, ...) {
  # Plot original data
  graphics::plot(x, y, pch = 4, ...)

  # Calculate SE
  fit <- PSplinesR::FitPSpline(x, y, lambda)
  SE <- stats::qnorm(0.5 * (CI + 1)) * sqrt(diag(fit$YVarEst))

  # Calculate fitted values with CI values
  toPlot <- list(mean = fit$YEst, lower = fit$YEst - SE, upper = fit$YEst + SE)

  ## CI
  graphics::polygon(c(x, rev(x)), c(toPlot$lower, rev(toPlot$upper)),
          col = "grey75", border = F)
  # Re-add points over polygon
  graphics::points(x, y, pch = 4)
  ## Fitted line
  graphics::lines(x, toPlot$mean)
}
