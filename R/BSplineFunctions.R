## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to generate B-spline basis functions
## ------------------------------------------------------------------------

# GetBasis() --------------------------------------------------------------
#' Get Basis functions
#'
#' Retrieves the \eqn{i}-th basis function mapping for some given range
#' \code{x}, degree and augmented knot set by means of recursion.
#'
#' @param x Range of values to define the function over.
#' @param deg Degree of the desired B-Spline.
#' @param KKnots Augmented knot set.
#' @param i Current level of basis function.
#'
#' @return The \eqn{i}-th basis function \eqn{\phi} of level \code{deg}
#' from B-Spline generative procedure
#'
GetBasis <- function(x, deg, KKnots, i) {
  # Return base case (degree = 0)
  if (deg == 0) {
    B <- ifelse(((x >= KKnots[i]) & (x < KKnots[i+1])), 1, 0)
    return(B)
  }

  if ((KKnots[i + deg] - KKnots[i]) != 0) {
    alpha1 <- (x - KKnots[i]) / (KKnots[i + deg] - KKnots[i])
  } else {
    alpha1 <- 0
  }

  if ((KKnots[i + deg + 1] - KKnots[i + 1]) != 0) {
    alpha2 <- (KKnots[i + deg + 1] - x) / (KKnots[i + deg + 1] - KKnots[i + 1])
  } else {
    alpha2 <- 0
  }

  phi <- (alpha1 * GetBasis(x, deg - 1, KKnots, i)) +
    (alpha2 * GetBasis(x, deg - 1, KKnots, i + 1))

  return(phi)
}


# GetBSpline() ------------------------------------------------------------
#' Get B-Spline Matrix
#'
#' Generates B-Spline functions over some parameters and places such functions
#' into columns of a \eqn{n} by \eqn{(m + deg + 1)} matrix
#'
#' @param x Range of values to define the function over.
#' @param deg Degree of the desired B-Spline.
#' @param IntKnots Interior knots that partially define the B-Spline.
#' @param  ExtKnots Exterior knots, often \code{ExtKnots = c(min(x),max(x))}.
#'
#' @return Matrix where \eqn{i,j}-th entry corresponds to \eqn{j}-th basis
#' function evaluated at \eqn{i}-th data point.
#'
#' @export
#'
GetBSpline <- function(x, deg = 3, IntKnots, ExtKnots) {
  # Augment exterior knots around interior knots
  AugKnots <- c(rep(ExtKnots[1], deg + 1), IntKnots, rep(ExtKnots[2], deg + 1))

  # Expect m+k basis functions, call this integer "NumF"
  NumF <- length(IntKnots) + (deg + 1)

  # Fill matrix columns with basis functions
  B <- matrix(0, length(x), NumF)
  for (i in 1:NumF) B[,i] <- PSplinesR::GetBasis(x, deg, AugKnots, i)

  # Manually add in boundary to final basis function
  if(any(x == ExtKnots[2])) B[x == ExtKnots[2], NumF] <- 1
  return(B)
}

# PlotBSpline() -----------------------------------------------------------
#' Plot B-Spline functions
#'
#' Takes some B-Spline basis matrix from \code{GetBSpline()} and plots the
#' output. *Note* \code{x} must match the \code{x} used in
#' \code{GetBSpline(x, ...)}.
#'
#' @param x Same \code{x} used in \code{GetBSpline(x, ...)}.
#' @param B Output from \code{GetBSpline()} call.
#' @param lty Line type. Can be single value or range of values.
#' @param col Line colour. Can be single value or range of values.
#'
#' @export
#'
PlotBSpline <- function(x, B, lty = 1:ncol(B), col = 1:ncol(B)) {
  # Possible to only pass a single graphical option
  if (length(lty == 1)) lty = rep(lty, ncol(B))
  if (length(col == 1)) col = rep(col, ncol(B))

  # Create new plot window
  plot(NULL, xlab = "x", ylab="y",
       xlim = c(min(x), max(x)), ylim = c(min(B),1.1 * max(B)))

  # Add basis functions to plot
  for (i in 1:ncol(B)) lines(x, B[, i], lty = lty[i], col = col[i])
}
