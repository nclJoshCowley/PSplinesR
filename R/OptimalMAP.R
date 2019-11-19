## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to select best smoothing parameter via MAP
## ------------------------------------------------------------------------

# postDensity() ------------------------------------------------------------
#' Posterior Density of Lambda
#'
#' Uses the equation given in "Efficient and automatic methods for flexible
#'   regression on spatiotemporal data, with applications to groundwater
#'   monitoring" by Evers et al. (2015).
#'
#' @param lambda Parameter of focus.
#' @param B Design matrix obtained via natural cubic spline.
#' @param D Difference matrix used in penalty term approximation.
#' @param a,b Prior information related to normal-inverse gamma distribution.
#' @param y Vector of observations.
#' @param PriorDens Prior density of penalisation parameter.
#'
#' @return Evaluation of posterior density function given in Evers et al.
#'
#' @export
#'
postDensity <- function(lambda, B, D, a, b, y, PriorDens) {
  # Use improper flat prior if prior density not given
  if (missing(PriorDens)) PriorDens <- function(x) 1

  # Sample size
  n = length(y)

  # Gram-determinants of the two matrices
  BGram <- t(B) %*% B
  DGram <- t(D) %*% D

  # Use of the Hat matrix
  HatMatr <- B %*% solve(BGram + (lambda * DGram), t(B))

  # Ans = Term1 * (Term2a/Term2b) * PriorDens
  Term1 <- lambda ^ (0.5 * Matrix::rankMatrix(DGram))
  Term2a <- det(BGram + (lambda * DGram))^(-0.5)
  Term2b <- ((2 * b) + (t(y) %*% (diag(n) - HatMatr) %*% y))^(a + 0.5 * n)

  return(Term1 * (Term2a / Term2b) * PriorDens(lambda))
}


# PlotPostDensity() -------------------------------------------------------
#' Plot Posterior Density for Smoothing Parameter in PSpline
#'
#' Plots the posterior density (to be used with MAP) when a natural cubic
#'   spline and second order differencing matrix approximation is used.
#'
#' @inheritParams postDensity
#' @param x Data points to be used as knots.
#' @param Lambdas Range of smoothing parameters to plot / scan.
#' @param add Logical. If true, plot is added to currently open graphical
#'   space.
#' @param ... Optional graphical parameters.
#'
#' @return Value of smoothing parameter that reuturned minimal density.
#'
#' @export
#'
PlotPostDensity <- function(x, y, a, b, Lambdas, PriorDens, add = F, ...) {
  # Matrices currently defined by PSpline method
  n <- length(x)
  B <- PSplinesR::GetBSpline(x, IntKnots = x[-c(1,n)], ExtKnots = c(x[1],x[n]))
  D <- PSplinesR::GetDiffMatrix(n)

  # Calculate posterior density for given arg Lambdas
  Dens <- vector(mode = "numeric", length = length(Lambdas))
  for (i in 1:length(Lambdas)) {
    Dens[i] <- postDensity(Lambdas[i], B, D, a, b, y, PriorDens)
  }

  # IF add is FALSE, create new plot space
  if (!add) {
    graphics::plot(Lambdas, Dens, type="l", ...)
  } else {
    graphics::lines(Lambdas, Dens, ...)
  }

  # Return Lambda value that minimised density
  MAPLambda <- Lambdas[which(Dens == max(Dens))]
  return(MAPLambda)
}

