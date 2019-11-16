## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to select best smoothing parameter via LOOCV
## ------------------------------------------------------------------------

#' Posterior Density of Lambda
#'
#' Uses the equation given in "Efficient and automatic methods for flexible
#'   regression on spatiotemporal data, with applications to groundwater
#'   monitoring" by Evers et al. (2015).
#'
#' @param lambda Parameter of focus.
#' @param B Design matrix obtained via natural cubic spline.
#' @param D Difference matrix used in penalty term approximation.
#' @param a,b Prior information.
#' @param y Vector of observations.
#' @param PriorDens Prior density of penalisation parameter.
#'
#' @export
#'
postDensity <- function(lambda, B, D, a, b, y, PriorDens) {
  # Temporary -- PriorDens <- 1
  if (missing(PriorDens)) PriorDens <- function(x) dlnorm(x)

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


PlotPostDensity <- function(x, y, Lambdas) {
  n <- length(x)
  B <- PSplinesR::GetBSpline(x, IntKnots = x[-c(1,n)], ExtKnots = c(x[1],x[n]))
  D <- PSplinesR::GetDiffMatrix(n)
  a <- 0.0001
  b <- 0.0001

  Dens <- vector(mode = "numeric", length = length(Lambdas))
  for (i in 1:length(Lambdas)) Dens[i] <- postDensity(Lambdas[i], B, D, a, b, y)
  plot(Lambdas, Dens, type="l")

  MAPLambda <- Lambdas[which(Dens == max(Dens))]
  return(MAPLambda)
}
