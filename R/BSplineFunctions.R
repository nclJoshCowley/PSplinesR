## ------------------------------------------------------------------------
## Josh Cowley. Last Modified 06-11-19
## Functions to generate B-spline basis functions
## ------------------------------------------------------------------------

# FUN -- Gets the i^th basis functions for some given deg / AugKnots.
GetBasis <- function(x, deg, KKnots, i) {
  # Return base case (degree = 0)
  if(deg == 0) {
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

  phi <- (alpha1 * GetBasis(x,deg-1,KKnots,i)) + ( alpha2 * GetBasis(x,deg-1,KKnots,i+1) )

  return(phi)
}

# FUN -- Create matrix with columns equal to B-Spline basis functions.
GetBSpline <- function(x, deg = 3, IntKnots, ExtKnots) {
  # Augment exterior knots around interior knots
  AugKnots <- c(rep(ExtKnots[1], deg + 1), IntKnots, rep(ExtKnots[2], deg + 1))

  # Expect m+k basis functions, call this integer "NumF"
  NumF = length(IntKnots) + (deg+1)

  # Fill matrix columns with basis functions
  B <- matrix(0, length(x), NumF)
  for (i in 1:NumF) B[,i] = PSplinesR::GetBasis(x, deg, AugKnots, i)

  # Manually add in boundary to final basis function
  if(any(x == ExtKnots[2])) B[x == ExtKnots[2], NumF] <- 1

  return(B)
}

# METHOD -- Plot the B-Spline basis matrix produced by GetBSPline()
ShowBSpline <- function(x, B, lty=1:ncol(B), col=1:ncol(B)) {
  # Possible to only pass a single graphical option
  if (length(lty == 1)) lty = rep(lty, ncol(B))
  if (length(col == 1)) col = rep(col, ncol(B))

  # Create new plot window
  plot(NULL, xlab="x", ylab="y", xlim = c(min(x),max(x)), ylim=c(min(B),1.1 * max(B)))

  # Add basis functions to plot
  for (i in 1:ncol(B)) lines(x, B[, i], lty=lty[i], col=col[i])
}

# FUN -- Estimate values given data (x,y) with dispersion matrix estimate
FittedVals <- function(xx, yy, lambda) {
  # Get natural cubic B-spline
  n <- length(xx)
  IK <- xx[-c(1,n)]
  EK <- xx[c(1,n)]
  B <- PSplinesR::GetBSpline(xx, deg=3, IK, EK)

  # Get difference matrix (see DiffMatrix.R)
  D <- PSplinesR::GetDiffMatrix(dim(B)[1],2)

  # Calculate Hat matrix
  H <- B %*% solve((t(B) %*% B) + lambda * (t(D) %*% D), t(B))

  # Using theory from Molinair (2014) and LaTeX doc.
  yEst <- H %*% yy
  sigEst <- sum((yy - -yEst)^2) / (n - sum(diag(H)))
  yVarEst <- sigEst * (H %*% t(H))

  return(list(yEst = as.vector(t(yEst)), yVarEst = yVarEst))
}

# FUN -- Leave one out cross validation
CalculateCV <- function(xx, yy, lambda, yEst = NULL) {
  n = length(xx)

  # Get fitted values if not passed
  if(is.null(yEst)) yEst <- PSplinesR::FittedVals(xx, yy, lambda)$yEst

  # Get CV error for each observation
  CVerror <- vector(mode="numeric", length = length(xx))
  for (i in 1:n) CVerror = yy[i] - yEst[i]

  # CV(lambda) is mean of CV error sqaured
  return(mean(CVerror^2))
}
