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
  for (i in 1:NumF) B[,i] = GetBasis(x, deg, AugKnots, i)

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
