## ------------------------------------------------------------------------
## Josh Cowley. Last Modified 06-11-19
## Functions to generate B-spline basis functions
## ------------------------------------------------------------------------

# FUN -- Gets the i^th basis functions for some given deg / AugKnots.
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

# FUN -- Create matrix with columns equal to B-Spline basis functions.
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

# METHOD -- Plot the B-Spline basis matrix produced by GetBSPline()
ShowBSpline <- function(x, B, lty = 1:ncol(B), col = 1:ncol(B)) {
  # Possible to only pass a single graphical option
  if (length(lty == 1)) lty = rep(lty, ncol(B))
  if (length(col == 1)) col = rep(col, ncol(B))

  # Create new plot window
  plot(NULL, xlab = "x", ylab="y",
       xlim = c(min(x), max(x)), ylim = c(min(B),1.1 * max(B)))

  # Add basis functions to plot
  for (i in 1:ncol(B)) lines(x, B[, i], lty = lty[i], col = col[i])
}

# FUN -- Estimate values given data (x,y) with dispersion matrix estimate
FittedVals <- function(xx, yy, lambda) {
  # Get natural cubic B-spline
  n <- length(xx)
  IK <- xx[-c(1, n)]
  EK <- xx[c(1, n)]
  B <- PSplinesR::GetBSpline(xx, deg = 3, IK, EK)

  # Get difference matrix (see DiffMatrix.R)
  D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2)

  # Calculate \hat{\bm{\alpha} and Hat matrix
  temp <-  solve((t(B) %*% B) + lambda * (t(D) %*% D), t(B))
  AlphaEst <-  temp %*% yy
  H <- B %*% temp

  # Using theory from Molinair (2014) and LaTeX doc.
  YEst <- H %*% yy
  Sig2Est <- sum((yy - YEst)^2) / (n - sum(diag(H)))
  YVarEst <- Sig2Est * (H %*% t(H))

  return(list(AlphaEst = AlphaEst,
              YEst = as.vector(t(YEst)),
              YVarEst = YVarEst))
}

# FUN -- Get alpha fitted values from (xx,yy) then use est's to get y* = f(x*)
PredictVals <- function(xx, yy, lambda, XNew) {
  # Get B-spline basis functions at XNew
  IK <- xx[-c(1, length(xx))]
  EK <- xx[c(1, length(xx))]
  BNew <- PSplinesR::GetBSpline(XNew, deg = 3, IK, EK)

  # y* = sum_{j=1,..,q} b_q(x*) alpha_q
  AlphaEst <- PSplinesR::FittedVals(xx, yy, lambda)$AlphaEst
  YNew = sum(BNew %*% AlphaEst)

  return(as.numeric(YNew))
}

# FUN -- Leave one out cross validation
CalculateCV <- function(xx, yy, lambda) {
  n <- length(xx)

  # Get CV error for each observation
  CVerror <- vector(mode="numeric", length = length(xx))
  for (i in 1:n) {
    TestYEst <- PredictVals(xx[-i], yy[-i], lambda, xx[i])
    CVerror[i] <- yy[i] - TestYEst
  }

  # CV(lambda) is mean of CV error sqaured
  return(mean(CVerror^2))
}

# FUN -- Returns lambda that minimises CV error
GetOptimalSmooth <- function(xx, yy, toPlot = T, Lambdas) {
  if(missing(Lambdas)) stop("You must provide a range of smoothing values")

  CVs <-  vector(mode = "numeric", length = length(Lambdas))

  for (i in 1:length(Lambdas)) {
    CVs[i] <- CalculateCV(xx, yy, Lambdas[i])
  }

  if (toPlot) {
    plot(Lambdas, CVs, type="l")
  }

  OptimalLambda = Lambdas[which(CVs == min(CVs))]
  return(OptimalLambda)
}

# FUN -- Plot some data then plot the fitted values
ShowPSplineFit <- function(xx, yy, lambda, CI = 0.95) {
  # Plot original data
  plot(xx, yy, pch = 4)

  # Calculate SE
  fit <- PSplinesR::FittedVals(xx, yy, lambda)
  SE <- qnorm(0.5 * (CI + 1)) * sqrt(diag(fit$YVarEst))

  # Calculate fitted values with CI values
  toPlot <- list(mean = fit$YEst, lower = fit$YEst - SE, upper = fit$YEst + SE)

  ## CI
  polygon(c(xx, rev(xx)), c(toPlot$lower, rev(toPlot$upper)),
          col = "grey75", border = F)
  # Re-add points over polygon
  points(xx, yy, pch = 4)
  ## Fitted line
  lines(xx, toPlot$mean)
}

# GetOptimalSmooth(EgData129$x,EgData129$y, T, seq(907.5,908.5,0.1))
