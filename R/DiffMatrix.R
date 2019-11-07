## ------------------------------------------------------------------------
## Josh Cowley. Last Modified 07-11-19
## Functions to generate difference matrix
## ------------------------------------------------------------------------

# FUN -- Obtain difference matrix
# ISSUE -- Only gets difference matrix of order 2 for now
# m <- No. of interior knots
# ord <- Order of difference matrix
GetDiffMatrix <- function(m, ord = 2) {
  # Set up difference matrix
  D <- matrix(0, nrow = m, ncol = m + ord)

  # Determine first row
  FRow <- c(c(1, -2, 1), rep(0,m + ord - 3))

  # Fill in diff. matrix
  for (nr in 1:nrow(D)) D[nr, ] <- CPerm(FRow,nr-1)

  return(D)
}

# FUN -- Move final element to start of vector (cyclic permutation)
# x <- Vector to be permuted
# i <- No. of cycles
CPerm <- function(x, i = 1) {
  if(i == 0) {
    return(x)
  } else {
    LastElement <- tail(x,1)
    CycledX <- c(LastElement, head(x,-1))
    return(CPerm(CycledX, i - 1))
  }
}
