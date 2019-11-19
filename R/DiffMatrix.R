## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to generate difference matrix
## ------------------------------------------------------------------------

# GetDiffMatrix() ---------------------------------------------------------
#' Get Difference Matrix
#'
#' Calculates the difference matrix of order 2 (ISSUE -- Only gets difference
#'   matrix of order 2 for now).
#'
#' @param nr Number of rows the difference matrix will have.
#' @param ord Order of the difference matrix.
#'
#' @return A difference matrix of order \code{ord} with dimensions \code{nr}
#'   by \eqn{\code{nr} + \code{ord}}
#'
#' @export
#'
GetDiffMatrix <- function(nr, ord = 2) {
  # Set up difference matrix
  D <- matrix(0, nrow = nr, ncol = nr + ord)

  # Determine first row
  FRow <- c(c(1, -2, 1), rep(0, nr + ord - 3))

  # Fill in diff. matrix
  for (i in 1:nr) D[i, ] <- CPerm(FRow, i - 1)

  return(D)
}

# CPerm() -----------------------------------------------------------------
#' Cyclic Permuatation
#'
#' Get cyclic permutation of vector. Can can the function multiple times
#'   through the use of the \code{i} argument.
#'
#' @param x Vector to be cycled through
#' @param i Number of cyclic permutations
#'
#' @return Vector with all elements moved \code{i} spaces to the left with
#'   end elements wrapping around.
#'
#' @export
#'
CPerm <- function(x, i = 1) {
  if(i == 0) {
    return(x)
  } else {
    LastElement <- utils::tail(x,1)
    CycledX <- c(LastElement, utils::head(x,-1))
    return(CPerm(CycledX, i - 1))
  }
}
