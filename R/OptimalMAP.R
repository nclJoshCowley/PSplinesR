## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to select best smoothing parameter via LOOCV
## ------------------------------------------------------------------------

#' Title
#'
#' Description.
#'
#' @param
#'
#' @export
#'
ff <- function() {
  a <- b <- 0.0001

  postDensity <- function(lambda, D, B, a, b, y, n =length(y)) {
    temp1 <- lambda ^ (0.5 * Matrix::rankMatrix(t(D) %*% D))
    temp2Num <- norm(t(B))
  }
}
