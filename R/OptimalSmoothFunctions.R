## ------------------------------------------------------------------------
## Josh Cowley.
## Functions to select best smoothing parameter
## ------------------------------------------------------------------------
#
# # FUN -- Leave one out cross validation
# CalculateCV <- function(xx, yy, lambda) {
#   n <- length(xx)
#
#   # Get CV error for each observation
#   CVerror <- vector(mode="numeric", length = length(xx))
#   for (i in 1:n) {
#     TestYEst <-   PredictPSpline(xx[-i], yy[-i], lambda, xx[i])
#     CVerror[i] <- yy[i] - TestYEst
#   }
#
#   # CV(lambda) is mean of CV error sqaured
#   return(mean(CVerror^2))
# }
#
# # FUN -- Returns lambda that minimises CV error
# GetOptimalSmooth <- function(xx, yy, toPlot = T, Lambdas) {
#   if(missing(Lambdas)) stop("You must provide a range of smoothing values")
#
#   CVs <-  vector(mode = "numeric", length = length(Lambdas))
#
#   for (i in 1:length(Lambdas)) {
#     CVs[i] <- CalculateCV(xx, yy, Lambdas[i])
#   }
#
#   if (toPlot) {
#     graphics::plot(Lambdas, CVs, type="l")
#   }
#
#   OptimalLambda = Lambdas[which(CVs == min(CVs))]
#   return(OptimalLambda)
# }
