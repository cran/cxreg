#' print a cross-validated classo object
#'
#' Print a summary of the results of cross-validation for a classo model.
#'
#' A two-row table is printed showing the optimal \code{lambda} values
#' selected by the minimum CV error (\code{lambda.min}) and the 1-standard-error
#' rule (\code{lambda.1se}), along with the corresponding CV error, its standard
#' error, and the number of nonzero coefficients.
#'
#' @aliases print.cv.classo
#' @param x fitted 'cv.classo' object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' 
#' @return The matrix above is silently returned
#' 
#' @seealso \code{classo}, \code{predict} and \code{coef} methods.
#' @keywords models regression
#' @method print cv.classo
#' @export
#' @export print.cv.classo
print.cv.classo <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  cat("\nCall: ", deparse(x$call), "\n\n")
  
  optlams <- c(x$lambda.min,x$lambda.1se)
  which <- match(optlams,x$lambda)
  mat <- with(x, cbind(optlams, which, cvm[which], cvsd[which], nzero[which]))
  dimnames(mat) <- list(c("min", "1se"), c("Lambda", "Index","Measure",
                                           "SE", "Nonzero"))
  cat("Measure:", x$name,"\n\n")
  
  mat <- data.frame(mat,check.names=FALSE)
  class(mat) <- c("anova",class(mat))
  print(mat,digits=digits)
}