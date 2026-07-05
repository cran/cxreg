#' plot coefficients from a "classo" object
#'
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{"classo"} object.
#'
#' A coefficient profile plot is produced.
#'
#' @aliases plot.classo
#' @param x fitted \code{"classo"} model
#' @param xvar What is on the X-axis. \code{"norm"} plots against the L1-norm
#' of the coefficients, \code{"lambda"} against the log-lambda sequence, and
#' \code{"dev"} against the percent deviance explained.
#' @param label If \code{TRUE}, label the curves with variable sequence
#' numbers.
#' @param \dots Other graphical parameters to plot
#' 
#' @return No return value, called for side effects (produces a plot).
#' 
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @seealso \code{classo}
#'
#' @method plot classo
#' @export
plot.classo <- function(x, xvar = c("norm", "lambda", "dev"), label = FALSE, ...) {
  
  xvar <- match.arg(xvar)
  # Pass complex beta directly; plotCoef splits into Re and Im panels internally
  plotCoef(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio,
           label = label, xvar = xvar, ...)
}