#' Draw error bars on a plot
#'
#' Internal helper used by \code{\link{plot.cv.classo}} to draw vertical error
#' bars with horizontal caps at each lambda value on the CV curve.
#'
#' @param x Numeric vector of x-axis positions (e.g. \code{log(lambda)}).
#' @param upper Numeric vector of upper bar endpoints (e.g. \code{cvm + cvsd}).
#' @param lower Numeric vector of lower bar endpoints (e.g. \code{cvm - cvsd}).
#' @param width Relative width of the horizontal end-caps as a fraction of the
#'   x-axis range. Default \code{0.02}.
#' @param \dots Additional graphical arguments passed to \code{segments},
#'   e.g. \code{col} or \code{lwd}.
#' @return \code{NULL}, invisibly. Called for its side effect of adding segments
#'   to the current plot.
#' @keywords internal
error.bars <- function(x, upper, lower, width = 0.02, ...) {
  
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  invisible(NULL)
}