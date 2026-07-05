#' Extract coefficients from a cv.classo object
#'
#' A convenience wrapper that extracts coefficients from the
#' \code{classo.fit} stored in a \code{"cv.classo"} object at a chosen
#' lambda value.
#'
#' @param object Fitted \code{"cv.classo"} object.
#' @param s Value(s) of the penalty parameter \code{lambda}. Default is
#'   \code{"lambda.1se"}. Alternatively \code{"lambda.min"} can be used,
#'   or a numeric value (or vector) of \code{lambda} directly.
#' @param \dots Additional arguments passed to \code{\link{coef.classo}}.
#' @return A complex-valued coefficient matrix; see \code{\link{predict.classo}}.
#' @seealso \code{\link{classo}}, \code{\link{coef.classo}},
#'   \code{\link{cv.classo}}.
#' @keywords models regression
#' @method coef cv.classo
#' @export
coef.cv.classo <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s      <- match.arg(s)
    lambda <- object[[s]]
    names(lambda) <- s     # name the lambda so downstream code can identify it
  } else {
    stop("Invalid form for s")
  }
  coef(object$classo.fit, s = lambda, ...)
}