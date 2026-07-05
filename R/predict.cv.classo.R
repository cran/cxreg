#' Make predictions from a "cv.classo" object
#'
#' This function makes predictions from a cross-validated classo model, using
#' the stored \code{"classo.fit"} object and a chosen lambda value.
#' @param object Fitted \code{"cv.classo"} object.
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#'   made. Must be a complex-valued matrix.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#'   predictions are required. Default is \code{"lambda.1se"}. Alternatively
#'   \code{"lambda.min"} can be used, or a numeric value (or vector) of
#'   \code{lambda} directly.
#' @param \dots Additional arguments passed to the \code{predict} method for
#'   \code{"classo"} objects.
#' @return The object returned by \code{\link{predict.classo}} for the chosen
#'   \code{s}: a complex-valued matrix of predicted values, or a coefficient
#'   matrix, depending on \code{type}.
#' @author Younghoon Kim, Navonil Deb, Sumanta Basu \cr Maintainer:
#' Younghoon Kim <ykim124@ua.edu>
#' @seealso \code{\link{classo}}, \code{\link{predict.classo}},
#'   \code{\link{coef.cv.classo}}, \code{\link{cv.classo}}.
#' @keywords models regression
#' @method predict cv.classo
#' @export
predict.cv.classo <- function(object, newx,
                              s = c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s      <- match.arg(s)
    lambda <- object[[s]]
    names(lambda) <- s
  } else {
    stop("Invalid form for s")
  }
  predict(object$classo.fit, newx, s = lambda, ...)
}