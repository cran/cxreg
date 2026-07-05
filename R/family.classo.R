#' Extract the family from a classo object
#'
#' Returns the model family for objects fitted by \code{\link{classo}} or
#' \code{\link{cv.classo}}. All \code{classo} models fit a complex-valued
#' Gaussian likelihood, so the family is always \code{"gaussian"}.
#'
#' @param object A fitted model object of class \code{"classo"},
#'   \code{"classofit"}, or \code{"cv.classo"}.
#' @param \dots Not used.
#' @return A character string: \code{"gaussian"}.
#' @seealso \code{\link{classo}}, \code{\link{cv.classo}}
#' @keywords models regression
#' @method family classo
#' @export
family.classo <- function(object, ...) {
  "gaussian"
}

#' @rdname family.classo
#' @method family classofit
#' @export
family.classofit <- function(object, ...) {
  "gaussian"
}

#' @rdname family.classo
#' @method family cv.classo
#' @export
family.cv.classo <- function(object, ...) {
  family(object$classo.fit)
}