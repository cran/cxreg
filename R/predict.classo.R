#' Make predictions from a "classo" object
#'
#' Similar to other predict methods, this function predicts fitted values,
#' coefficients and more from a fitted \code{"classo"} object.
#'
#' \code{coef(...)} is equivalent to \code{predict(type="coefficients",...)}.
#' \code{"link"} and \code{"response"} are treated identically and both return
#' the fitted complex-valued predictions \eqn{\hat{y} = X\hat{\beta}}.
#'
#' @aliases coef.classo predict.classo
#' @param object Fitted \code{"classo"} model object.
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#'   made. Must be a complex-valued matrix. Not used for
#'   \code{type = c("coefficients", "nonzero")}.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#'   predictions are required. Default is the entire sequence used to create
#'   the model.
#' @param type Type of prediction required. \code{"link"} and \code{"response"}
#'   give the fitted values \eqn{X\hat{\beta}}. \code{"coefficients"} returns
#'   the coefficient matrix at the requested values of \code{s}.
#'   \code{"nonzero"} returns a list of indices of nonzero coefficients for
#'   each value of \code{s}.
#' @param exact Relevant only when predictions are made at values of \code{s}
#'   different from those used in fitting. If \code{exact = FALSE} (default),
#'   linear interpolation is used. If \code{exact = TRUE}, the model is refit
#'   at the new \code{s} values; the original \code{x} and \code{y} (and
#'   \code{weights} if used) must be supplied as additional named arguments.
#' @param \dots Additional arguments, e.g. \code{x=} and \code{y=} when
#'   \code{exact = TRUE}.
#' @return The object returned depends on \code{type}. \code{"coefficients"}
#'   returns a matrix; \code{"nonzero"} returns a list of index vectors;
#'   \code{"link"} and \code{"response"} return a complex matrix of
#'   predictions.
#' @author Younghoon Kim, Navonil Deb, and Sumanta Basu \cr Maintainer:
#' Younghoon Kim <ykim124@ua.edu>
#' @seealso \code{classo}, and \code{print}, and \code{coef} methods, and
#' \code{cv.classo}.
#' @references Deb, N., Kuceyeski, A. and Basu, S. (2024)
#' \emph{Regularized estimation of sparse spectral precision matrices},
#' \url{https://arxiv.org/abs/2401.11128}.
#' @keywords models regression
#' @method predict classo
#' @export
#' @export predict.classo
predict.classo <- function(object, newx, s = NULL,
                           type = c("link", "response", "coefficients", "nonzero"),
                           exact = FALSE, ...) {
  
  type <- match.arg(type)
  
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE)) {
      stop("You need to supply a value for 'newx'")
    }
  }
  
  if (exact && !is.null(s)) {
    lambda <- object$lambda
    which  <- match(s, lambda, FALSE)
    if (!all(which > 0)) {
      lambda <- unique(rev(sort(c(s, lambda))))
      object <- update(object, lambda = lambda, ...)
    }
  }
  
  nbeta <- object$beta
  
  if (!is.null(s)) {
    vnames        <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda        <- object$lambda
    lamlist       <- lambda.interp(lambda, s)
    nbeta <- (nbeta[, lamlist$left,  drop = FALSE] %*% diag(lamlist$frac, length(lamlist$frac)) +
                nbeta[, lamlist$right, drop = FALSE] %*% diag(1 - lamlist$frac, length(lamlist$frac)))
    namess <- names(s)
    if (is.null(namess)) {
      namess <- paste0("s", seq_along(s))
    }
    dimnames(nbeta) <- list(vnames, namess)
  }
  
  if (type == "coefficients") {
    return(nbeta)
  }
  
  if (type == "nonzero") {
    # Only drop the first row when an intercept was fitted
    coef_mat <- if (!is.null(object$a0)) nbeta[-1, , drop = FALSE] else nbeta
    return(nonzeroCoef(coef_mat, bystep = TRUE))
  }
  
  # type == "link" or "response": compute fitted values
  if (inherits(newx, "sparseMatrix")) {
    newx <- as(newx, "CsparseMatrix")
  }
  p <- object$dim[1]
  if (is.null(dim(newx))) {
    newx <- matrix(newx, 1, byrow = TRUE)
  }
  if (ncol(newx) != p) {
    stop(paste0("The number of variables in newx must be ", p))
  }
  
  as.matrix(newx %*% nbeta)
}