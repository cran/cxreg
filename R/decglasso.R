#' Debias the cglasso estimator
#'
#' Computes the one-step debiased (desparsified) estimator of the spectral
#' precision matrix from a fitted \code{cglasso} object. The debiasing
#' correction follows the formula
#' \deqn{\tilde{\Theta} = 2\hat{\Theta} - \hat{\Theta}\hat{f}\hat{\Theta},}
#' where \eqn{\hat{\Theta}} is the \code{cglasso} estimate and \eqn{\hat{f}}
#' is the smoothed periodogram at the target frequency.
#'
#' @param object A fitted object of class \code{"cglassofit"} returned by
#'   \code{\link{cglasso}}, or a \eqn{p \times p} complex-valued matrix
#'   giving \eqn{\hat{\Theta}} directly.
#' @param fhat A \eqn{p \times p} complex-valued matrix. The smoothed spectral
#'   density (periodogram) estimate at the target Fourier frequency,
#'   typically the output of \code{fhat_at()}.
#' @param index Integer. When \code{object} is a \code{"cglassofit"} object,
#'   \code{index} selects which element of \code{object$Theta_list} to use as
#'   \eqn{\hat{\Theta}}. Defaults to \code{object$min_index}, i.e. the lambda
#'   chosen by the stopping criterion.
#' @return A list of class \code{"decglasso"} with the following components:
#' \describe{
#'   \item{\code{Theta_tilde}}{The \eqn{p \times p} complex-valued debiased
#'     spectral precision matrix estimate.}
#'   \item{\code{Theta_hat}}{The \eqn{p \times p} complex-valued \code{cglasso}
#'     estimate used as input to the debiasing step.}
#'   \item{\code{fhat}}{The smoothed spectral density matrix supplied via
#'     \code{fhat}.}
#' }
#'
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @references Deb, N., Kim Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords models complex-valued precision matrix
#'
#' @seealso \code{\link{cglasso}}, \code{\link{cglasso.path}}
#'
#' @examples
#' library(mvtnorm)
#' p <- 30
#' n <- 500
#' C <- diag(0.7, p)
#' C[row(C) == col(C) + 1] <- 0.3
#' C[row(C) == col(C) - 1] <- 0.3
#' Sigma <- solve(C)
#' set.seed(1010)
#' m <- floor(sqrt(n)); j <- 1
#' X_t <- rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
#' dft  <- dft.all(X_t)
#' fhat <- fhat_at(dft, j, m)
#' fit  <- cglasso(S = fhat, m = m)
#'
#' res <- decglasso(object = fit, fhat = fhat)
#' Theta_tilde <- res$Theta_tilde
#'
#' @export decglasso
decglasso <- function(object,
                      fhat,
                      index = NULL) {
  
  # Input validation 
  # Resolve Theta_hat from 'object'
  if (inherits(object, "cglassofit")) {
    if (is.null(index)) {
      index <- object$min_index
    }
    if (!is.numeric(index) || length(index) != 1L ||
        index < 1 || index > length(object$Theta_list)) {
      stop("'index' must be a single integer in 1:", length(object$Theta_list))
    }
    Theta_hat <- object$Theta_list[[index]]
  } else if (is.matrix(object) && is.complex(object)) {
    Theta_hat <- object
  } else {
    stop("'object' must be a 'cglassofit' object returned by cglasso(), ",
         "or a complex-valued matrix.")
  }
  
  p <- nrow(Theta_hat)
  
  # Validate fhat
  if (!is.matrix(fhat) || !is.complex(fhat)) {
    stop("'fhat' must be a complex-valued matrix.")
  }
  if (nrow(fhat) != p || ncol(fhat) != p) {
    stop("'fhat' must have the same dimensions as the precision matrix (",
         p, " x ", p, ").")
  }
  
  # One-step debiasing 
  #   Theta_tilde = 2 * Theta_hat  -  Theta_hat %*% fhat %*% Theta_hat
  Theta_tilde <- 2 * Theta_hat - Theta_hat %*% fhat %*% Theta_hat
  
  # Return 
  output <- list(
    Theta_tilde = Theta_tilde,
    Theta_hat   = Theta_hat,
    fhat        = fhat
  )
  class(output) <- "decglasso"
  return(output)
}