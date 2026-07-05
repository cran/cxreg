#' Interpolate between lambda values
#'
#' Internal function used by \code{\link{predict.classo}} to linearly
#' interpolate coefficient paths at new values of \code{lambda}.
#'
#' Given the decreasing sequence \code{lambda} from the fitted model and a
#' target vector \code{s}, the function locates the two adjacent lambda values
#' that bracket each element of \code{s} and returns the indices and
#' interpolation fraction so that the caller can compute
#' \eqn{\mathrm{frac} \cdot \beta_{\mathrm{left}} + (1 - \mathrm{frac}) \cdot \beta_{\mathrm{right}}}.
#'
#' @param lambda Decreasing numeric vector of lambda values from the fitted
#'   model.
#' @param s Numeric vector of target lambda values at which predictions are
#'   required.
#' @return A list with three elements:
#' \describe{
#'   \item{\code{left}}{Integer vector of left-bracket indices into
#'     \code{lambda}.}
#'   \item{\code{right}}{Integer vector of right-bracket indices into
#'     \code{lambda}.}
#'   \item{\code{frac}}{Numeric vector of interpolation fractions in
#'     \eqn{[0, 1]}. The interpolated value is
#'     \eqn{\mathrm{frac} \cdot \beta_{\mathrm{left}} +
#'     (1 - \mathrm{frac}) \cdot \beta_{\mathrm{right}}}.}
#' }
#' @keywords internal
lambda.interp <- function(lambda, s) {
  
  if (length(lambda) == 1) {
    # Degenerate case: only one lambda value
    nums  <- length(s)
    left  <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
    
  } else {
    
    k <- length(lambda)
    
    # Normalise both lambda and s to [0, 1] (lambda[1] is largest)
    sfrac  <- (lambda[1] - s)        / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)   / (lambda[1] - lambda[k])
    
    # Clamp s to the observed range to avoid extrapolation
    sfrac[sfrac < min(lambda)] <- min(lambda)
    sfrac[sfrac > max(lambda)] <- max(lambda)
    
    # Map normalised s to a continuous position in [1, k]
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left  <- floor(coord)
    right <- ceiling(coord)
    
    # Compute interpolation fraction between the two bracketing lambdas.
    # Result is used as: frac * beta[left] + (1 - frac) * beta[right]
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    
    # At exact grid points left == right: fraction is 1 (use left only)
    sfrac[left == right] <- 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
  }
  
  list(left = left, right = right, frac = sfrac)
}