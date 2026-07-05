#' Variance-covariance estimation for the debiased spectral precision matrix
#'
#' Estimates the asymptotic variance and pseudovariance of the debiased
#' (desparsified) \code{cglasso} estimator \eqn{\tilde{\Theta}} at a target
#' Fourier frequency \eqn{\omega_j}, then decomposes them into the
#' variance of the real part, variance of the imaginary part, and their
#' cross-covariance. Two estimators are available: a plug-in estimator
#' based on the theoretical sandwich formula, and a heteroscedasticity- and
#' autocorrelation-consistent (HAC) estimator.
#'
#' For a \eqn{p \times p} complex-valued estimator \eqn{\tilde{\Theta}},
#' the asymptotic distribution involves both the variance matrix
#' \eqn{V^{(ab)} = \mathrm{Var}(\tilde{\Theta}_{ab})} and the
#' pseudovariance matrix \eqn{PV^{(ab)} = \mathrm{PVar}(\tilde{\Theta}_{ab})}.
#' From these, the variances and covariance of the real and imaginary parts are
#' recovered via
#' \deqn{\mathrm{Var}(\mathrm{Re}\,\tilde{\Theta}_{ab}) = \frac{1}{2}\mathrm{Re}(V + PV),}
#' \deqn{\mathrm{Var}(\mathrm{Im}\,\tilde{\Theta}_{ab}) = \frac{1}{2}\mathrm{Re}(V - PV),}
#' \deqn{\mathrm{Cov}(\mathrm{Re}\,\tilde{\Theta}_{ab},\,\mathrm{Im}\,\tilde{\Theta}_{ab})
#'   = \frac{1}{2}\mathrm{Im}(PV).}
#'
#' @param Theta A \eqn{nvar \times nvar} complex-valued matrix. The debiased spectral
#'   precision matrix estimate \eqn{\tilde{\Theta}} at frequency \eqn{\omega_j},
#'   typically the \code{Theta_tilde} component returned by \code{\link{decglasso}}.
#' @param X A numeric matrix of size \eqn{nobs \times nvar} (time points by
#'   variables). The DFT is computed internally via \code{\link{dft.all}}.
#' @param j An integer frequency index (1-based) identifying the target Fourier
#'   frequency \eqn{\omega_j = 2\pi j / nobs}.
#' @param m A positive integer. The half-bandwidth used in the smoothed
#'   periodogram estimator; must satisfy \eqn{2m + 1 < nobs}.
#' @param type Character string selecting the variance estimator. Either
#'   \code{"plug-in"} (default) for the sandwich plug-in estimator, or
#'   \code{"HAC"} for the heteroscedasticity- and autocorrelation-consistent
#'   estimator.
#'
#' @return A list of class \code{"varcov"} with the following components:
#' \describe{
#'   \item{\code{Vre}}{A \eqn{nvar \times nvar} real matrix. Estimated variance of
#'     the real part of \eqn{\tilde{\Theta}}.}
#'   \item{\code{Vim}}{A \eqn{nvar \times nvar} real matrix. Estimated variance of
#'     the imaginary part of \eqn{\tilde{\Theta}}.}
#'   \item{\code{Cov}}{A \eqn{nvar \times nvar} real matrix. Estimated covariance
#'     between the real and imaginary parts of \eqn{\tilde{\Theta}}.}
#'   \item{\code{V}}{A \eqn{nvar \times nvar} complex matrix. The raw variance
#'     matrix before decomposition.}
#'   \item{\code{PV}}{A \eqn{nvar \times nvar} complex matrix. The raw
#'     pseudovariance matrix before decomposition.}
#'   \item{\code{type}}{Character string indicating which estimator was used.}
#' }
#'
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @references Deb, N., Kim Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords models complex-valued precision matrix
#' @examples
#' library(mvtnorm)
#' p <- 10
#' n <- 200
#' set.seed(42)
#' X <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))
#' m <- floor(sqrt(n)); j <- 1
#' dft  <- dft.all(X)
#' fhat <- fhat_at(dft, j, m)
#' fit  <- cglasso(S = fhat, m = m)
#' res  <- decglasso(object = fit, fhat = fhat)
#'
#' # Plug-in variance estimator:
#' vc_plug <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m,
#'                    type = "plug-in")
#'
#' # HAC variance estimator:
#' vc_hac  <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m,
#'                    type = "HAC")
#' @export var.cov
var.cov <- function(Theta,
                    X,
                    j,
                    m,
                    type = c("plug-in", "HAC")) {
  
  type <- match.arg(type, choices = c("plug-in", "HAC"))
  
  # Input validation
  if (!is.matrix(Theta) || !is.complex(Theta)) {
    stop("'Theta' must be a complex-valued matrix.")
  }
  p <- nrow(Theta)
  if (ncol(Theta) != p) {
    stop("'Theta' must be square.")
  }
  if (!is.matrix(X)) {
    stop("'X' must be a matrix.")
  }
  n <- nrow(X)
  if (ncol(X) != p) {
    stop("'X' must have the same number of columns as the dimension of 'Theta'.")
  }
  if (!is.numeric(j) || length(j) != 1L || j < 1 || j > n) {
    stop("'j' must be a single integer in 1, ..., n.")
  }
  if (!is.numeric(m) || length(m) != 1L || m < 1 || m != round(m)) {
    stop("'m' must be a positive integer.")
  }
  if (2 * m + 1 >= n) {
    stop("Bandwidth 2*m+1 must be smaller than n.")
  }
  
  # Compute DFT and dispatch to the chosen estimator
  dft <- dft.all(X)
  
  if (type == "plug-in") {
    est <- var_plug(Theta = Theta, dft = dft, j = j, n = n, m = m)
  } else {
    est <- var_hac(Theta = Theta, dft = dft, j = j, n = n, m = m)
  }
  
  # Decompose into real/imaginary parts
  dcmp <- var_re_im(est$V, est$PV)
  
  output <- list(
    Vre  = dcmp$Vre,
    Vim  = dcmp$Vim,
    Cov  = dcmp$Cov,
    V    = est$V,
    PV   = est$PV,
    type = type
  )
  class(output) <- "varcov"
  return(output)
}

####################################################################
# Decompose variance/pseudovariance into real-imaginary covariances.
var_re_im <- function(V, PV) {
  Vre <- 0.5 * Re(V + PV)
  Vim <- 0.5 * Re(V - PV)
  Cov <- 0.5 * Im(PV)
  return(list(Vre = Vre, Vim = Vim, Cov = Cov))
}

####################################################################
# Plug-in sandwich variance estimator.
# Returns a list with V, PV, Phi_array, Xi_array, k_grid, W_B.
var_plug <- function(Theta, dft, j, n, m) {
  p  <- nrow(Theta)
  bw <- 2*m + 1
  k_grid <- (j + (-m:m)) %% n + 1
  
  Phi_array <- array(0 + 0i, dim = c(p, p, bw))
  Xi_array  <- array(0 + 0i, dim = c(p, p, bw))
  
  for (r in seq_len(bw)) {
    f_k <- fhat_at(dft, k_grid[r], m)
    
    ## Phi_{j,k}^{(a,b)} = (Theta^* f_k Theta)_{a,b}
    Phi_array[, , r] <- Theta %*% f_k %*% Theta
    
    ## Xi_{j,k}^{(a,b)} = (Theta^T f_k Theta)_{a,b}
    Xi_array[, , r] <- t(Theta) %*% f_k %*% Theta
  }
  W_B <- ((outer(k_grid, k_grid, "+") %% n) == 0)
  
  V  <- matrix(0 + 0i, p, p)
  PV <- matrix(0 + 0i, p, p)
  
  for (a in seq_len(p)) {
    for (b in seq_len(p)) {
      
      ## sigma term 1: sum_k Phi_k^{aa} Phi_k^{bb}
      V_A <- sum(Phi_array[a, a, ] * Phi_array[b, b, ])
      
      ## sigma term 2: sum_{k1,k2} Xi_{k1}^{ab} Conj(Xi_{k2}^{ba}) I(k1+k2 in nZ)
      V_B <- 0 + 0i
      
      ## delta term 2: sum_{k1,k2} Conj(Xi_{k1}^{aa}) Xi_{k2}^{bb} I(k1+k2 in nZ)
      PV_B <- 0 + 0i
      
      for (r1 in seq_len(bw)) {
        for (r2 in seq_len(bw)) {
          if (W_B[r1, r2]) {
            V_B <- V_B +
              Xi_array[a, b, r1] * Conj(Xi_array[b, a, r2])
            
            PV_B <- PV_B +
              Conj(Xi_array[a, a, r1]) * Xi_array[b, b, r2]
          }
        }
      }
      ## delta term 1: sum_k (Phi_k^{ab})^2
      PV_A <- sum(Phi_array[a, b, ]^2)
      
      V[a, b]  <- (V_A + V_B) / bw
      PV[a, b] <- (PV_A + PV_B) / bw
    }
  }
  
  list(V = V, PV = PV, Phi_array = Phi_array,
       Xi_array = Xi_array, k_grid = k_grid, W_B = W_B)
}

####################################################################
# HAC variance estimator.
# Returns a list with V, PV, k_grid, W.
var_hac <- function(Theta, dft, j, n, m) {
  p  <- nrow(Theta)
  bw <- 2*m + 1
  h  <- -m:m
  
  W <- diag(1, bw)
  k_grid <- (j + h) %% n
  W <- W + (outer(k_grid, k_grid, "+") %% n == 0) * 1
  
  # Compute cached Y-matrix: each row r is vec(Y_r), where
  # Y_r = Conj(dr_k) %o% (dr_k^T Theta) - Theta f_k Theta
  Y_mat <- matrix(0+0i, bw, p*p)
  for (r in seq_len(bw)) {
    dr_k <- dft[((k_grid[r] %% n) + 1), , drop = TRUE]
    Yr_vec <- as.vector(t(Conj(dr_k)) %*% Theta)
    fr_k <- fhat_at(dft, k_grid[r], m)
    Br_k <- Theta %*% fr_k %*% Theta
    Yr <- outer(Conj(Yr_vec), Yr_vec, "*") - Br_k
    Y_mat[r, ] <- as.vector(Yr)
  }
  
  # V[c]  = sum_r Y_mat[r,c] * (W %*% Conj(Y_mat))[r,c]  / bw
  Z_sig <- W %*% Conj(Y_mat)
  sigma_vec <- colSums(Y_mat * Z_sig) / bw
  
  # PV[c] = sum_r Y_mat[r,c] * (W %*% Y_mat)[r,c]  / bw
  Z_del <- W %*% Y_mat
  delta_vec <- colSums(Y_mat * Z_del) / bw
  
  list(
    V      = matrix(sigma_vec, p, p),
    PV     = matrix(delta_vec, p, p),
    k_grid = k_grid,
    W      = W
  )
}