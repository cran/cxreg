#' Confidence regions and test statistics for spectral precision matrix entries
#'
#' Computes entry-wise Z-statistics, half-widths of marginal confidence
#' intervals, Mahalanobis chi-squared statistics, and joint elliptical
#' confidence region areas for the real and imaginary parts of the debiased
#' spectral precision matrix estimator \eqn{\tilde{\Theta}}.
#'
#' For each entry \eqn{(a,b)}, the joint asymptotic distribution of
#' \eqn{(\mathrm{Re}\,\tilde{\Theta}_{ab},\, \mathrm{Im}\,\tilde{\Theta}_{ab})}
#' is bivariate normal with covariance \eqn{\Sigma_{ab}/\mathrm{bw}}, where
#' \eqn{\Sigma_{ab}} is the \eqn{2 \times 2} matrix built from the
#' \code{Vre}, \code{Vim}, and \code{Cov} components of \code{varcov}. Diagonal
#' entries are treated as real-valued (one degree of freedom); off-diagonal
#' entries use a 2-df chi-squared statistic. The default null hypothesis is
#' \eqn{H_0 : \Theta_{ab} = 0} (i.e. \code{Truth = 0}), which corresponds to
#' hypothesis testing. A non-zero \code{Truth} can be supplied for coverage
#' evaluation in simulation studies.
#'
#' @param Est A \eqn{nvar \times nvar} complex-valued matrix. The debiased spectral
#'   precision matrix estimate \eqn{\tilde{\Theta}}, typically the
#'   \code{Theta_tilde} component of a \code{\link{decglasso}} result.
#' @param varcov A \code{"varcov"} object returned by \code{\link{var.cov}},
#'   supplying the estimated \code{Vre}, \code{Vim}, and \code{Cov} matrices.
#' @param m A positive integer. The half-bandwidth used to compute
#'   \eqn{\hat{f}}; the full bandwidth \eqn{\mathrm{bw} = 2m+1} is used to
#'   scale the asymptotic covariance.
#' @param alpha A number in \eqn{(0,1)}. Nominal significance level for
#'   confidence intervals and confidence regions.
#' @param Truth A \eqn{nvar \times nvar} complex-valued matrix giving the true
#'   (or hypothesised) value of \eqn{\Theta}. Default is
#'   \code{matrix(0 + 0i, nvar, nvar)}, which corresponds to testing
#'   \eqn{H_0 : \Theta_{ab} = 0} for all entries.
#' @param eps A small positive number used to regularise near-zero variances.
#'   Default \code{1e-8}.
#'
#' @return A list of class \code{"spectest"} with the following
#'   \eqn{nvar \times nvar} matrices:
#' \describe{
#'   \item{\code{Chi_sq}}{Mahalanobis chi-squared test statistics. For
#'     off-diagonal entries these are asymptotically \eqn{\chi^2(2)} under
#'     the null; for diagonal entries \eqn{\chi^2(1)}.}
#'   \item{\code{area}}{Area of the joint \eqn{(1-\alpha)} confidence ellipse
#'     for \eqn{(\mathrm{Re}\,\tilde{\Theta}_{ab},\,\mathrm{Im}\,\tilde{\Theta}_{ab})}.}
#'   \item{\code{Z_re}}{Marginal Z-statistics for the real parts.}
#'   \item{\code{wing_re}}{Half-widths of marginal \eqn{(1-\alpha)} confidence
#'     intervals for the real parts.}
#'   \item{\code{Z_im}}{Marginal Z-statistics for the imaginary parts.}
#'   \item{\code{wing_im}}{Half-widths of marginal \eqn{(1-\alpha)} confidence
#'     intervals for the imaginary parts.}
#' }
#'
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @references Deb, N., Kim Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords models complex-valued precision matrix hypothesis test
#' @examples
#' library(mvtnorm)
#' p <- 10; n <- 200
#' set.seed(42)
#' X    <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))
#' m    <- floor(sqrt(n)); j <- 1; alpha <- 0.05
#' dft  <- dft.all(X)
#' fhat <- fhat_at(dft, j, m)
#' fit  <- cglasso(S = fhat, m = m)
#' res  <- decglasso(object = fit, fhat = fhat)
#' vc   <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m,
#'                 type = "plug-in")
#'
#' # Hypothesis test (H0: Theta = 0):
#' st <- spec.test(Est = res$Theta_tilde, varcov = vc, m = m, alpha = alpha)
#' st$Chi_sq[1:3, 1:3]
#' @export spec.test
spec.test <- function(Est,
                      varcov,
                      m,
                      alpha,
                      Truth = NULL,
                      eps   = 1e-8) {
  
  # Input validation
  if (!is.matrix(Est) || !is.complex(Est)) {
    stop("'Est' must be a complex-valued matrix.")
  }
  p <- nrow(Est)
  if (ncol(Est) != p) {
    stop("'Est' must be square.")
  }
  if (!inherits(varcov, "varcov")) {
    stop("'varcov' must be a 'varcov' object returned by var.cov().")
  }
  if (!is.numeric(m) || length(m) != 1L || m < 1 || m != round(m)) {
    stop("'m' must be a positive integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a number in (0, 1).")
  }
  if (is.null(Truth)) {
    Truth <- matrix(0 + 0i, p, p)
  }
  if (!is.matrix(Truth) || nrow(Truth) != p || ncol(Truth) != p) {
    stop("'Truth' must be a ", p, " x ", p, " matrix.")
  }
  
  bw <- 2 * m + 1
  
  # Dispatch to internal function
  result <- stat_cmp(Est      = Est,
                     Truth    = Truth,
                     var_dcmp = varcov,
                     bw       = bw,
                     alpha    = alpha,
                     eps      = eps)
  class(result) <- "spectest"
  return(result)
}

#' FDR-controlled hypothesis testing for spectral precision matrix entries
#'
#' Applies an asymptotic FDR control procedure to the matrix of chi-squared
#' test statistics produced by \code{\link{spec.test}}, returning the
#' rejection threshold, the binary decision matrix, and (optionally) the
#' empirical FDR and power when the true support is known.
#'
#' The threshold \eqn{\hat{\tau}} is the infimum over a grid of values
#' \eqn{t \in (0, 2\log q]} of
#' \deqn{\frac{q \, e^{-t/2}}{\max(1,\, |\{T_{ab} \geq t\}|)} \leq \alpha,}
#' where \eqn{q} is the number of upper-triangular entries tested and
#' \eqn{e^{-t/2}} approximates the tail probability of a \eqn{\chi^2(2)}
#' distribution. Entries with \eqn{T_{ab} \geq \hat{\tau}} are rejected.
#'
#' @param Chi_sq A \eqn{nvar \times nvar} matrix of chi-squared test statistics,
#'   typically the \code{Chi_sq} component of a \code{\link{spec.test}} result
#'   evaluated at \code{Truth = 0}.
#' @param alpha A number in \eqn{(0,1)}. Target FDR level.
#' @param diag Logical. If \code{FALSE} (default), diagonal entries are
#'   excluded from the multiple testing procedure. If \code{TRUE}, they are
#'   included.
#' @param Truth An optional \eqn{nvar \times nvar} logical or binary matrix
#'   indicating the true support (\code{TRUE}/\code{1} for non-zero entries).
#'   When supplied, the empirical FDR and power are computed. Default
#'   \code{NULL} (only the threshold and decision matrix are returned).
#' @param grid_size Integer. Number of points in the threshold grid.
#'   Default \code{100}.
#'
#' @return A list of class \code{"specfdr"} with the following components:
#' \describe{
#'   \item{\code{Decision}}{A \eqn{nvar \times nvar} binary matrix with \code{1}
#'     for rejected entries.}
#'   \item{\code{tau}}{The selected threshold \eqn{\hat{\tau}}, or \code{NA}
#'     if no threshold satisfies the FDR constraint.}
#'   \item{\code{FDR}}{Empirical FDR (only present when \code{Truth} is
#'     supplied).}
#'   \item{\code{Power}}{Empirical power (only present when \code{Truth} is
#'     supplied).}
#' }
#'
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @references Deb, N., Kim Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords models multiple testing FDR spectral precision matrix
#' @examples
#' library(mvtnorm)
#' p <- 10; n <- 200
#' set.seed(42)
#' X    <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))
#' m    <- floor(sqrt(n)); j <- 1; alpha <- 0.05
#' dft  <- dft.all(X)
#' fhat <- fhat_at(dft, j, m)
#' fit  <- cglasso(S = fhat, m = m)
#' res  <- decglasso(object = fit, fhat = fhat)
#' vc   <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m,
#'                 type = "plug-in")
#' st   <- spec.test(Est = res$Theta_tilde, varcov = vc, m = m, alpha = alpha)
#'
#' # FDR-controlled test (no truth known):
#' fdr_res <- spec.fdr(Chi_sq = st$Chi_sq, alpha = alpha, diag = FALSE)
#' fdr_res$tau
#' fdr_res$Decision
#' @export spec.fdr
spec.fdr <- function(Chi_sq,
                     alpha,
                     diag      = FALSE,
                     Truth     = NULL,
                     grid_size = 100) {
  
  # Input validation
  if (!is.matrix(Chi_sq)) {
    stop("'Chi_sq' must be a matrix.")
  }
  p <- nrow(Chi_sq)
  if (ncol(Chi_sq) != p) {
    stop("'Chi_sq' must be square.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a number in (0, 1).")
  }
  if (!is.null(Truth)) {
    if (!is.matrix(Truth) || nrow(Truth) != p || ncol(Truth) != p) {
      stop("'Truth' must be a ", p, " x ", p, " matrix.")
    }
  }
  
  # Threshold and decision matrix
  tau_hat <- compute_tau_inf(Chi_sq, alpha, diag, grid_size = grid_size)
  
  Decision <- matrix(0L, p, p)
  if (!is.na(tau_hat)) {
    Decision[Chi_sq >= tau_hat] <- 1L
  }
  
  output <- list(Decision = Decision, tau = tau_hat)
  
  # Empirical FDR and power (only when Truth is supplied)
  if (!is.null(Truth)) {
    fdr_pw <- FDR_power(Est = Decision, Truth = Truth, diag = diag)
    output$FDR   <- fdr_pw$FDR
    output$Power <- fdr_pw$Power
  }
  
  class(output) <- "specfdr"
  return(output)
}

####################################################################
# Compute entry-wise Z, chi-squared, CI half-widths, and ellipse areas.
stat_cmp <- function(Est, Truth, var_dcmp, bw, alpha, eps = 1e-8) {
  
  p <- nrow(Est)
  Re_diff <- Re(Est) - Re(Truth)
  Im_diff <- Im(Est) - Im(Truth)
  
  Z_re <- Z_im <- Chi_sq <- wing_re <- wing_im <- area <- matrix(NA, p, p)
  
  z_alpha   <- qnorm(1 - alpha / 2)
  chi_alpha <- qchisq(1 - alpha, 2)
  
  for (a in seq_len(p)) {
    for (b in seq_len(p)) {
      
      Vre <- var_dcmp$Vre[a, b]
      Vim <- var_dcmp$Vim[a, b]
      Cov <- var_dcmp$Cov[a, b]
      
      # Completely degenerate: both variances negligible
      if (abs(Vre) < eps && abs(Vim) < eps) {
        Z_re[a, b] <- Z_im[a, b] <- 0
        wing_re[a, b] <- wing_im[a, b] <- 0
        Chi_sq[a, b] <- 0
        area[a, b]   <- 0
        next
      }
      
      # Diagonal entry: imaginary part is identically zero (real-valued self-coherence)
      if (Vim < eps && a == b) {
        var_re_scaled  <- Vre / bw
        Z_re[a, b]     <- Re_diff[a, b] / sqrt(var_re_scaled)
        Z_im[a, b]     <- 0
        wing_re[a, b]  <- z_alpha * sqrt(var_re_scaled)
        wing_im[a, b]  <- 0
        Chi_sq[a, b]   <- Z_re[a, b]^2
        area[a, b]     <- pi * chi_alpha * (Vre / bw)
        next
      }
      
      # General off-diagonal entry: bivariate normal, 2-df chi-squared
      cov_mat    <- matrix(c(Vre, Cov, Cov, Vim), 2, 2)
      cov_scaled <- cov_mat / bw
      
      # Regularise: floor negative eigenvalues at eps
      eig <- eigen(cov_scaled, symmetric = TRUE)
      eig$values[eig$values < eps] <- eps
      cov_reg <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
      
      vec <- c(Re_diff[a, b], Im_diff[a, b])
      
      # Marginal Z-statistics and CI half-widths
      Z_re[a, b]    <- vec[1] / sqrt(cov_reg[1, 1])
      Z_im[a, b]    <- vec[2] / sqrt(cov_reg[2, 2])
      wing_re[a, b] <- z_alpha * sqrt(cov_reg[1, 1])
      wing_im[a, b] <- z_alpha * sqrt(cov_reg[2, 2])
      
      # Mahalanobis chi-squared (2 df)
      Chi_sq[a, b] <- drop(t(vec) %*% solve(cov_reg) %*% vec)
      
      # Area of joint (1-alpha) confidence ellipse: pi * chi_alpha * sqrt(det(Sigma/bw))
      det_cov <- det(cov_mat)
      if (det_cov < eps) {
        area[a, b] <- 0
      } else {
        area[a, b] <- pi * chi_alpha / bw * sqrt(det_cov)
      }
    }
  }
  
  list(Chi_sq  = Chi_sq,  area    = area,
       Z_re    = Z_re,    wing_re = wing_re,
       Z_im    = Z_im,    wing_im = wing_im)
}

####################################################################
# Find the infimum threshold tau that controls FDR at level alpha.
# Uses the chi-squared tail approximation: P(T >= t) ~ exp(-t/2) for T ~ chi^2(2).
compute_tau_inf <- function(Chi_mat, alpha, diag, grid_size = 100) {
  
  p      <- nrow(Chi_mat)
  idx    <- which(upper.tri(Chi_mat, diag = diag), arr.ind = TRUE)
  T_vals <- Chi_mat[idx]
  q      <- length(T_vals)
  t_max  <- 2 * log(q)
  
  t_grid <- seq(1e-8, t_max, length.out = grid_size)
  G_vals <- exp(-t_grid / 2)
  counts <- sapply(t_grid, function(t) sum(T_vals >= t))
  ratio  <- (G_vals * q) / pmax(1, counts)
  
  valid_idx <- which(ratio <= alpha)
  if (length(valid_idx) == 0) {
    return(NA)
  }
  return(t_grid[min(valid_idx)])
}

####################################################################
# Compute empirical FDR and power from binary decision and truth matrices.
FDR_power <- function(Est, Truth, diag) {
  
  idx      <- upper.tri(Truth, diag = diag)
  true_vec <- Truth[idx]
  est_vec  <- Est[idx]
  
  TP <- sum(est_vec != 0 & true_vec != 0)
  FP <- sum(est_vec != 0 & true_vec == 0)
  FN <- sum(est_vec == 0 & true_vec != 0)
  
  FDR   <- ifelse((TP + FP) == 0, 0, FP / (TP + FP))
  Power <- ifelse(sum(true_vec != 0) == 0, NA, TP / sum(true_vec != 0))
  
  list(FDR = FDR, Power = Power)
}