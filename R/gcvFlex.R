#' Select the optimal bandwidth \code{m} for spectral density estimation
#'
#' Searches over a grid of half-bandwidths \code{m} and selects the one that
#' minimises a generalised cross-validation (GCV) criterion derived from the
#' gamma deviance of the periodogram, as proposed by Ombao et al. (2001).
#' For each candidate \code{m}, the smoothed diagonal periodogram
#' \eqn{\hat{f}_{jj}(\omega_j)} is compared to the raw diagonal periodogram
#' \eqn{I_{jj}(\omega_j)} via the average gamma deviance
#' \deqn{D(m) = \frac{1}{|\mathcal{J}|}\sum_{j \in \mathcal{J}} w_j
#'   \left[-\log\!\left(\frac{I_{jj}}{f_{jj}}\right)
#'   + \frac{I_{jj} - f_{jj}}{f_{jj}}\right],}
#' and the GCV score is \eqn{\mathrm{GCV}(m) = D(m) / (1 - 1/(2m+1))^2}.
#'
#' @param X A numeric matrix of size \eqn{nobs \times nvar} (time points by variables).
#' @param m_grid An integer vector of candidate half-bandwidths to search over.
#'   Default is a data-driven grid of multiples of \eqn{\sqrt{nobs}},
#'   filtered to satisfy \eqn{m \geq 1} and \eqn{2m+1 < nobs}.
#' @param freq_idx An integer vector of frequency indices (0-based) over which
#'   the GCV criterion is averaged. Default is \code{0:(nobs-1)}, i.e. all
#'   Fourier frequencies. The DC component (\code{j = 0}) and, for even
#'   \code{n}, the Nyquist frequency (\code{j = nobs/2}) are automatically
#'   down-weighted to 0.5 via \code{q_weights}.
#' @param q_weights A numeric vector of non-negative weights of the same
#'   length as \code{freq_idx}, controlling the relative contribution of each
#'   frequency to the GCV criterion. Default weights are 1 for all frequencies,
#'   except 0.5 for the DC component (\code{j = 0}) and for the Nyquist
#'   frequency (\code{j = nobs/2}) when \code{nobs} is even. Weights are internally
#'   normalised to have mean 1.
#' @param lambda_fun A function of the form \code{function(m, p)} returning the
#'   regularisation parameter \eqn{\lambda} to pass to \code{\link{cglasso}}
#'   for a given half-bandwidth \code{m} and number of variables \code{p}. Default is
#'   \code{function(m, p) sqrt(log(p) / (2*m + 1))}.
#' @param eps A small positive number used as a lower bound when clamping
#'   diagonal periodogram values to avoid \code{log(0)}. Default \code{1e-10}.
#' @param verbose Logical. If \code{TRUE} (default), prints a progress line
#'   for each candidate \code{m}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{m_opt}}{The optimal half-bandwidth minimising the GCV score.}
#'   \item{\code{bw_opt}}{The corresponding full bandwidth \code{2*m_opt + 1}.}
#'   \item{\code{lambda_opt}}{The \eqn{\lambda} value returned by
#'     \code{lambda_fun} at \code{m_opt}.}
#'   \item{\code{gcv_min}}{The minimum GCV score.}
#'   \item{\code{gcv_table}}{A \code{data.frame} with columns \code{m},
#'     \code{bw}, \code{lambda}, \code{deviance}, \code{gcv}, and
#'     \code{gcv_denom}, one row per candidate in \code{m_grid}.}
#'   \item{\code{all}}{A list of full output from \code{gcv_theta_one_m} for
#'     each candidate \code{m}.}
#' }
#'
#'
#' @references Ombao, H. C., Raz, J. A., Strawderman, R. L., Von Sachs, R. (2001).
#' \emph{A simple generalised crossvalidation method of span selection for
#' periodogram smoothing.}
#' Biometrika, \strong{88}(4), 1186--1192.
#' \doi{10.1093/biomet/88.4.1186}.
#'
#' @seealso \code{\link{cglasso}}, \code{\link{fhat_at}}, \code{\link{dft.all}}
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @keywords models bandwidth selection spectral density
#' @examples
#' p <- 10
#' n <- 200
#' set.seed(42)
#' X <- matrix(rnorm(n * p), n, p)
#' res <- select_m(X)
#' res$m_opt
#' res$gcv_table
#' @export select_m
####################################################################
select_m <- function(X,
                     m_grid = NULL,
                     freq_idx = NULL,
                     q_weights = NULL,
                     lambda_fun = NULL,
                     eps = 1e-10,
                     verbose = TRUE) {
  
  dft <- dft.all(X)
  n <- nrow(dft)
  p <- ncol(dft)
  
  if (is.null(m_grid)) {
    m_grid <- unique(floor(c(0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3) * sqrt(n)))
    m_grid <- m_grid[m_grid >= 1]
    m_grid <- m_grid[(2 * m_grid + 1) < n]
  }
  
  if (is.null(freq_idx)) {
    freq_idx <- 0:(n - 1)
  }
  
  if (is.null(lambda_fun)) {
    lambda_fun <- function(m, p) {
      sqrt(log(p) / (2 * m + 1))
    }
  }
  
  out_list <- vector("list", length(m_grid))
  
  for (ii in seq_along(m_grid)) {
    m <- m_grid[ii]
    lambda <- lambda_fun(m, p)
    
    if (verbose) {
      cat(sprintf("GCV diagonal bandwidth search: m = %d, bw = %d\n",
                  m, 2 * m + 1))
    }
    
    out_list[[ii]] <- gcv_theta_one_m(
      dft = dft,
      m = m,
      lambda = lambda,
      freq_idx = freq_idx,
      q_weights = q_weights,
      eps = eps
    )
  }
  
  gcv_table <- data.frame(
    m        = sapply(out_list, `[[`, "m"),
    bw       = sapply(out_list, `[[`, "bw"),
    lambda   = sapply(out_list, `[[`, "lambda"),
    deviance = sapply(out_list, `[[`, "deviance"),
    gcv      = sapply(out_list, `[[`, "gcv"),
    gcv_denom = sapply(out_list, `[[`, "gcv_denom")
  )
  
  best_id <- which.min(gcv_table$gcv)
  
  output <- list(
    m_opt      = gcv_table$m[best_id],
    bw_opt     = gcv_table$bw[best_id],
    lambda_opt = gcv_table$lambda[best_id],
    gcv_min    = gcv_table$gcv[best_id],
    gcv_table  = gcv_table,
    all        = out_list
  )
  class(output) <- "selectm"
  return(output)
}

####################################################################
# Evaluate the GCV criterion at a single bandwidth m.
# Returns a list with m, bw, lambda, deviance, gcv, gcv_denom, freq_idx.
gcv_theta_one_m <- function(dft,
                            m,
                            lambda = NULL,
                            freq_idx = NULL,
                            q_weights = NULL,
                            eps = 1e-10) {
  n <- nrow(dft)
  p <- ncol(dft)
  bw <- 2 * m + 1
  
  if (bw >= n) {
    stop("Bandwidth 2*m+1 must be smaller than n.")
  }
  
  if (is.null(freq_idx)) {
    freq_idx <- 0:(n - 1)
  }
  
  if (is.null(q_weights)) {
    q_weights <- rep(1, length(freq_idx))
    q_weights[freq_idx == 0] <- 0.5
    if (n %% 2 == 0) {
      q_weights[freq_idx == n / 2] <- 0.5
    }
  }
  
  if (length(q_weights) != length(freq_idx)) {
    stop("q_weights must have the same length as freq_idx.")
  }
  
  q_weights <- q_weights / mean(q_weights)
  
  dev_vals <- numeric(length(freq_idx))
  
  for (ii in seq_along(freq_idx)) {
    j <- freq_idx[ii]
    
    I_diag <- diag_periodogram.j(dft, j, eps = eps)
    f_diag <- diag_fhat.j(dft, j, m, eps = eps)
    
    dev_vals[ii] <- q_weights[ii] * mean(-log(I_diag / f_diag) + (I_diag - f_diag) / f_diag)
  }
  
  dev <- mean(dev_vals)
  
  # GCV denominator: (1 - 1/bw)^2 = (2m/(2m+1))^2
  gcv_denom <- (2 * m / (2 * m + 1))^2
  gcv <- dev / gcv_denom
  
  list(
    m         = m,
    bw        = bw,
    lambda    = lambda,
    deviance  = dev,
    gcv       = gcv,
    gcv_denom = gcv_denom,
    freq_idx  = freq_idx
  )
}

####################################################################
# Diagonal of the smoothed periodogram at frequency index j.
# Returns a real vector of length p, clamped below by eps.
diag_fhat.j <- function(dft, j, m, eps = 1e-10) {
  n <- nrow(dft)
  ind <- fixm(j + (-m):m, n)
  Z <- dft[ind, , drop = FALSE]
  f_diag <- colMeans(Mod(Z)^2)
  return(pmax(Re(f_diag), eps))
}

####################################################################
# Diagonal of the raw periodogram at frequency index j.
# Returns a real vector of length p, clamped below by eps.
diag_periodogram.j <- function(dft, j, eps = 1e-10) {
  n <- nrow(dft)
  z <- dft[fixm(j, n), , drop = TRUE]
  return(pmax(Mod(z)^2, eps))
}