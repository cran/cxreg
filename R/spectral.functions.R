####################################################################
#' Discrete Fourier Transform of matrix X at a fixed frequency omega_j
#' 
#' Computes the (normalized) discrete Fourier transform (DFT) of a matrix \code{X} row-wise using \code{mvfft}, and extracts a window of frequencies centered at a target index.
#'
#' @param X A numeric matrix of size \eqn{nobs \times nvar}, where DFT is applied across the rows (time points).
#' @param j An integer index in \eqn{1,\ldots,nobs} around which the frequency window is centered.
#' @param m A non-negative integer specifying the window half-width. The function returns \code{2m + 1} DFT components centered around \code{j}.
#'
#' @return A complex-valued matrix of dimension \eqn{(2m + 1) \times nvar} representing selected DFT components of the original matrix.
#'
#' @export
dft.j <- function(X, j, m){
  n <- nrow(X)
  dft <- mvfft(X)/sqrt(2*pi*n)
  vj <- j + (-m):m
  ind <- fixm(vj, n)
  Z <- dft[ind, ]
  return(Z)
}

####################################################################
#' Discrete Fourier Transform of matrix X
#' 
#' Computes the (normalized) discrete Fourier transform (DFT) of a matrix \code{X} row-wise using \code{mvfft} over all Fourier frequencies.
#' @param X A numeric matrix of size \eqn{nobs \times nvar}, where DFT is applied across the rows (time points).
#' 
#' @return A complex-valued matrix of dimension \eqn{nobs \times nvar} representing all DFT components of the original matrix.
#'
#' @export
dft.all <- function(X) {
  return(mvfft(X) / sqrt(2*pi*nrow(X)))
}

####################################################################
#' @noRd
fixm <- function(v, n){
  v[v<1] <- v[v<1]+n
  v[v>n] <- v[v>n]-n
  return(v)
}

####################################################################
#' Smoothed periodogram at a fixed frequency
#'
#' Extracts a window of \code{2m + 1} DFT components centered at index \code{j}
#' from a pre-computed full DFT matrix and returns their averaged outer product,
#' i.e. the smoothed periodogram (spectral density) estimate at frequency
#' \eqn{\omega_j = 2\pi j / n}:
#' \deqn{\hat{f}(\omega_j) = \frac{1}{2m+1} \sum_{k=j-m}^{j+m} d_k d_k^*,}
#' where \eqn{d_k} is the \eqn{k}-th row of the DFT matrix (as a column vector)
#' and \eqn{*} denotes the conjugate transpose.
#'
#' @param dft A complex-valued matrix of dimension \eqn{nobs \times nvar},
#'   typically the output of \code{\link{dft.all}}.
#' @param j An integer index in \eqn{1,\ldots,nobs} at which the smoothed
#'   periodogram is evaluated.
#' @param m A non-negative integer specifying the window half-width. Frequency
#'   indices outside \eqn{[1, nobs]} are wrapped circularly.
#'
#' @return A \eqn{nvar \times nvar} complex-valued Hermitian matrix giving the
#'   smoothed periodogram estimate at frequency \eqn{\omega_j}.
#'
#' @seealso \code{\link{dft.all}}, \code{\link{dft.j}}
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @references Deb, N., Kim Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords models spectral density
#' @examples
#' p <- 30
#' n <- 500
#' set.seed(1010)
#' X_t <- matrix(rnorm(n * p), n, p)
#' dft  <- dft.all(X_t)
#' m    <- floor(sqrt(n)); j <- 1
#' fhat <- fhat_at(dft, j, m)
#' @export
fhat_at <- function(dft, j, m){
  n <- nrow(dft)
  vj <- j + (-m):m
  ind <- fixm(vj, n)
  Z <- dft[ind, ]
  return(t(Z) %*% Conj(Z) / (2*m + 1))
}