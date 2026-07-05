#' Example data for cglasso
#'
#' A synthetic dataset used to illustrate the complex-valued graphical Lasso
#' (\code{\link{cglasso}}). The data consist of a smoothed periodogram matrix
#' at a single Fourier frequency, computed from a simulated multivariate
#' time series with a known sparse spectral precision matrix structure.
#'
#' The dataset is generated from a \eqn{p = 30}-dimensional VMA(3) process
#' with \eqn{n = 500} observations. The half-bandwidth \code{m} is chosen
#' as \code{floor(sqrt(n))} and the smoothed periodogram \code{f_hat} is
#' evaluated at frequency index \code{j = 1}.
#'
#' @docType data
#' @name cglasso_example
#' @usage data(cglasso_example)
#' @format A list with the following components:
#' \describe{
#'   \item{\code{f_hat}}{A \eqn{p \times p} complex-valued Hermitian matrix.
#'     The smoothed periodogram (spectral density) estimate at Fourier
#'     frequency index \code{j}, computed via \code{\link{fhat_at}} with
#'     half-bandwidth \code{m = floor(sqrt(n))}.}
#'   \item{\code{n}}{A positive integer. The number of time series
#'     observations used to generate the data. The half-bandwidth is
#'     \code{m = floor(sqrt(n))}.}
#'   \item{\code{p}}{A positive integer. The number of variables (dimension
#'     of the time series).}
#'   \item{\code{j}}{A positive integer. The Fourier frequency index at which
#'     \code{f_hat} was evaluated.}
#' }
#' @seealso \code{\link{cglasso}}, \code{\link{fhat_at}}, \code{\link{dft.all}}
#' @references Deb, N., Kim, Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords datasets
"cglasso_example"

#' Example data for classo
#'
#' A synthetic dataset used to illustrate the complex-valued Lasso
#' (\code{\link{classo}}). The data consist of a complex-valued design matrix
#' and response vector generated from a sparse linear model.
#'
#' The dataset is generated with \eqn{n = 1000} observations and
#' \eqn{p = 200} variables. The true coefficient vector \eqn{\beta} has
#' two nonzero entries. Each column of \code{x} is normalised to have unit
#' complex magnitude. The response \code{y} is \eqn{X\beta + \varepsilon}
#' with complex Gaussian noise \eqn{\varepsilon}.
#'
#' @docType data
#' @name classo_example
#' @usage data(classo_example)
#' @format A list with the following components:
#' \describe{
#'   \item{\code{x}}{A complex-valued matrix of dimension
#'     \eqn{n \times p} (\eqn{1000 \times 200}). Each column has been
#'     normalised so that \eqn{\|x_j\| / \sqrt{n} = 1}.}
#'   \item{\code{y}}{A complex-valued vector of length \eqn{n}. The response
#'     \eqn{y = X\beta + \varepsilon} where \eqn{\varepsilon} is complex
#'     Gaussian noise.}
#'   \item{\code{beta}}{A complex-valued vector of length \eqn{p}. The true
#'     sparse coefficient vector, with nonzero entries at positions 1 and 2.}
#'   \item{\code{n}}{A positive integer. The number of observations.}
#'   \item{\code{p}}{A positive integer. The number of variables.}
#' }
#' @seealso \code{\link{classo}}, \code{\link{cv.classo}}
#' @references Deb, N., Kuceyeski, A., Basu, S. (2024)
#' \emph{Regularized Estimation of Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2401.11128}.
#' @keywords datasets
"classo_example"