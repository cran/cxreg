#' Simulated data for the cxreg vignette
#'
#' Simple simulated data, used to demonstrate the features of cxreg
#'
#' These datasets are artificial, and are used to test out some of the
#' features of cxreg.
#' @name cxreg
#' @aliases x y cxreg
#' @format Data objects used to demonstrate features in the cxreg vignette
#' @keywords datasets
#' @useDynLib cxreg
#' @import methods
#' @import Matrix
#' @importFrom utils packageDescription
#' @importFrom graphics abline axis matplot points segments text par plot
#' @importFrom stats approx coef median predict qnorm qchisq runif weighted.mean family rnorm gaussian
#' @importFrom grDevices rainbow
#' @importFrom mvtnorm rmvnorm
#' @importFrom gdata upperTriangle
#' @examples
#' \donttest{ 
#' data(classo_example)
#' x <- classo_example$x
#' y <- classo_example$y
#' classo(x, y)
#' 
#' data(cglasso_example)
#' f_hat <- cglasso_example$f_hat
#' m     <- floor(sqrt(cglasso_example$n))
#' cglasso(S = f_hat, m = m, type = "I")
#' cglasso(S = f_hat, m = m, type = "II")
#' }
NULL

#' Internal classo functions
#'
#' @description
#' These are not intended for use by users. \code{lambda.interp} does linear
#' interpolation of the lambdas to obtain a prediction at a new point s.
#' \code{nonzeroCoef} determines in an efficient manner which variables are
#' nonzero in each fit.
#'
#' @name classo-internal
#' @aliases cvtype cvstats
#' @author Younghoon Kim
#' @keywords internal
NULL

#' Complex-valued Lasso, graphical Lasso, and spectral inference
#'
#' This package fits a complex-valued Lasso for regression using coordinate
#' descent. The algorithm is extremely fast and exploits sparsity in the input
#' \code{x} matrix where it exists. A variety of predictions can be made from
#' the fitted models.
#'
#' This package also provides fitting for a complex-valued graphical Lasso
#' (CGLASSO) using coordinate descent. The function is built upon
#' \code{classo} with covariate updates, just as the regular real-valued
#' coordinate descent algorithm for graphical Lasso.
#'
#' Beyond estimation, the package implements a full inference pipeline for
#' high-dimensional sparse spectral precision matrices:
#' \itemize{
#'   \item \code{\link{select_m}}: data-driven bandwidth selection via
#'     generalised cross-validation (GCV) on the diagonal periodogram.
#'   \item \code{\link{decglasso}}: one-step debiasing (desparsification)
#'     of the CGLASSO estimator.
#'   \item \code{\link{var.cov}}: asymptotic variance and pseudovariance
#'     estimation (plug-in or HAC) with real/imaginary decomposition.
#'   \item \code{\link{spec.test}}: entry-wise Z-statistics, Mahalanobis
#'     chi-squared statistics, and confidence ellipse areas.
#'   \item \code{\link{spec.fdr}}: FDR-controlled multiple testing for
#'     spectral precision matrix support recovery.
#' }
#'
#' \tabular{ll}{
#' Package: \tab cxreg \cr
#' Type: \tab Package \cr
#' Version: \tab 1.1.0 \cr
#' Date: \tab 2026-06-01 \cr
#' License: \tab MIT + file LICENSE \cr
#' }
#'
#' @name cxreg-package
#' @author Younghoon Kim, Navonil Deb, Sumanta Basu \cr Maintainer:
#' Younghoon Kim <ykim124@ua.edu>
#' @references Deb, N., Kuceyeski, A., Basu, S. (2024)
#' \emph{Regularized Estimation of Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2401.11128}.
#'
#' Deb, N., Kim, Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords models regression package
#' @examples
#' \donttest{
#' ## --- classo: complex-valued lasso regression ---
#' set.seed(1234)
#' n <- 100; p <- 20
#' x <- array(rnorm(n*p), c(n,p)) + (1+1i) * array(rnorm(n*p), c(n,p))
#' for (j in 1:p) x[,j] <- x[,j] / sqrt(mean(Mod(x[,j])^2))
#' e <- rnorm(n) + (1+1i) * rnorm(n)
#' b <- c(1, -1, rep(0, p-2)) + (1+1i) * c(-0.5, 2, rep(0, p-2))
#' y <- x %*% b + e
#' fit <- classo(x, y)
#' predict(fit, newx = x[1:5, ], s = c(0.01, 0.005))
#' predict(fit, type = "coef")
#' plot(fit, xvar = "lambda")
#' }
#'
#' \donttest{
#' ## --- cglasso: complex-valued graphical lasso ---
#' library(mvtnorm)
#' p <- 30; n <- 500
#' C <- diag(0.7, p)
#' C[row(C) == col(C) + 1] <- 0.3
#' C[row(C) == col(C) - 1] <- 0.3
#' Sigma <- solve(C)
#' set.seed(1010)
#' X_t <- rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
#' m <- floor(sqrt(n)); j <- 1
#' d_j <- dft.j(X_t, j, m)
#' f_j_hat <- t(d_j) %*% Conj(d_j) / (2*m + 1)
#' fit <- cglasso(S = f_j_hat, m = m, type = "I")
#' }
#'
#' \donttest{
#' ## --- Inference pipeline: debiasing, variance, testing ---
#' library(mvtnorm)
#' p <- 10; n <- 200
#' set.seed(42)
#' X <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))
#'
#' ## Bandwidth selection
#' bw_sel <- select_m(X)
#' m <- bw_sel$m_opt
#'
#' ## Spectral density estimate and cglasso fit
#' j <- floor(n / 4)
#' dft  <- dft.all(X)
#' fhat <- fhat_at(dft, j, m)
#' fit  <- cglasso(S = fhat, m = m)
#'
#' ## Debiasing
#' res  <- decglasso(object = fit, fhat = fhat)
#'
#' ## Variance estimation
#' vc <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m, type = "plug-in")
#'
#' ## Confidence regions and test statistics (H0: Theta = 0)
#' st <- spec.test(Est = res$Theta_tilde, varcov = vc, m = m, alpha = 0.05)
#'
#' ## FDR-controlled support recovery
#' fdr_res <- spec.fdr(Chi_sq = st$Chi_sq, alpha = 0.05, diag = FALSE)
#' fdr_res$tau
#' fdr_res$Decision
#' }
NULL