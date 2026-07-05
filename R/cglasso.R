#' fit a complex-valued graphical lasso
#'
#' Fit a complex-valued graphical lasso for spectral precision matrix (inverse spectral density matrix) via a complex variable-wise coordinate descent algorithm for classo with covariates update.
#'
#' The sequence of models implied by \code{lambda} is fit by coordinate descent.
#'
#' @param S nvar x nvar-dimensional symmetric spectral density (or spectral coherence) matrix. S is considered as being computed by average smoothed periodogram (the bandwidth is computed by using the given nobs).
#' @param D The nvar x nvar-dimensional diagonal matrix to produce the partial spectral coherence matrix from S. Default is \code{NULL}. If D is not specified, the function automatically takes the inverse square root of the diagonals of S.
#' @param type A logical flag to choose the formulation to solve. Default is \code{I}. If type is \code{I}, the algorithm solves CGLASSO-I in the reference, 
#' \deqn{ D^{-1/2} \left( \arg\min_{\Theta} \mathrm{Tr} \left[ \hat{R} \hat{\Theta} \right] - \log \det \Theta + \sum_{i \ne j} \left| \Theta_{ij} \right| \right) D^{-1/2} }
#' for the given D. If type is \code{II}, the algorithm solves CGLASSO-II in the reference. It is for each iterative classo with covariate update, the squared-root of scale matrix \eqn{ D^{-1/2}} is multiplied. Please refer to Algorithm 2 in the reference for the details.
#' @param m Window size. This quantity is need to compute the Fourier frequency, extended BIC, and bandwidth for the average smoothed periodogram.
#' @param lambda A user supplied \code{lambda} sequence.
#' Typical usage is to have the program compute its own \code{lambda} sequence based on
#' \code{nlambda} and \code{lambda.min.ratio}.
#' Supplying a value of \code{lambda} overrides this.
#' WARNING: use with care. Avoid supplying a single value for \code{lambda}
#' @param nlambda The number of \code{lambda} values - default is 50.
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a fraction of
#' \code{lambda.max}, the (data derived) entry value (i.e. the smallest value
#' for which all coefficients are zero). ???
#' @param W.init Logical flag whether the initially estimated spectral density matrix is given. Default is \code{NULL}.
#' @param stopping_rule Logical flag if the algorithm is terminated by stopping rule. If the algorithm is early terminated ,not all estimates for initially designated lambdas are explored. 
#' @param stop_criterion Stopping criterion for early termination. Default is \code{EBIC} (Extended BIC). Alternatively, \code{AIC} (AIC) and \code{RMSE} (root mean squared error between two consecutive estimates) can be used.
#' @param maxit Maximum number of iterations of both outer and inner loops. Default 500.
#' @param thresh Convergence threshold for coordinate descent. Default is 1e-4.
#' @param trace.it If \code{trace.it=1}, then a progress bar is displayed;
#' @param \dots Other arguments that can be passed to \code{cglasso}
#' useful for big models that take a long time to fit.
#' 
#' @return An object with class \code{"cglassofit"} and \code{"cglasso"}.
#' \describe{
#' \item{stop_arr}{Sequence of values of information criterion for a fixed lambda.}
#' \item{stop_criterion}{Stopping criterion used.}
#' \item{min_index}{The index for lambda that minimizes the value of the information criterion.}
#' \item{lambda_grid}{Sequence of lambdas used.}
#' \item{Theta_list}{Estimated inverse spectral matrix for each fixed lambda. It is provided in the list.}
#' \item{type}{Type of the formulation used, either CGLASSO-I or CGLASSO-II.}
#' \item{D}{Used scale diagonal matrix.}
#' }
#' 
#' @importFrom mvtnorm rmvnorm
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @references Deb, N., Kim, Y., Basu, S. (2026)
#' \emph{Inference for High-Dimensional Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2606.07986}.
#' @keywords models complex-valued precision matrix
#' @examples
#' p <- 30
#' n <- 500
#' C <- diag(0.7, p)
#' C[row(C) == col(C) + 1] <- 0.3  
#' C[row(C) == col(C) - 1] <- 0.3  
#' Sigma <- solve(C)
#' set.seed(1010)
#' m <- floor(sqrt(n)); j <- 1
#' X_t <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
#' d_j <- dft.j(X_t,j,m)
#' f_j_hat <- t(d_j) %*% Conj(d_j) / (2*m+1)
#' fit <- cglasso(S=f_j_hat, m=floor(sqrt(n)))
#' @export cglasso
cglasso <- function(S,
                    D = NULL,
                    type = c("I","II"),
                    m,
                    lambda = NULL,
                    nlambda = 50,
                    lambda.min.ratio = 1e-2,
                    W.init = NULL,
                    stopping_rule = TRUE,
                    stop_criterion = c("EBIC","AIC","RMSE"),
                    maxit = 500,
                    thresh = 1e-4,
                    trace.it = 0,...){
  
  this.call <- match.call()
  type <- match.arg(type, choices = c("I","II"))
  stop_criterion <- match.arg(stop_criterion, choices = c("EBIC", "AIC", "RMSE"))
  
  # check for NAs in S
  if (any(is.na(S))) {
    stop("S has missing values; S needs to be complete")
  } else {
    p1 <- dim(S)[1]; p2 <- dim(S)[2]
    if (p1 != p2) {
      stop("S must be square")
    } else {
      if (!isSymmetric(S)) {
        stop("S must be symmetric")
      } else {
        p <- p1
      }
    }
  }
  
  # check that S is complex-valued
  if (!is.complex(S)) {
    stop("S must be a complex-valued matrix; real-valued inputs are not supported")
  }
  
  # check that m is provided
  if (is.null(m)) {
    stop("m is missing; the window size must be given.")
  }
  
  # check that m is a positive integer
  if (!is.numeric(m) || length(m) != 1 || m <= 0 || m != round(m)) {
    stop("m must be a positive integer")
  }
  
  # check if D is properly given under the different settings
  if (!is.null(D)) {
    if (is.vector(D)) {
      if (length(D) != p) {
        stop("D must have the same dimension as S")
      } else {
        D <- diag(D, p)
      }
    } else if (is.matrix(D)) {
      if ((dim(D)[1] != p) | (dim(D)[2] != p)) {
        stop("D must be nvar by nvar matrix")
      } else if (!all(D[row(D) != col(D)] == 0)) {
        stop("D must be diagonal")
      }
    } else {
      stop("D must be either nvar length of vector or nvar by nvar diagonal matrix")
    }
  }
  
  fit_cglasso <- cglasso.path(S,
                              D,
                              type,
                              m,
                              lambda,
                              nlambda,
                              lambda.min.ratio,
                              W.init,
                              stopping_rule,
                              stop_criterion,
                              maxit = maxit,
                              thresh = thresh,
                              trace.it = trace.it)
  
  fit_cglasso$call <- this.call
  fit_cglasso
}