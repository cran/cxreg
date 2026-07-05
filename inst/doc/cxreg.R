## ----include=FALSE------------------------------------------------------------
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

## ----eval=FALSE, message=FALSE------------------------------------------------
# library(devtools)
# devtools::install_github("yk748/cxreg")

## -----------------------------------------------------------------------------
library(cxreg)

## -----------------------------------------------------------------------------
data(classo_example)
x <- classo_example$x
y <- classo_example$y

## -----------------------------------------------------------------------------
fit <- classo(x,y)

## -----------------------------------------------------------------------------
plot(fit, xvar="lambda", label=TRUE)

## -----------------------------------------------------------------------------
plot(fit, xvar="norm", label=TRUE)

## -----------------------------------------------------------------------------
plot(fit, xvar="dev", label=TRUE)

## -----------------------------------------------------------------------------
any(fit$lambda == 0.1)
coef(fit, s=0.1, exact=FALSE)

## -----------------------------------------------------------------------------
coef(fit, s=0.1, exact=TRUE, x=x, y=y)

## -----------------------------------------------------------------------------
set.seed(29)
nx <- array(rnorm(5*20), c(5,20)) + (1+1i) * array(rnorm(5*20), c(5,20))
for (j in 1:20) {
  nx[,j] <- nx[,j] / sqrt(mean(Mod(nx[,j])^2))
}
predict(fit, newx = nx, s = c(0.1, 0.05), type="response")

## -----------------------------------------------------------------------------
predict(fit, newx = nx, s = c(0.1, 0.05), type="coefficients")

## -----------------------------------------------------------------------------
predict(fit, newx = nx, s = c(0.1, 0.05), type="nonzero")

## -----------------------------------------------------------------------------
cvfit <- cv.classo(x,y,trace.it = 1)

## -----------------------------------------------------------------------------
print(cvfit)

## -----------------------------------------------------------------------------
plot(cvfit)

## -----------------------------------------------------------------------------
cvfit$lambda.min

## -----------------------------------------------------------------------------
coef(cvfit, s = "lambda.min")

## -----------------------------------------------------------------------------
predict(cvfit, newx = x[1:5,], s = "lambda.min")

## -----------------------------------------------------------------------------
data(cglasso_example)
f_hat <- cglasso_example$f_hat
n     <- cglasso_example$n
m     <- floor(sqrt(n))   # half-bandwidth used to compute f_hat

## -----------------------------------------------------------------------------
fit_cglasso_I <- cglasso(S=f_hat, m=m, type="I")

## -----------------------------------------------------------------------------
fit_cglasso_II <- cglasso(S=f_hat, m=m, type="II", nlambda=30, stop_criterion = "AIC")

## -----------------------------------------------------------------------------
fit_cglasso_another <- cglasso(S=f_hat, m=m, type="II", stopping_rule = FALSE)
fit_cglasso_another$lambda_grid

## -----------------------------------------------------------------------------
plot(fit_cglasso_I, index=fit_cglasso_I$min_index, type="mod",       label=TRUE)
plot(fit_cglasso_I, index=fit_cglasso_I$min_index, type="both",      label=TRUE)
plot(fit_cglasso_II, index=fit_cglasso_II$min_index, type="real",    label=FALSE)
plot(fit_cglasso_II, index=fit_cglasso_II$min_index, type="imaginary", label=FALSE)

## -----------------------------------------------------------------------------
library(mvtnorm)
set.seed(1010)
p <- 10    # number of variables
n <- 500   # sample size

# sparse precision matrix (tridiagonal)
C <- diag(0.7, p)
C[row(C) == col(C) + 1] <- 0.3
C[row(C) == col(C) - 1] <- 0.3
Sigma <- solve(C)

# generate multivariate Gaussian time series
X_t <- rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)

## -----------------------------------------------------------------------------
# compute DFT at a specific Fourier frequency j
m <- floor(sqrt(n))   # bandwidth
j <- 1                # frequency index
d_j <- dft.j(X_t, j, m)

# estimate spectral density via smoothed periodogram
f_j_hat <- t(d_j) %*% Conj(d_j) / (2*m + 1)

## -----------------------------------------------------------------------------
fit_sim <- cglasso(S = f_j_hat, m = m, type = "I")

## -----------------------------------------------------------------------------
plot(fit_sim,
     index = fit_sim$min_index,
     type  = "mod",
     label = TRUE)

## -----------------------------------------------------------------------------
# selected lambda
fit_sim$lambda_grid[fit_sim$min_index]

# number of nonzero off-diagonal entries
Theta_est <- fit_sim$Theta_list[[fit_sim$min_index]]
sum(Mod(Theta_est[row(Theta_est) != col(Theta_est)]) > 1e-6)

