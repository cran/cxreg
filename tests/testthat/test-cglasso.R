test_that("cglasso works", {
  
  rmvnorm_manual <- function(n, p, Sigma) {
    D <- chol(Sigma)
    return(matrix(rnorm(n * p), ncol = p) %*% D)
  }
  
  p <- 40
  n <- 500
  C <- diag(0.8, p)
  C[row(C) == col(C) + 1] <- -0.2
  C[row(C) == col(C) - 1] <- -0.2
  Sigma <- solve(C)
  
  set.seed(1234)
  X_t <- rmvnorm_manual(n, p, Sigma)
  m <- floor(sqrt(n))
  j <- 1
  d_j <- dft.j(X_t, j, m)
  f_j_hat <- t(d_j) %*% Conj(d_j) / (2 * m + 1)
  
  fit <- cglasso(S = f_j_hat, m = m, nlambda = 20, stop_criterion = "AIC",
                 stopping_rule = FALSE)
  
  expect_length(fit$stop_arr, 20)
  expect_equal(fit$stop_criterion, "AIC")
  expect_true(is.numeric(fit$min_index))
  expect_length(fit$lambda_grid, 20)
  expect_true(is.list(fit$Theta_list))
  expect_equal(fit$type, "I")
  expect_true(is.matrix(fit$D))
})