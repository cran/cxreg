test_that("classo works", {
  
  set.seed(1234)
  n <- 500
  p <- 100
  x <- array(rnorm(n * p), c(n, p)) + (1 + 1i) * array(rnorm(n * p), c(n, p))
  for (j in 1:p) x[, j] <- x[, j] / sqrt(mean(Mod(x[, j])^2))
  e <- rnorm(n) + (1 + 1i) * rnorm(n)
  b <- c(1, -1, rep(0, p - 2)) + (1 + 1i) * c(-0.5, 2, rep(0, p - 2))
  y <- x %*% b + e
  fit.test <- classo(x, y, nlambda = 50)
  
  newn <- 10
  newx <- array(rnorm(newn * p), c(newn, p)) + (1 + 1i) * array(rnorm(newn * p), c(newn, p))
  for (j in 1:p) newx[, j] <- newx[, j] / sqrt(mean(Mod(newx[, j])^2))
  resp  <- predict(fit.test, newx, s = c(0.5), type = "response")
  coef  <- predict(fit.test, newx, s = c(0.5), type = "coefficients")
  nzero <- predict(fit.test, newx, s = c(0.5), type = "nonzero")
  
  expect_equal(fit.test$nobs, n)
  expect_length(fit.test$lambda, 50)
  expect_true(is.matrix(fit.test$beta))
  expect_equal(dim(fit.test$beta), c(100, 50))
  expect_length(fit.test$df, 50)
  expect_equal(fit.test$dim, c(100, 50))
  expect_length(fit.test$dev.ratio, 50)
  expect_true(is.numeric(fit.test$nulldev))
  expect_equal(fit.test$nobs, n)
  
  expect_equal(dim(resp), c(newn, 1))
  expect_equal(dim(coef), c(p, 1))
  expect_equal(class(nzero), "list")
})