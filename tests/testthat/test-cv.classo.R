test_that("cv works", {
  
  set.seed(1234)
  n <- 200
  p <- 50
  x <- array(rnorm(n * p), c(n, p)) + (1 + 1i) * array(rnorm(n * p), c(n, p))
  for (j in 1:p) x[, j] <- x[, j] / sqrt(mean(Mod(x[, j])^2))
  e <- rnorm(n) + (1 + 1i) * rnorm(n)
  b <- c(1, -1, 2, rep(0, p - 3)) + (1 + 1i) * c(-0.5, 2, 3, rep(0, p - 3))
  y <- x %*% b + e
  
  cvfit.test <- cv.classo(x, y, trace.it = 0)
  
  expect_length(cvfit.test$cvm, 100)
  expect_length(cvfit.test$cvsd, 100)
  expect_length(cvfit.test$lambda, 100)
  expect_length(cvfit.test$cvup, 100)
  expect_length(cvfit.test$cvlo, 100)
  expect_length(cvfit.test$nzero, 100)
  expect_true(is.numeric(cvfit.test$lambda.min))
  expect_true(is.numeric(cvfit.test$lambda.1se))
  expect_equal(dim(cvfit.test$index), c(2, 1))
})