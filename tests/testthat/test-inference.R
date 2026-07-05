# Shared setup: small multivariate Gaussian time series
# p=10, n=200 keeps tests fast while exercising all code paths

make_setup <- function(p = 10, n = 200, seed = 42) {
  set.seed(seed)
  # AR-like sparse precision structure
  C <- diag(0.8, p)
  C[row(C) == col(C) + 1] <- -0.2
  C[row(C) == col(C) - 1] <- -0.2
  Sigma <- solve(C)
  D_chol <- chol(Sigma)
  X <- matrix(rnorm(n * p), n, p) %*% D_chol
  
  m    <- floor(sqrt(n))
  j    <- floor(n / 4)
  dft  <- dft.all(X)
  fhat <- fhat_at(dft, j, m)
  fit  <- cglasso(S = fhat, m = m)
  res  <- decglasso(object = fit, fhat = fhat)
  vc_plug <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m,
                     type = "plug-in")
  vc_hac  <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m,
                     type = "HAC")
  
  list(p = p, n = n, m = m, j = j, X = X,
       fhat = fhat, fit = fit, res = res,
       vc_plug = vc_plug, vc_hac = vc_hac,
       Sigma = Sigma, C = C)
}

###################################################################
test_that("decglasso returns correct structure and class", {
  s <- make_setup()
  
  res <- s$res
  p   <- s$p
  
  expect_s3_class(res, "decglasso")
  expect_named(res, c("Theta_tilde", "Theta_hat", "fhat"))
  
  # Theta_tilde is a p x p complex matrix
  expect_true(is.matrix(res$Theta_tilde))
  expect_true(is.complex(res$Theta_tilde))
  expect_equal(dim(res$Theta_tilde), c(p, p))
  
  # Theta_hat is a p x p complex matrix
  expect_true(is.matrix(res$Theta_hat))
  expect_true(is.complex(res$Theta_hat))
  expect_equal(dim(res$Theta_hat), c(p, p))
  
  # fhat stored matches what was passed in
  expect_equal(dim(res$fhat), c(p, p))
  
  # Debiasing formula: Theta_tilde = 2*Theta_hat - Theta_hat %*% fhat %*% Theta_hat
  Theta_hat <- res$Theta_hat
  fhat      <- res$fhat
  expected  <- 2 * Theta_hat - Theta_hat %*% fhat %*% Theta_hat
  expect_equal(res$Theta_tilde, expected)
})

###################################################################
test_that("decglasso accepts a raw complex matrix as object", {
  s <- make_setup()
  
  Theta_hat <- s$fit$Theta_list[[s$fit$min_index]]
  res2 <- decglasso(object = Theta_hat, fhat = s$fhat)
  
  expect_s3_class(res2, "decglasso")
  expect_equal(dim(res2$Theta_tilde), c(s$p, s$p))
})

###################################################################
test_that("decglasso respects the index argument", {
  s <- make_setup()
  
  res1 <- decglasso(object = s$fit, fhat = s$fhat, index = 1)
  res2 <- decglasso(object = s$fit, fhat = s$fhat, index = 2)
  
  # Different lambda -> different Theta_hat -> different Theta_tilde
  expect_false(isTRUE(all.equal(res1$Theta_tilde, res2$Theta_tilde)))
})

###################################################################
test_that("var.cov plug-in returns correct structure and class", {
  s <- make_setup()
  
  vc <- s$vc_plug
  p  <- s$p
  
  expect_s3_class(vc, "varcov")
  expect_named(vc, c("Vre", "Vim", "Cov", "V", "PV", "type"))
  expect_equal(vc$type, "plug-in")
  
  # All matrices are p x p
  for (nm in c("Vre", "Vim", "Cov")) {
    expect_true(is.matrix(vc[[nm]]))
    expect_equal(dim(vc[[nm]]), c(p, p),
                 info = paste("Dimension check for", nm))
    # Real-valued decomposition components
    expect_true(is.numeric(vc[[nm]]),
                info = paste(nm, "should be real"))
  }
  
  # Raw V and PV are complex
  expect_true(is.complex(vc$V))
  expect_true(is.complex(vc$PV))
  expect_equal(dim(vc$V),  c(p, p))
  expect_equal(dim(vc$PV), c(p, p))
  
  # Decomposition identities: Vre = 0.5*Re(V+PV), Vim = 0.5*Re(V-PV)
  expect_equal(vc$Vre, 0.5 * Re(vc$V + vc$PV), tolerance = 1e-10)
  expect_equal(vc$Vim, 0.5 * Re(vc$V - vc$PV), tolerance = 1e-10)
  expect_equal(vc$Cov, 0.5 * Im(vc$PV),         tolerance = 1e-10)
  
  # Variance matrices should be non-negative on the diagonal
  expect_true(all(diag(vc$Vre) >= 0))
})

###################################################################
test_that("var.cov HAC returns correct structure and class", {
  s <- make_setup()
  
  vc <- s$vc_hac
  p  <- s$p
  
  expect_s3_class(vc, "varcov")
  expect_equal(vc$type, "HAC")
  expect_equal(dim(vc$Vre), c(p, p))
  expect_equal(dim(vc$V),   c(p, p))
  expect_true(is.complex(vc$V))
  expect_true(is.numeric(vc$Vre))
})

###################################################################
test_that("spec.test returns correct structure and class", {
  s     <- make_setup()
  alpha <- 0.05
  p     <- s$p
  
  st <- spec.test(Est    = s$res$Theta_tilde,
                  varcov = s$vc_plug,
                  m      = s$m,
                  alpha  = alpha)
  
  expect_s3_class(st, "spectest")
  expect_named(st, c("Chi_sq", "area", "Z_re", "wing_re", "Z_im", "wing_im"))
  
  for (nm in names(st)) {
    expect_true(is.matrix(st[[nm]]),
                info = paste(nm, "should be a matrix"))
    expect_equal(dim(st[[nm]]), c(p, p),
                 info = paste("Dimension check for", nm))
    expect_true(is.numeric(st[[nm]]),
                info = paste(nm, "should be real-valued"))
  }
  
  # Chi_sq values are non-negative
  expect_true(all(st$Chi_sq >= 0, na.rm = TRUE))
  
  # CI half-widths are non-negative
  expect_true(all(st$wing_re >= 0, na.rm = TRUE))
  expect_true(all(st$wing_im >= 0, na.rm = TRUE))
  
  # Ellipse areas are non-negative
  expect_true(all(st$area >= 0, na.rm = TRUE))
  
  # Diagonal entries: imaginary part of Theta is zero -> wing_im = 0
  expect_true(all(st$wing_im[cbind(seq_len(p), seq_len(p))] == 0))
})

###################################################################
test_that("spec.test with Truth = 0 gives same result as default", {
  s     <- make_setup()
  alpha <- 0.05
  p     <- s$p
  
  st_default <- spec.test(Est    = s$res$Theta_tilde,
                          varcov = s$vc_plug,
                          m      = s$m,
                          alpha  = alpha)
  
  st_explicit <- spec.test(Est    = s$res$Theta_tilde,
                           varcov = s$vc_plug,
                           m      = s$m,
                           alpha  = alpha,
                           Truth  = matrix(0 + 0i, p, p))
  
  expect_equal(st_default$Chi_sq,  st_explicit$Chi_sq)
  expect_equal(st_default$wing_re, st_explicit$wing_re)
})

###################################################################
test_that("spec.fdr returns correct structure without Truth", {
  s     <- make_setup()
  alpha <- 0.05
  
  st  <- spec.test(Est    = s$res$Theta_tilde,
                   varcov = s$vc_plug,
                   m      = s$m,
                   alpha  = alpha)
  
  fdr <- spec.fdr(Chi_sq = st$Chi_sq, alpha = alpha, diag = FALSE)
  
  expect_s3_class(fdr, "specfdr")
  expect_named(fdr, c("Decision", "tau"))
  
  # Decision is a binary integer matrix
  expect_true(is.matrix(fdr$Decision))
  expect_equal(dim(fdr$Decision), c(s$p, s$p))
  expect_true(all(fdr$Decision %in% c(0L, 1L)))
  
  # tau is either NA or a positive number
  if (!is.na(fdr$tau)) {
    expect_true(is.numeric(fdr$tau))
    expect_true(fdr$tau > 0)
  }
})

###################################################################
test_that("spec.fdr with Truth returns FDR and Power", {
  s     <- make_setup()
  alpha <- 0.05
  p     <- s$p
  
  st    <- spec.test(Est    = s$res$Theta_tilde,
                     varcov = s$vc_plug,
                     m      = s$m,
                     alpha  = alpha)
  
  # True support: off-diagonal entries of C that are nonzero
  Truth <- (s$C != 0) * 1L
  
  fdr <- spec.fdr(Chi_sq = st$Chi_sq, alpha = alpha,
                  diag = FALSE, Truth = Truth)
  
  expect_s3_class(fdr, "specfdr")
  expect_named(fdr, c("Decision", "tau", "FDR", "Power"))
  
  # FDR in [0, 1]
  expect_gte(fdr$FDR, 0)
  expect_lte(fdr$FDR, 1)
  
  # Power in [0, 1] or NA (if no true positives)
  if (!is.na(fdr$Power)) {
    expect_gte(fdr$Power, 0)
    expect_lte(fdr$Power, 1)
  }
})

###################################################################
test_that("spec.fdr diag=TRUE includes diagonal entries", {
  s     <- make_setup()
  alpha <- 0.05
  
  st <- spec.test(Est    = s$res$Theta_tilde,
                  varcov = s$vc_plug,
                  m      = s$m,
                  alpha  = alpha)
  
  fdr_off  <- spec.fdr(Chi_sq = st$Chi_sq, alpha = alpha, diag = FALSE)
  fdr_diag <- spec.fdr(Chi_sq = st$Chi_sq, alpha = alpha, diag = TRUE)
  
  # Including diagonal can only increase or maintain rejections
  expect_gte(sum(fdr_diag$Decision), sum(fdr_off$Decision))
})

###################################################################
test_that("spec.test input validation catches bad arguments", {
  s <- make_setup()
  p <- s$p
  
  # Non-complex Est
  expect_error(
    spec.test(Est    = Re(s$res$Theta_tilde),
              varcov = s$vc_plug, m = s$m, alpha = 0.05),
    "'Est' must be a complex-valued matrix"
  )
  
  # Wrong class for varcov
  expect_error(
    spec.test(Est    = s$res$Theta_tilde,
              varcov = list(Vre = matrix(0, p, p)),
              m      = s$m, alpha = 0.05),
    "'varcov' must be a 'varcov' object"
  )
  
  # alpha out of range
  expect_error(
    spec.test(Est    = s$res$Theta_tilde,
              varcov = s$vc_plug, m = s$m, alpha = 1.5),
    "'alpha' must be a number in"
  )
})

###################################################################
test_that("var.cov input validation catches bad arguments", {
  s <- make_setup()
  
  # Non-complex Theta
  expect_error(
    var.cov(Theta = Re(s$res$Theta_tilde), X = s$X, j = s$j, m = s$m),
    "'Theta' must be a complex-valued matrix"
  )
  
  # j out of range
  expect_error(
    var.cov(Theta = s$res$Theta_tilde, X = s$X, j = s$n + 1, m = s$m),
    "'j' must be a single integer"
  )
  
  # Bandwidth too large
  expect_error(
    var.cov(Theta = s$res$Theta_tilde, X = s$X, j = s$j, m = s$n),
    "Bandwidth 2\\*m\\+1 must be smaller than n"
  )
})