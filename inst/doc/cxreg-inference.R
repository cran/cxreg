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

## -----------------------------------------------------------------------------
library(cxreg)
library(mvtnorm)

## -----------------------------------------------------------------------------
set.seed(1010)
p <- 10     # number of variables
n <- 500    # sample size

# True precision matrix: tridiagonal
C <- diag(0.7, p)
C[row(C) == col(C) + 1] <- 0.3
C[row(C) == col(C) - 1] <- 0.3
Sigma <- solve(C)

# True SPM: constant across all frequencies
Theta_true <- 2*pi*C

# Generate white noise observations
X <- rmvnorm(n = n, mean = rep(0,p), sigma=Sigma)
cat("Data dimensions:", dim(X), "\n")
cat("True SPM (1:4, 1:4):\n")
print(round(Re(Theta_true[1:4, 1:4]), 3))

## -----------------------------------------------------------------------------
bw_sel <- select_m(X, verbose = FALSE)
m <- bw_sel$m_opt
cat("Selected bandwidth: m =", m, "(full bandwidth:", 2*m+1,")\n")
bw_sel$gcv_table

## -----------------------------------------------------------------------------
j <- floor(n / 4)
dft <- dft.all(X) # full DFT matrix: n x p
fhat <- fhat_at(dft, j = j, m = m)  # smoothed periodogram: p x p

cat("Frequency index j =", j, "(omega =", round(2*pi*j/n, 3), ")\n")

# For white noise: fhat should be close to Sigma / (2*pi)
cat("Max deviation of Re(fhat) from Sigma/(2pi):", round(max(abs(Re(fhat) - Sigma / (2*pi))), 4), "\n")

## -----------------------------------------------------------------------------
# Fit CGLASSO
fit <- cglasso(S = fhat, m = m, type = "II")
cat("EBIC-selected lambda index:", fit$min_index, "\n")
cat("Lambda selected:", round(fit$lambda_grid[fit$min_index], 4), "\n")

## -----------------------------------------------------------------------------
# One-step debiasing to produce deCGLASSO
deb <- decglasso(object = fit, fhat = fhat)
cat("Class of deCGLASSO output:", class(deb), "\n")
Theta_tilde <- deb$Theta_tilde

## -----------------------------------------------------------------------------
vc_plug <- var.cov(Theta = Theta_tilde, X = X, j = j, m = m, type = "plug-in")
vc_hac <- var.cov(Theta = Theta_tilde, X = X, j = j, m = m, type = "HAC")
cat("Fields:", names(vc_plug), "\n")

## -----------------------------------------------------------------------------
cat("Entry (1,2) variance (Re part):\n")
cat("Plug-in Vre:", round(vc_plug$Vre[1,2], 6), "\n")
cat("HAC Vre:", round(vc_hac$Vre[1,2],  6), "\n")

## -----------------------------------------------------------------------------
alpha <- 0.2
st_plug <- spec.test(Est = Theta_tilde, varcov = vc_plug, m = m, alpha = alpha)
st_hac <- spec.test(Est = Theta_tilde, varcov = vc_hac, m = m, alpha = alpha)

cat("Chi-squared statistic at (1,2) [true edge]:\n")
cat("Plug-in:", round(st_plug$Chi_sq[1,2], 3), "\n")
cat("HAC:", round(st_hac$Chi_sq[1,2],  3), "\n")
cat("Chi-sq(2, 0.95) critical value:", round(qchisq(0.95, 2), 3), "\n")

cat("\n Chi-squared statistic at (1,3) [true null edge]:\n")
cat("Plug-in:", round(st_plug$Chi_sq[1,3], 3), "\n")
cat("HAC:", round(st_hac$Chi_sq[1,3],  3), "\n")

## -----------------------------------------------------------------------------
cat("80% CI half-width at (1,2):\n")
cat("Re part (plug-in):", round(st_plug$wing_re[1,2], 4), "\n")
cat("Im part (plug-in):", round(st_plug$wing_im[1,2], 4), "\n")

## -----------------------------------------------------------------------------
# Check if true value falls inside the 95% CI
true_re_12 <- Re(Theta_true[1,2])
est_re_12 <- Re(Theta_tilde[1,2])
hw_re_12 <- st_plug$wing_re[1,2]
cat("True Re(Theta_12):", round(true_re_12, 4), "\n")
cat("80% CI: [",round(est_re_12 - hw_re_12, 4), ",", round(est_re_12 + hw_re_12, 4),"]\n")
cat("True value inside CI:", abs(est_re_12 - true_re_12) <= hw_re_12, "\n")

## -----------------------------------------------------------------------------
fdr_plug <- spec.fdr(Chi_sq = st_plug$Chi_sq, alpha = alpha, diag = FALSE)
fdr_hac <- spec.fdr(Chi_sq = st_hac$Chi_sq, alpha = alpha, diag = FALSE)

cat("Number of significant edges (plug-in):", sum(fdr_plug$Decision[upper.tri(fdr_plug$Decision)]), "\n")
cat("Number of significant edges (HAC):", sum(fdr_hac$Decision[upper.tri(fdr_hac$Decision)]), "\n")
cat("FDR threshold tau (plug-in):", round(fdr_plug$tau, 3), "\n")

## -----------------------------------------------------------------------------
# True support: nonzero off-diagonal entries of C (same as Theta_true up to scale)
true_support <- (C!=0)*1L
diag(true_support) <- 0

# Confusion matrix (plug-in)
dec <- fdr_plug$Decision
TP <- sum(dec == 1 & true_support == 1 & upper.tri(dec))
FP <- sum(dec == 1 & true_support == 0 & upper.tri(dec))
FN <- sum(dec == 0 & true_support == 1 & upper.tri(dec))

cat("True edges in upper triangle:", sum(true_support[upper.tri(true_support)]), "\n")
cat("True Positives:", TP, "\n")
cat("False Positives:", FP, "\n")
cat("False Negatives:", FN, "\n")
cat("Empirical FDR:", round(FP / max(TP + FP, 1), 3), "\n")
cat("Power:", round(TP / max(TP + FN, 1), 3), "\n")

