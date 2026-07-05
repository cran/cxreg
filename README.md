# cxreg: Complex-valued Lasso and Spectral Inference

`cxreg` provides efficient estimation and inference procedures for
complex-valued regression and spectral precision matrix problems.

- **CLASSO**: complex-valued lasso for penalised regression via pathwise
  coordinate descent.
- **CGLASSO**: complex-valued graphical lasso for estimating sparse spectral
  precision matrices (inverse spectral density matrices).
- **Spectral inference pipeline**: bandwidth selection, one-step debiasing,
  asymptotic variance estimation, entry-wise confidence regions, and
  FDR-controlled hypothesis testing for high-dimensional spectral precision
  matrices.

Details can be found in Deb, Kuceyeski, and Basu (2024) and Deb, Kim, and
Basu (2026).

---

## Installation

### Windows

Before installing via `devtools`, ensure you have the correct version of
**Rtools** installed. Rtools is required to compile C/C++/Fortran code in R
packages on Windows.

- Download Rtools from: <https://cran.r-project.org/bin/windows/Rtools/>
- Match the Rtools version to your R version (e.g. Rtools44 for R 4.4.x)

### macOS

Install GCC (which includes `gfortran`) via Homebrew:

``` bash
brew install gcc
```

Check whether `~/.R` exists:

``` bash
ls -ld ~/.R
```

If not, create it and open `~/.R/Makevars`:

``` bash
mkdir ~/.R
vim ~/.R/Makevars
```

Add the following lines, replacing `15.0.1` with the version returned by
`gfortran --version`:

``` bash
FC    = /opt/homebrew/Cellar/gcc/15.0.1/bin/gfortran
F77   = /opt/homebrew/Cellar/gcc/15.0.1/bin/gfortran
FLIBS = -L/opt/homebrew/Cellar/gcc/15.0.1/lib/gcc/15
```

### Install from GitHub

``` r
# install.packages("devtools")
devtools::install_github("yk748/cxreg")
```

---

## Quick start

### CLASSO: complex-valued lasso regression

``` r
library(cxreg)

data(classo_example)
x <- classo_example$x   # complex matrix, 1000 x 200
y <- classo_example$y   # complex vector, length 1000

# Fit lasso path
fit <- classo(x, y)
print(fit)
plot(fit, xvar = "lambda")

# Cross-validation
cvfit <- cv.classo(x, y, trace.it = 1)
print(cvfit)
plot(cvfit)

# Coefficients and predictions at optimal lambda
coef(cvfit, s = "lambda.min")
predict(cvfit, newx = x[1:5, ], s = "lambda.min")
```

### CGLASSO: complex-valued graphical lasso

``` r
data(cglasso_example)
f_hat <- cglasso_example$f_hat   # p x p complex spectral density matrix
m     <- cglasso_example$m       # half-bandwidth

fit <- cglasso(S = f_hat, m = m, type = "I")
print(fit)
plot(fit, index = fit$min_index, type = "mod")
```

### Spectral inference pipeline

``` r
library(mvtnorm)
set.seed(42)
p <- 10; n <- 200
X <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))

# 1. Data-driven bandwidth selection
bw_sel <- select_m(X)
m <- bw_sel$m_opt

# 2. Smoothed periodogram and cglasso fit
j    <- floor(n / 4)
dft  <- dft.all(X)
fhat <- fhat_at(dft, j, m)
fit  <- cglasso(S = fhat, m = m)

# 3. One-step debiasing
res  <- decglasso(object = fit, fhat = fhat)

# 4. Asymptotic variance estimation
vc   <- var.cov(Theta = res$Theta_tilde, X = X, j = j, m = m,
                type = "plug-in")

# 5. Entry-wise test statistics and confidence regions (H0: Theta = 0)
st   <- spec.test(Est = res$Theta_tilde, varcov = vc, m = m, alpha = 0.05)

# 6. FDR-controlled support recovery
fdr  <- spec.fdr(Chi_sq = st$Chi_sq, alpha = 0.05, diag = FALSE)
fdr$tau
fdr$Decision
```

---

## Function reference

| Category | Function | Description |
|---|---|---|
| Regression | `classo()` | Fit complex lasso path |
| Regression | `cv.classo()` | k-fold cross-validation for classo |
| Regression | `coef.classo()` | Extract coefficients |
| Regression | `predict.classo()` | Predict fitted values or coefficients |
| Regression | `print.classo()` | Print coefficient path summary |
| Regression | `plot.classo()` | Plot coefficient paths (Re and Im panels) |
| Graphical lasso | `cglasso()` | Fit complex graphical lasso path |
| Graphical lasso | `plot.cglasso()` | Heatmap of estimated precision matrix |
| Spectral utilities | `dft.all()` | Full normalised DFT of a time series matrix |
| Spectral utilities | `dft.j()` | DFT windowed around a single frequency |
| Spectral utilities | `fhat_at()` | Smoothed periodogram at a frequency index |
| Bandwidth selection | `select_m()` | GCV-based half-bandwidth selection |
| Inference | `decglasso()` | One-step debiased spectral precision estimator |
| Inference | `var.cov()` | Plug-in or HAC variance/pseudovariance estimation |
| Inference | `spec.test()` | Z-statistics, chi-squared statistics, CI half-widths |
| Inference | `spec.fdr()` | FDR-controlled multiple testing |

---

## Parallel cross-validation

`cv.classo()` supports parallel fold fitting via the `parallel = TRUE`
argument. Register a backend before calling:

``` r
library(doParallel)
registerDoParallel(cores = 4)
cvfit <- cv.classo(x, y, parallel = TRUE)
```

Any `foreach`-compatible backend (`doMC`, `doFuture`, etc.) is supported.

---

## Version history

| Version | Date | Notes |
|---|---|---|
| 1.1.0 | 2026-06-01 | Added spectral inference pipeline: `select_m`, `decglasso`, `var.cov`, `spec.test`, `spec.fdr` |
| 1.0.0 | 2025-07-01 | Initial release: CLASSO, CGLASSO, cross-validation, plotting |

---

## Reporting issues

If you encounter a bug or have a feature request, please open an issue at
<https://github.com/yk748/cxreg/issues> or email Younghoon Kim at
<ykim124@ua.edu>.

Please include:

- A minimal reproducible example
- Session info (`sessionInfo()`)
- R version and OS

---

## References

Navonil Deb, Amy Kuceyeski, and Sumanta Basu. 2024. "Regularized
Estimation of Sparse Spectral Precision Matrices." *arXiv preprint
arXiv:2401.11128.* <https://arxiv.org/abs/2401.11128>.

Navonil Deb, Younghoon Kim, and Sumanta Basu. 2026. "Inference for
High-Dimensional Sparse Spectral Precision Matrices." *arXiv preprint
arXiv:2606.07986.* <https://arxiv.org/abs/2606.07986>.

Jerome Friedman, Trevor Hastie, Holger Hofling, and Robert Tibshirani.
2007. "Pathwise Coordinate Optimization." *The Annals of Applied
Statistics* 1(2): 302–332.
<https://doi.org/10.1214/07-AOAS131>

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. 2008. "Sparse
Inverse Covariance Estimation with the Graphical Lasso." *Biostatistics*
9(3): 432–441.
<https://doi.org/10.1093/biostatistics/kxm045>

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. 2010.
"Regularization Paths for Generalized Linear Models via Coordinate
Descent." *Journal of Statistical Software* 33(1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>