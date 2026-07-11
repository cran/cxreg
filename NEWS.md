# cxreg 1.1.4 (2026-07-09)

## New vignette

* Added `vignettes/cxreg-inference.Rmd`: **"Inference for Sparse Spectral
  Precision Matrices with `cxreg`"**. Covers the full inference pipeline
  from bandwidth selection through FDR-controlled hypothesis testing,
  illustrated on a multivariate Gaussian white noise example where the
  true spectral precision matrix is known in closed form.

---

# cxreg 1.1.2 (2026-06-27)

## Bug fixes

* `vignettes/cxreg.Rmd`: corrected several API calls that had become
  stale relative to the v1.1.0 changes:
  - `cglasso_example$n` used to derive `m`; the stored object has no
    `$m` field, so `m` is now computed as `floor(sqrt(cglasso_example$n))`.
  - `cglasso(S = f_hat, nobs = n, ...)` → `cglasso(S = f_hat, m = m, ...)`
    throughout the CGLASSO example section.
  - `plot(fit$Theta_list, index = ..., type = ...)` → `plot(fit, ...)`
    (plot dispatch now takes the full `cglassofit` object, not `$Theta_list`).
  - `dft.X(X_t, j, m)` → `dft.j(X_t, j, m)` (old function name).
  - `predict(fit, ..., type = "coefficient")` → `type = "coefficients"`.
  - `\textrm{...}` in LaTeX math environments → `\mathrm{...}` to avoid
    `amsmath` dependency in the PDF build.

## Internal

* `src/classo_init.cpp`: removed all `extern "C"` forward declarations
  of Fortran subroutines. The declarations used `Rcomplex*` argument
  types while the Fortran compiler sees `double complex*`, causing
  `-Wlto-type-mismatch` warnings under LTO on CRAN's Linux checks.
  Symbols are now resolved at runtime via `R_useDynamicSymbols(info, TRUE)`,
  which is the standard approach for packages with Fortran routines.

---

# cxreg 1.1.1 (2026-06-26)

## Bug fixes

* `vignettes/cxreg.pdf`: compacted using
  `tools::compactPDF(gs_quality = "ebook")`, reducing file size from
  457 Kb to 168 Kb to address a CRAN PDF size warning.

---

# cxreg 1.1.0 (2026-06-01)

## New features

* **Spectral inference pipeline** — complete end-to-end workflow for
  high-dimensional sparse spectral precision matrices (Deb, Kim, and
  Basu 2026):
  - `select_m()`: data-driven bandwidth selection via generalised
    cross-validation (GCV) on the diagonal periodogram (Ombao et al. 2001).
  - `decglasso()`: one-step debiased (desparsified) spectral precision
    estimator. Returns an object of class `"decglasso"`.
  - `var.cov()`: asymptotic variance and pseudovariance estimation with
    real/imaginary decomposition. Supports plug-in and HAC estimators.
    Returns an object of class `"varcov"`.
  - `spec.test()`: entry-wise Z-statistics, Mahalanobis chi-squared
    statistics, CI half-widths, and joint confidence ellipse areas.
    Returns an object of class `"spectest"`.
  - `spec.fdr()`: FDR-controlled multiple testing for support recovery
    of the spectral precision matrix. Returns an object of class
    `"specfdr"`.
* `fhat_at()`: smoothed periodogram matrix at a single Fourier frequency,
  added to `spectral_functions.R`.
* **Parallel cross-validation** — `cv.classo(parallel = TRUE)` is now
  fully operational via `foreach` / `%dopar%`. Requires a registered
  backend (e.g. `doParallel::registerDoParallel()`). A graceful fallback
  warning is issued when no backend is registered.

## Bug fixes

* `classo.control()` was ignoring all arguments and always returning
  hardcoded factory defaults, making `cv.classo(trace.it = 1)` have no
  effect. Replaced with an environment-backed settings store so that
  changes (e.g. `classo.control(itrace = 1)`) persist for the R session.
* `print.classo()`: `%Dev` column was reading from the non-existent
  `$dev` field; corrected to `$dev.ratio`.
* `cv.classo()`: stop message incorrectly said `cv.glasso`; corrected
  to `cv.classo`. CV-only arguments (`alignment`, `parallel`) were not
  stripped before forwarding the call to `classo()`.
* `cv.classo.raw()`: `standardized = TRUE` (wrong argument name) changed
  to `standardize = TRUE` in all `classo()` fold calls.
* `cv_classofit()`: removed a `glmnet`-inherited `family$initialize`
  block that was incompatible with complex regression. MSE was computed
  via `abs()` (MAE) rather than `Mod()^2`; corrected.
* `buildPredmat()`: S3 generic was defined after its method; reordered.
  Prediction matrix initialised as real `NA` instead of `NA_complex_`,
  causing silent type-coercion. The alignment `switch()` was commented
  out and bypassed; reactivated.
* `dev_comp()` and `get_start()` in `classoFlex.R`: residuals were split
  into real and imaginary parts incorrectly. Fixed to
  `Mod(y - x %*% beta)^2`. `weighted.mean()` on complex `y` (unsupported
  in base R) replaced with `sum(w * y) / sum(w)`. `t(rv) %*% x` for
  `lambda_max` lacked conjugation; corrected to `Conj(t(rv)) %*% x`.
* `plot.cglasso()`: `x[[index]]` accessed a named field of the fit object
  by position rather than `x$Theta_list[[index]]`. `is.integer(index)`
  rejected all user-supplied plain integers (e.g. `index = 1`); replaced
  with `is.numeric` + round check. Matrix orientation was wrong (`S[,
  nrow(S):1]` reverses columns, not rows); corrected to `t(S)[, p:1]`.
  The imaginary panel in `"both"` mode used `zlim = z_re` instead of
  `zlim = z_im`. `par()` was restored twice (explicitly and via
  `on.exit`).
* `predict.classo()`: `cbind2()` (a sparse-matrix S4 generic) replaced
  with plain matrix multiplication. `"link"` type documented but missing
  from `match.arg` choices. `nonzeroCoef()` always dropped row 1
  assuming an intercept; now conditional on `object$a0`.
* `family.classo()` / `family.classofit()`: were mapping `glmnet`-style
  class names (`elnet`, `lognet`, etc.) that do not exist in `cxreg`,
  returning `NA`. Both now return `"gaussian"` directly.
* `coef.cv.classo()`: missing `names(lambda) <- s` when `s` is a
  character string, making the output unnamed; added to match
  `predict.cv.classo()`.
* `lambda.interp()`: top-level `=` assignment; two index assignments
  using `=`; comment formula `sfrac*left+(1-sfrac*right)` missing
  parentheses. All fixed.
* `plotCoef()`: `switch(length(which) + 1, "0" = ...)` — the `"0"` case
  was unreachable because `length(NULL) + 1 = 1`. Fixed to
  `switch(length(active), "0" = ...)`. The `which` variable was reused
  for original indices and then overwritten inside each panel, causing
  labels to show local rather than original variable indices.
* `cvtype()`: `subclass.ch = c(1, 2, 5)` produced an `NA` entry (only
  2 elements in `type.measures`); changed to `c(1, 2)`. Both entries
  labelled `"Mean-Absolute Error"`; corrected to `"Mean-Squared Error"`.
* `error.bars()`: returned `range(upper, lower)` visibly instead of
  `invisible(NULL)`.

## Internal changes

* Distinct S3 classes assigned to all return objects: `"decglasso"`,
  `"varcov"`, `"selectm"`, `"spectest"`, `"specfdr"` (previously some
  returned plain lists or had the wrong class).
* `classo.path()`: duplicate intercept/`cbind` block removed (would have
  appended the intercept column twice when `intercept = TRUE`). Orphaned
  `classofit` list object removed. `classo.control()` called only once.
* `cxreg-package.R`: removed unused `@import foreach`, `@importFrom
  fields image.plot`, and `@importFrom Rcpp sourceCpp`. Added
  `@importFrom stats qnorm qchisq`.

---

# cxreg 1.0.5

* A minor fix was made during the revision for the Journal of Open
  Statistical Software.

# cxreg 1.0.0

* First stable release: CLASSO and CGLASSO estimation, cross-validation,
  coefficient and prediction methods, regularization path plots, and
  precision matrix heatmaps.

# cxreg 0.8.0

* Pathwise estimation function for complex-valued graphical lasso updated.

# cxreg 0.1.1

* `plot.classo()`: changed to display real and imaginary parts in separate
  panels.

# cxreg 0.1.0

* Initial development version.