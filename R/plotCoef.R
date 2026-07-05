#' Plot coefficient paths from a classo fit
#'
#' Internal function called by \code{\link{plot.classo}} to produce the
#' coefficient profile plot. Draws two panels: one for the real parts and one
#' for the imaginary parts of the complex-valued coefficient paths.
#'
#' @param beta A complex-valued coefficient matrix of dimension
#'   \code{nvars x nlambda}, as stored in a \code{"classofit"} object.
#' @param norm Optional numeric vector of L1 norms to use as the x-axis when
#'   \code{xvar = "norm"}. If missing, computed as \code{colSums(Mod(beta))}.
#' @param lambda Numeric vector of lambda values from the fitted model.
#' @param df Integer vector of nonzero coefficient counts at each lambda.
#' @param dev Numeric vector of deviance ratios at each lambda.
#' @param label Logical. If \code{TRUE}, annotate the curves with the original
#'   variable indices at the end of the path. Default \code{FALSE}.
#' @param xvar Character string selecting the x-axis scale: \code{"norm"} for
#'   the L1 norm, \code{"lambda"} for \code{log(lambda)}, or \code{"dev"} for
#'   the fraction of deviance explained.
#' @param xlab X-axis label. Defaults to a label matching \code{xvar}.
#' @param \dots Additional graphical arguments passed to \code{matplot}.
#' @return \code{NULL}, invisibly. Called for its side effect of producing a
#'   two-panel plot.
#' @keywords internal
plotCoef <- function(beta, norm, lambda, df, dev, label = FALSE,
                     xvar = c("norm", "lambda", "dev"),
                     xlab = iname, ...) {
  
  # Identify variables with at least one nonzero coefficient across the path
  active <- nonzeroCoef(beta)
  
  # Guard: nothing to plot
  switch(length(active),   # length 0 -> "0"; length 1 -> "1"; >1 -> falls through
         "0" = {
           warning("No plot produced since all coefficients zero")
           return(invisible(NULL))
         },
         "1" = warning("1 or less nonzero coefficients; classo plot is not meaningful")
  )
  
  # Subset to active variables; preserve original indices for labelling
  beta_sub <- as.matrix(beta[active, , drop = FALSE])
  
  xvar <- match.arg(xvar)
  switch(xvar,
         "norm" = {
           index    <- if (missing(norm)) apply(Mod(beta_sub), 2, sum) else norm
           iname    <- "L1 Norm (Modulus)"
           approx.f <- 1
         },
         "lambda" = {
           index    <- log(lambda)
           iname    <- expression(Log(lambda))
           approx.f <- 0
         },
         "dev" = {
           index    <- dev
           iname    <- "Fraction Deviance Explained"
           approx.f <- 1
         }
  )
  
  # Default type to "l" (line) if not supplied by caller
  dotlist <- list(...)
  if (is.null(dotlist$type)) {
    dotlist$type <- "l"
  }
  
  # Save and restore graphical parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
  
  # Shared helper: add Df tick labels on top axis
  add_df_axis <- function() {
    atdf    <- pretty(index)
    prettydf <- approx(x = index, y = df, xout = atdf, rule = 2,
                       method = "constant", f = approx.f)$y
    axis(3, at = atdf, labels = prettydf, tcl = NA)
  }
  
  # Shared helper: add variable labels at end of path
  add_labels <- function(coef_panel) {
    xpos <- if (xvar == "lambda") min(index) else max(index)
    pos  <- if (xvar == "lambda") 2L else 4L
    ypos <- abs(coef_panel[, ncol(coef_panel)])
    text(rep(xpos, length(active)), ypos, paste(active), cex = 0.5, pos = pos)
  }
  
  # Real part 
  do.call(matplot, c(list(x = index, y = t(Re(beta_sub)),
                          lty = 1, xlab = xlab,
                          ylab = "Coefficients (Re)"),
                     dotlist))
  add_df_axis()
  if (label) add_labels(Re(beta_sub))
  
  # Imaginary part 
  do.call(matplot, c(list(x = index, y = t(Im(beta_sub)),
                          lty = 1, xlab = xlab,
                          ylab = "Coefficients (Im)"),
                     dotlist))
  add_df_axis()
  if (label) add_labels(Im(beta_sub))
  
  invisible(NULL)
}