#' Find nonzero coefficients in a coefficient matrix
#'
#' Internal function used by \code{\link{predict.classo}} to identify which
#' variables have nonzero coefficients, either across the entire path or at
#' each lambda step individually.
#'
#' @param beta A coefficient matrix of dimension \code{nvars x nlambda}.
#'   May be complex-valued; a coefficient is considered nonzero when
#'   \code{abs()} (equivalently \code{Mod()} for complex values) is positive.
#' @param bystep Logical. If \code{FALSE} (default), returns the indices of
#'   variables that are nonzero at \emph{any} lambda. If \code{TRUE}, returns
#'   a named list of length \code{nlambda}, each element giving the integer
#'   indices of nonzero variables at that lambda step (or \code{integer(0)}
#'   if all are zero).
#' @return When \code{bystep = FALSE}, an integer vector of variable indices
#'   (or \code{NULL} if all are zero). When \code{bystep = TRUE}, a named
#'   list of length \code{nlambda}.
#' @keywords internal
nonzeroCoef <- function(beta, bystep = FALSE) {
  
  nr <- nrow(beta)
  
  if (nr == 1) {
    # Degenerate case: single variable
    if (bystep) {
      lapply(seq_len(ncol(beta)), function(k) {
        if (abs(beta[1, k]) > 0) 1L else integer(0)
      })
    } else {
      if (any(abs(beta) > 0)) 1L else NULL
    }
    
  } else {
    
    # Logical matrix: TRUE where coefficient is nonzero
    beta_nz <- abs(beta) > 0
    idx     <- seq_len(nr)
    ones    <- rep(1, ncol(beta))
    
    # Variables that are nonzero in at least one lambda step
    nz  <- as.vector((beta_nz %*% ones) > 0)
    idx <- idx[nz]
    
    if (bystep) {
      # Always return a named list — one element per lambda step.
      # Using lapply avoids apply()'s unpredictable return type.
      dn <- dimnames(beta)[[2]]
      if (length(idx) > 0) {
        beta_nz_sub <- as.matrix(beta_nz[idx, , drop = FALSE])
        result <- lapply(seq_len(ncol(beta_nz_sub)), function(k) {
          hits <- beta_nz_sub[, k]
          if (any(hits)) idx[hits] else integer(0)
        })
      } else {
        result <- lapply(seq_len(ncol(beta)), function(k) integer(0))
      }
      if (!is.null(dn)) names(result) <- dn
      result
      
    } else {
      idx
    }
  }
}