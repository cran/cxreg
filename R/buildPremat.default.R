#' Build Prediction Matrix
#'
#' These are not intended for use by users. Constructs a prediction matrix for
#' cross-validation folds using the coefficient paths in \code{outlist}.
#' Inspired by internal functions in the \pkg{glmnet} package.
#' @param outlist A list of fitted \code{classo} models, one per fold.
#' @param lambda A vector of penalty values from the master fit.
#' @param x The full predictor matrix (complex-valued).
#' @param foldid Integer vector of fold assignments, length \code{nrow(x)}.
#' @param alignment Alignment method: \code{"lambda"} (default) aligns fold
#'   predictions to the master lambda sequence; \code{"fraction"} aligns by
#'   fraction of path progress.
#' @param \dots Other arguments passed to \code{predict}.
#' @return A complex matrix of predicted values, \code{nrow(x)} by
#'   \code{length(lambda)}. Rows correspond to observations held out in each
#'   fold; all other entries are \code{NA}.
#' @export
buildPredmat <- function(outlist, lambda, x, foldid, alignment, ...) {
  UseMethod("buildPredmat")
}

#' @rdname buildPredmat
#' @export
buildPredmat.default <- function(outlist, lambda, x, foldid, alignment, ...) {
  
  nlambda <- length(lambda)
  nfolds  <- max(foldid)
  # Initialise as complex so the matrix can hold complex predictions
  predmat <- matrix(NA_complex_, nrow(x), nlambda)
  
  for (i in seq(nfolds)) {
    which  <- foldid == i
    fitobj <- outlist[[i]]
    
    preds <- switch(alignment,
                    lambda   = predict(fitobj, x[which, , drop = FALSE], s = lambda, ...),
                    fraction = predict(fitobj, x[which, , drop = FALSE], ...)
    )
    
    nlami <- min(ncol(preds), nlambda)
    predmat[which, seq(nlami)] <- preds[, seq(nlami)]
    
    # Fill remaining columns with last available prediction (forward-fill)
    if (nlami < nlambda) {
      predmat[which, seq(from = nlami, to = nlambda)] <- preds[, nlami]
    }
  }
  
  rn <- rownames(x)
  sn <- paste0("s", seq(0, length.out = nlambda))
  dimnames(predmat) <- list(rn, sn)
  predmat
}