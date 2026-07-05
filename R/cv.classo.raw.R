#' @importFrom foreach foreach %dopar%
cv.classo.raw <- function(x, y, weights, lambda, type.measure, nfolds, foldid,
                          alignment, keep, parallel, trace.it, classo.call, cv.call, ...) {
  
  if (isTRUE(trace.it == 1)) {
    cat("Training\n")
  }
  
  classo.object <- classo(x, y,
                          weights     = weights,
                          lambda      = lambda,
                          standardize = TRUE,
                          intercept   = FALSE,
                          maxit       = 100000,
                          trace.it    = 0, ...)
  classo.object$call <- classo.call
  subclass <- class(classo.object)[[1]]

  # cvtype.R
  type.measure <- cvtype(type.measure, subclass)
  
  lambda  <- classo.object$lambda
  nz      <- sapply(predict(classo.object, type = "nonzero"), length)
  outlist <- as.list(seq(nfolds))
  N       <- nrow(x)

  if (parallel) {
    
    # Guard: require a registered parallel backend
    if (foreach::getDoParWorkers() == 1) {
      warning(paste(
        "parallel=TRUE but no parallel backend is registered.",
        "Register one first, e.g. doParallel::registerDoParallel().",
        "Falling back to sequential execution."
      ))
      parallel <- FALSE
    }
  }
  
  if (parallel) {
    
    outlist <- foreach(i = seq(nfolds), .packages = "cxreg") %dopar% {
      which   <- foldid == i
      y_sub   <- if (length(dim(y)) > 1) y[!which, ] else y[!which]
      x_sub   <- x[!which, , drop = FALSE]
      w_sub   <- weights[!which]
      classo(x           = x_sub,
             y           = y_sub,
             weights     = w_sub,
             lambda      = lambda,
             standardize = TRUE,
             intercept   = FALSE,
             maxit       = 100000,
             trace.it    = 0)
    }
    
  } else {
    
    for (i in seq(nfolds)) {
      if (isTRUE(trace.it == 1)) {
        cat(sprintf("Fold: %d/%d\n", i, nfolds))
      }
      which   <- foldid == i
      y_sub   <- if (length(dim(y)) > 1) y[!which, ] else y[!which]
      x_sub   <- x[!which, , drop = FALSE]
      w_sub   <- weights[!which]
      outlist[[i]] <- classo(x           = x_sub,
                             y           = y_sub,
                             weights     = w_sub,
                             lambda      = lambda,
                             standardize = TRUE,
                             intercept   = FALSE,
                             maxit       = 100000,
                             trace.it    = 0)
    }
  }
  
  # buildPredmat.default.R
  class(outlist) <- paste0(subclass, "list")
  predmat <- buildPredmat(outlist, lambda, x, foldid, alignment)
  
  # Compute the cross-validated measure for classofit objects
  cvstuff <- cv.classofit(predmat, y, type.measure, weights, foldid)
  
  # cvstats.R
  out    <- cvstats(cvstuff, foldid, nfolds, lambda)
  cvname <- names(cvstuff$type.measure)
  names(cvname) <- cvstuff$type.measure
  out <- c(out, list(
    call       = cv.call,
    name       = cvname,
    classo.fit = classo.object,
    nzero      = nz
  ))
  
  if (keep) {
    out <- c(out, list(fit.preval = predmat, foldid = foldid))
  }
  
  # getOptcv.classo.R
  lamin <- with(out, getOptcv.classo(lambda, cvm, cvsd, cvname))
  obj   <- c(out, as.list(lamin))
  class(obj) <- "cv.classo"
  obj
}