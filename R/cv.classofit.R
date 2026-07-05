cv.classofit <- function(predmat, y, type.measure, weights, foldid) {
  
  y <- drop(y)
  N <- nrow(predmat) - apply(is.na(predmat), 2, sum)
  
  cvraw <- switch(type.measure,
                  mse      = Mod(y - predmat)^2,
                  deviance = Mod(y - predmat)^2
  )
  
  list(cvraw = cvraw, weights = weights, N = N, type.measure = type.measure)
}