#' Plot heatmap from a "cglasso" object
#'
#' Produces a heatmap of the estimated spectral precision matrix for a fitted
#' \code{"cglasso"} object. An inverse spectral matrix heatmap is produced.
#'
#' @aliases plot.cglasso
#' @param x Fitted \code{"cglasso"} model, i.e. an object of class
#'   \code{"cglassofit"} returned by \code{\link{cglasso}}.
#' @param index Integer index selecting which element of \code{x$Theta_list}
#'   to plot (i.e. which lambda value). Must be in
#'   \eqn{1, \ldots,} \code{length(x$Theta\_list)}. Default is \code{1}.
#' @param type Which component to plot: \code{"real"} for the real part,
#'   \code{"imaginary"} for the imaginary part, \code{"mod"} for the modulus
#'   (absolute value), or \code{"both"} for real and imaginary side by side.
#'   Default is \code{"mod"}.
#' @param label If \code{TRUE}, label the axes with variable names.
#' @param \dots Other graphical parameters passed to \code{image.plot}.
#'
#' @return No return value; called for its side effect of producing a plot.
#'
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{ykim124@ua.edu}
#' @seealso \code{\link{cglasso}}
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom fields image.plot
#' @method plot cglasso
#' @export
plot.cglasso <- function(x, index = 1,
                         type = c("mod", "real", "imaginary", "both"),
                         label = FALSE, ...) {
  
  type <- match.arg(type)
  
  # Validate index
  nlams <- length(x$Theta_list)
  if (!is.numeric(index) || length(index) != 1L ||
      index != round(index) || index < 1 || index > nlams) {
    stop(paste0("'index' must be a single integer in 1:", nlams, "."))
  }
  
  S <- x$Theta_list[[index]]
  p <- nrow(S)
  
  # Grid breaks for image.plot (one more than number of cells)
  x_coords <- seq(0, 1, length.out = p + 1)
  y_coords <- seq(0, 1, length.out = p + 1)
  
  # Mid-points for axis labels
  x_mids <- (x_coords[-length(x_coords)] + x_coords[-1]) / 2
  y_mids <- (y_coords[-length(y_coords)] + y_coords[-1]) / 2
  
  vnames <- if (label) colnames(S) else NULL
  
  # Helper: rotate matrix so row 1 appears at top-left in image.plot
  # image(x, y, z): z[i,j] is drawn at (x[i], y[j]), y increases bottom-to-top.
  # To display S[1,1] at top-left: transpose then reverse rows -> t(S)[, p:1]
  rotate <- function(M) t(M)[, p:1, drop = FALSE]
  
  # Helper: draw axis labels
  add_axes <- function() {
    if (!is.null(vnames)) {
      axis(side = 2, at = rev(y_mids), labels = vnames, las = 1)
      axis(side = 1, at = x_mids,      labels = vnames, las = 1)
    }
  }
  
  if (type != "both") {
    
    if (type == "real") {
      palette <- colorRampPalette(c("white", "blue"))
      mat     <- Re(S)
    } else if (type == "imaginary") {
      palette <- colorRampPalette(c("white", "red"))
      mat     <- Im(S)
    } else {
      palette <- colorRampPalette(c("white", "black"))
      mat     <- Mod(S)
    }
    
    z_scale <- range(mat)
    image.plot(x = x_coords, y = y_coords, z = rotate(mat),
               col = palette(100), zlim = z_scale,
               xaxt = "n", yaxt = "n", xlab = "", ylab = "",
               legend.shrink = 0.6, ...)
    add_axes()
    
  } else {
    
    palette_re <- colorRampPalette(c("white", "blue"))
    palette_im <- colorRampPalette(c("white", "red"))
    S_re <- Re(S)
    S_im <- Im(S)
    z_re <- range(S_re)
    z_im <- range(S_im)
    
    oldpar <- par(mfrow = c(2, 1))
    on.exit(par(oldpar), add = TRUE)
    
    image.plot(x = x_coords, y = y_coords, z = rotate(S_re),
               col = palette_re(100), zlim = z_re,
               xaxt = "n", yaxt = "n", xlab = "", ylab = "",
               legend.shrink = 0.6, ...)
    add_axes()
    
    image.plot(x = x_coords, y = y_coords, z = rotate(S_im),
               col = palette_im(100), zlim = z_im,
               xaxt = "n", yaxt = "n", xlab = "", ylab = "",
               legend.shrink = 0.6, ...)
    add_axes()
  }
  
  invisible(NULL)
}