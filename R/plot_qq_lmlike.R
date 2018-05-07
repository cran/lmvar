#' @title QQ-plot for an object of class 'lm' or class 'lmvar'
#'
#' @description Function produces QQ-plots for an object of class 'lm' or 'lmvar'. This function is
#' called by \code{\link{plot_qq.lm}} and
#' \code{\link{plot_qq.lmvar}}. It is not intended to be called directly.
#'
#' @param object_1 Object of class 'lm' or class 'lmvar'
#' @param object_2 Object of class 'lm' or class 'lmvar'
#' @param name_1 Character string, the name of \code{object_1}
#' @param name_2 Character string, the name of \code{object_2}
#'
#' @details If \code{object_2} is specified, a QQ-plot for \code{object_1} and one for \code{object_2} will be
#' combined in the same plot.
#'
#' The string \code{name_1} and (optionally) \code{name_2} are used in the legend of the plot as names for
#' \code{object_1} and (optionally) \code{object_2}.
#'
#' All inputs of class 'lm' must contain the response vector \eqn{y}. I.e., one must run \code{\link[stats]{lm}}
#' with argument \code{y = TRUE}.
#'


plot_qq_lmlike <- function( object_1, object_2 = NULL, name_1, name_2 = NULL){

  # check input
  if (!("y" %in% names(object_1))){
    stop("Response vector y must be in object_1. Please run lm with y = TRUE.")
  }
  if (!is.null(object_2) & !("y" %in% names(object_2))){
    stop("Response vector y must be in object_2. Please run lm with y = TRUE.")
  }

  # set colors
  col = list("black", grDevices::rgb( 0, 0.5, 0, 0.3))

  if (inherits( object_1, "lm")){
    mu = stats::fitted(object_1)
    sigma = summary(object_1)$sigma
  }
  else {
    mu = fitted.lmvar( object_1, sigma = FALSE)
    sigma = fitted.lmvar( object_1, mu = FALSE)
  }
  z = (object_1$y - mu) / sigma

  legend_text = name_1
  legend_col = col[[1]]

  # postpone plot if there is a second object
  args = list( z, main = "QQ plot", xlab = "z-scores from sample quantile", ylab = "z-scores from fit",
               cex.lab = 1.2, col = col[[1]], pch = 20, cex = 0.5)
  if (is.null(object_2)){
    p = do.call( stats::qqnorm, args)
  }
  else {
    p = stats::qqnorm( z, plot.it = FALSE)
  }

  x_min = min(p$x)
  x_max = max(p$x)
  y_min = min(p$y)
  y_max = max(p$y)
  x_legend = 1 * x_min + 0 * x_max
  y_legend = 0 * y_min + 1 * y_max

  if (!is.null(object_2)){
    if (inherits( object_2, "lm")){
      mu = stats::fitted(object_2)
      sigma = summary(object_2)$sigma
    }
    else {
      mu = fitted.lmvar( object_2, sigma = FALSE)
      sigma = fitted.lmvar( object_2, mu = FALSE)
    }
    z = (object_2$y - mu) / sigma
    p = stats::qqnorm( z, plot.it = FALSE)

    # set minimum and maximum of axes
    if (min(p$x) < x_min | max(p$x) > x_max){
      x_min = min( x_min, p$x)
      x_max = max( x_max, p$x)
      d = x_max - x_min
      x_min = x_min - 0.05 * d
      x_max = x_max + 0.05 * d
      args$xlim = c( x_min, x_max)
    }
    if (min(p$y) < y_min | max(p$y) > y_max){
      y_min = min( y_min, p$y)
      y_max = max( y_max, p$y)
      d = y_max - y_min
      y_min = y_min - 0.05 * d
      y_max = y_max + 0.05 * d
      args$ylim = c( y_min, y_max)
    }

    do.call( stats::qqnorm, args)
    graphics::points( p$x, p$y, col = col[[2]], pch = 20, cex = 0.5)

    legend_text = c( name_1, name_2)
    legend_col = c( legend_col, col[[2]])
  }
  graphics::abline( 0, 1, col = "red")

  graphics::legend( x_legend, y_legend, legend = legend_text, col = legend_col, lty = c( 3, 3), box.lty = 0)
}
