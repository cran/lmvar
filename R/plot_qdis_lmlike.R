#' @title Plot of the distribution of quantiles for objects of class 'lm' or 'lmvar'
#'
#' @description Function produces a histogram of quantiles for  objects of
#' class 'lm' or 'lmvar'. This function is called by \code{\link{plot_qdis.lm}} and
#' \code{\link{plot_qdis.lmvar}}. It is not intended to be called directly.
#'
#' @param object_1 Object of class 'lm' or 'lmvar'
#' @param object_2 Object of class 'lm' or class 'lmvar'
#' @param name_1 Character string, the name of \code{object_1}
#' @param name_2 Character string, the name of \code{object_2}
#'
#' @details If \code{object_2} is specified, a plot for \code{object_1} and one for \code{object_2} will be
#' combined in the same plot.
#'
#' The string \code{name_1} and (optionally) \code{name_2} are used in the legend of the plot as names for
#' \code{object_1} and (optionally) \code{object_2}.
#'
#' All inputs of class 'lm' must contain the response vector \eqn{y}. I.e., one must run \code{\link[stats]{lm}}
#' with argument \code{y = TRUE}.
#'

plot_qdis_lmlike <- function( object_1, object_2 = NULL, name_1, name_2 = NULL){

  # check input
  if (!("y" %in% names(object_1))){
    stop("Response vector y must be in object_1. Please run lm with y = TRUE.")
  }
  if (!is.null(object_2) & !("y" %in% names(object_2))){
    stop("Response vector y must be in object_2. Please run lm with y = TRUE.")
  }

  # set colors
  col = list("black", grDevices::rgb( 0, 0.5, 0, 0.3))

  if (inherits( object_1, 'lm')){
    mu = stats::fitted(object_1)
    sigma = summary(object_1)$sigma
  }
  else {
    mu = fitted.lmvar( object_1, sigma = FALSE)
    sigma = fitted.lmvar( object_1, mu = FALSE)
  }
  quant = stats::pnorm( object_1$y, mean = mu, sd = sigma)

  legend_text = name_1
  legend_col = col[[1]]

  # postpone plot if there is a second object (to determine maximum y-value)
  args = list( quant, main = "Distribution of quantiles", xlab = "Quantile", freq = FALSE,
               border = col[[1]], density = 10)

  if (is.null(object_2)){
    p = do.call(graphics::hist, args)
  }
  else {
    p = graphics::hist( quant, plot = FALSE)
  }

  y_max = max(p$density)
  x_legend = 0
  y_legend = max(p$density)

  if (!is.null(object_2)){
    if (inherits( object_2, "lm")){
      mu = stats::fitted(object_2)
      sigma = summary(object_2)$sigma
    }
    else {
      mu = fitted.lmvar( object_2, sigma = FALSE)
      sigma = fitted.lmvar( object_2, mu = FALSE)
    }
    quant = stats::pnorm( object_2$y, mean = mu, sd = sigma)
    p = graphics::hist( quant, plot = FALSE)

    args$ylim = c( 0, max( y_max, p$density))
    y_legend = args$ylim[2]

    p = do.call( graphics::hist, args)
    graphics::hist( quant, freq = FALSE, add = TRUE, col = col[[2]], border = col[[2]])

    legend_text = c( name_1, name_2)
    legend_col = c( legend_col, col[[2]])
  }
  graphics::abline( 1, 0, col = "black", lty = 3)

  graphics::legend( x_legend, y_legend, legend = legend_text, col = legend_col, lty = c( 1, 1), bty = "n")
}
