#' @title QQ-plot
#'
#' @description Function produces QQ-plots for one or more model fits
#'
#' @param object_1 Object which contains model fit
#' @param object_2 Object which contains model fit
#' @param ... Other arguments.
#'
#' @details If \code{object_2} is specified, a QQ-plot for \code{object_1} and one for \code{object_2} will be
#' combined in the same plot.
#'
#' @seealso \code{\link{plot_qq.lm}} and \code{\link{plot_qq.lmvar}}
#'
#' @export
#'

plot_qq <- function( object_1, object_2 = NULL, ...){
  UseMethod( "plot_qq")
}
