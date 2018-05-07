#' @title Plot of the distribution of quantiles
#'
#' @description Function produces plot of the distribution of quantiles for one or more model fits
#'
#' @param object_1 Object which contains model fit
#' @param object_2 Object which contains model fit
#' @param ... Other arguments.
#'
#' @details If \code{object_2} is specified, a plot for the distribution of quantiles for \code{object}
#' and one for \code{object_2} will be combined in the same plot.
#'
#' @seealso \code{\link{plot_qdis.lm}} and \code{\link{plot_qdis.lmvar}}
#'
#' @export
#'

plot_qdis <- function( object_1, object_2 = NULL, ...){
  UseMethod( "plot_qdis")
}
