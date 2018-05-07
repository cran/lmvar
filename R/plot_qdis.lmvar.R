#' @title Plot of the distribution of quantiles for an object of class 'lmvar'
#'
#' @description Function produces plot of the distribution of quantiles for an object of
#' class 'lmvar' and, optionally, for another object of class 'lm' or 'lmvar'.
#'
#' @param object_1 Object of class 'lmvar'
#' @param object_2 Object of class 'lm' or class 'lmvar'
#' @param ... for compatibility with \code{\link{plot_qdis}} generic.
#'
#' @details If \code{object_2} is specified, a plot for \code{object_1} and one for \code{object_2} will be
#' combined in the same plot.
#'
#' If \code{object_2} is of class 'lm', it must contain the response vector \eqn{y}. I.e., one must run \code{\link[stats]{lm}}
#' with argument \code{y = TRUE}.
#'
#' @seealso \code{\link{plot_qdis}}
#'
#' @export
#'
#' @example R/examples/plot_qdis.lmvar_examples.R
#'

plot_qdis.lmvar <- function( object_1, object_2 = NULL, ...){

  name_1 = deparse(substitute(object_1))
  name_2 = deparse(substitute(object_2))

  plot_qdis_lmlike( object_1, object_2, name_1, name_2)
}
