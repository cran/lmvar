#' @title Number of observations for an object of class 'lmvar'
#'
#' @description The number of observations in an object of class 'lmvar'.
#'
#' @param object Object of class 'lmvar'
#' @param ... For compatibility with \code{\link[stats]{nobs}} generic
#'
#' @return Integer containing the number of observations in the model in \code{object}.
#'
#' @importFrom stats nobs
#'
#' @export
#'
#' @example R/examples/nobs_examples.R

nobs.lmvar <- function( object, ...){
  return(length(object$y))
}
