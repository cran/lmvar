#' @title Residuals from an 'lmvar' object
#'
#' @description Calculates residuals from an 'lmvar' object. This object can be a fit to either a response vector or the
#' logarithm of the response vector.
#'
#' @param object Object of class 'lmvar'
#' @param log Boolean, specifies whether \code{object} is a fit to a response-variable \eqn{Y} or to its logarithm \eqn{\log Y}
#' In both cases, \code{residuals.lmvar} returns residuals for \eqn{Y} itself.
#' @param ... For compatibility with \code{\link[stats]{residuals}} generic
#'
#' @return A numeric vector with the residual for each observation in \code{object}.
#'
#' @details In case \code{log = FALSE}, the residual of an observation is defined as \eqn{y - \mu}, where \eqn{y} is the value of the observation and \eqn{\mu} its expected
#' value.
#'
#' In case \code{log = TRUE}, the residual of an observation is defined as \eqn{e^y - \mu}, where \eqn{\mu} is the expected
#' value of \eqn{e^y}.
#'
#' @seealso \code{\link{fitted.lmvar}} for the expected values in an object of class 'lmvar'.
#'
#' @export
#'
#' @example R/examples/residuals_examples.R
#'

residuals.lmvar <- function( object, log = FALSE, ...){

  if (log){
    res = exp(object$y) - fitted.lmvar(object, sigma = FALSE, log = TRUE)
  }
  else {
    res = object$y - fitted.lmvar(object, sigma = FALSE)
  }

  return(res)
}
