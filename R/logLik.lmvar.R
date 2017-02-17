#' Log-likelihood for an object of class 'lmvar'
#'
#' @param object Object of class 'lmvar'
#' @param ... For compatibility with \code{\link[stats]{logLik}} generic
#'
#' @return 'logLik' object, a number containing the log-likelihood with an attribute 'df' containing the degrees of freedom
#'
#' @export
#'
#' @seealso \code{\link{dfree}} for the degrees of freedom for an object of class 'lmvar'.
#'
#' @example R/examples/logLik_examples.R
#'
logLik.lmvar <- function( object, ...){

  logl = object$logLik
  attr(logl, "df") = dfree(object)
  class(logl) = "logLik"

  return(logl)
}
