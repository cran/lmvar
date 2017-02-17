#' @title AIC for an object of class 'lmvar'
#'
#' @description AIC (Aikaike's 'An Information Criterion') for an object of class 'lmvar'
#'
#' @param object Object of class 'lmvar'
#' @param ... For compatibility with \code{\link[stats]{AIC}} generic
#' @param k Numeric, the penalty per parameter to be used. The default k = 2 is the classical AIC.
#'
#' @return the AIC of the object
#'
#' @export
#'
#' @example R/examples/AIC_examples.R
#'
AIC.lmvar <- function(object, ..., k = 2){
  df = dfree(object)
  return(k*df - 2 * logLik.lmvar(object)[1])
}
