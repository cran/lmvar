#' @title Residual degrees of freedom for an object of class 'lmvar'
#'
#' @description Residual degrees of freedom for an object of class 'lmvar'. The residual degrees of freedom are defined
#' as the number of observations minus the degrees of freedom of the model.
#'
#' @param object Object of class 'lmvar'
#' @param ... For compatibility with \code{\link[stats]{df.residual}} generic
#'
#' @return Residual degees of freedom for \code{object}.
#'
#' @export
#'
#' @seealso \code{\link{dfree}} for the degrees of freedom of an object of class 'lmvar'
#'
#' \code{\link{nobs.lmvar}} for the number of observations in an object of class 'lmvar'
#'
#' @example R/examples/df.residual_examples.R
#'

df.residual.lmvar <- function(object, ...){
  df = nobs(object) - dfree(object)
  return(df)
}
