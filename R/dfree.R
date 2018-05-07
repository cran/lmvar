#' @title Degrees of freedom for an object of class 'lmvar'
#'
#' @description Degrees of freedom for the model in an object of class 'lmvar'. The degrees of freedom are defined as the rank of the
#' model matrix \eqn{X_\mu} for the expectation values, plus the rank of the model matrix \eqn{X_\sigma} for the standard deviations.
#'
#' @param object Object of class 'lmvar_no_fit' (hence it can also be of class 'lmvar')
#' @param mu Boolean, specifies whether the degrees of freedom for the model for the expectation values must be included.
#' @param sigma Boolean, specifies whether the degrees of freedom for the model for the standard deviations must be included.
#' @param ... Additional arguments, not used in the current implementation
#'
#' @return An integer containing the degrees of freedom for the model in \code{object}.
#'
#' @details If \code{mu = TRUE} and \code{sigma = TRUE}, the function returns the rank of the model-matrix \eqn{X_\mu} plus the
#' rank of the model matrix \eqn{X_\sigma}.
#'
#' If \code{mu = TRUE} and \code{sigma = FALSE}, the function returns the rank of the model-matrix \eqn{X_\mu}.
#'
#' If \code{mu = FALSE} and \code{sigma = TRUE}, the function returns the rank of the model-matrix \eqn{X_\sigma}.
#'
#' Both model matrices contain a column corresponding to an intercept term. This column is added by \code{\link{lmvar}}.
#' See also the vignette 'Intro'.
#'
#' @export
#'
#' @example R/examples/dfree_examples.R
#'

dfree <- function( object, mu = TRUE, sigma = TRUE, ...){

  if (!(inherits( object, 'lmvar_no_fit'))){
    stop("Object must be an 'lmvar_no_fit' (or 'lmvar') object")
  }

  if (mu & sigma){
    return( ncol(object$X_mu) + ncol(object$X_sigma))
  }
  else if (mu){
    return( ncol(object$X_mu))
  }
  else {
    return( ncol(object$X_sigma))
  }

}
