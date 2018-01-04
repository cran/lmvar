#'
#' @title Aliased coefficients in an 'lmvar' object
#'
#' @description Returns the columns present in the user-specified model-matrices \eqn{X_\mu} and \eqn{X_\sigma} that were removed by
#' \code{lmvar} to make the matrices full-rank.
#'
#' @param object Object of class 'lmvar_no_fit' (hence it can also be of class 'lmvar')
#' @param mu Boolean, specifies whether the aliased columns from the model matrix \eqn{X_\mu} must be returned
#' @param sigma Boolean, specifies whether the aliased columns from the model matrix \eqn{X_\sigma} must be returned
#' @param ... Additional arguments, not used in the current implementation
#'
#' @return A character vector containing the names of the aliased columns
#'
#' @details If \code{mu = TRUE} and \code{sigma = TRUE}, the function returns the aliased columns of both \eqn{X_\mu}
#' and \eqn{X_\sigma}. The string "_s" is appended to the aliased column names from \eqn{X_\sigma} if at least one of those
#' names also appears in \eqn{X_\mu}
#'
#'If \code{mu = TRUE} and \code{sigma = FALSE}, the function returns the aliased columns of \eqn{X_\mu}.
#'
#'If \code{mu = FALSE} and \code{sigma = TRUE}, the function returns the aliased columns of \eqn{X_\sigma}.
#'
#' @export
#'
#' @importFrom stats alias
#'
#' @example R/examples/alias_examples.R
#'
alias.lmvar_no_fit <- function( object, mu = TRUE, sigma = TRUE, ...){

  aliased_mu = character()
  if (mu){
    aliased_mu = names(object$aliased_mu)[object$aliased_mu]
  }

  aliased_sigma = character()
  if (sigma){
    aliased_sigma = names(object$aliased_sigma)[object$aliased_sigma]
  }

  if (mu & !sigma){
    return(aliased_mu)
  }
  else if (!mu & sigma){
    return(aliased_sigma)
  }
  else if (mu & sigma){
    a_names = c( aliased_mu, beta_sigma_names( names(object$aliased_mu), aliased_sigma))
    return(as.character(a_names))
  }
  else {
    return(character())
  }
}
