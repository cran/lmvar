#'
#' @title Extracts coefficients from an 'lmvar' object.
#'
#' @description Extracts maximum-likelihood estimators for \eqn{\beta_\mu} and \eqn{\beta_\sigma} from an 'lmvar' object.
#'
#' @param object Object of class 'lmvar'
#' @param mu Boolean, specifies whether or not to return the maximum-likelihood estimator for \eqn{\beta_\mu}
#' @param sigma Boolean, specifies whether or not to return the maximum-likelihood estimator for \eqn{\beta_\sigma}
#' @param ... For compatibility with \code{\link[stats]{coef}} generic
#'
#' @return When \code{mu = TRUE} and \code{sigma = TRUE}, a named numeric vector with the elements of \eqn{\beta_\mu},
#' followed by the elements of \eqn{\beta_\sigma}.
#'
#' When \code{mu = TRUE} and \code{sigma = FALSE}, a named numeric vector with the elements of \eqn{\beta_\mu}.
#'
#' When \code{mu = FALSE} and \code{sigma = TRUE}, a named numeric vector with the elements of \eqn{\beta_\sigma}.
#'
#' @details When both \code{mu = TRUE} and \code{sigma = TRUE}, the names of the
#' coefficients in \eqn{\beta_\sigma} are adapted to distinguish them from the names in \eqn{\beta_\mu}, if needed.
#'
#' @seealso \code{\link{beta_sigma_names}} for the adaptation of the names of the coefficients in \eqn{\beta_\sigma}.
#'
#' \code{\link[stats]{confint}} for the calculation of confidence intervals of \eqn{\beta_\mu} and \eqn{\beta_\sigma}.
#'
#' @export
#'
#' @example R/examples/coef_examples.R
#'

coef.lmvar <- function( object, mu = TRUE, sigma = TRUE, ...){

  beta_mu = numeric()
  beta_sigma = numeric()
  beta_mu_names = character()
  beta_sigma_names = character()

  if (mu){
    beta_mu = object$coefficients_mu
    beta_mu_names = names(object$coefficients_mu)
  }

  if (sigma){
    beta_sigma = object$coefficients_sigma
    if (mu){
      beta_sigma_names = lmvar::beta_sigma_names( names(object$coefficients_mu), names(object$coefficients_sigma))
    }
    else {
      beta_sigma_names = names(object$coefficients_sigma)
    }
  }

  betas = c( beta_mu, beta_sigma)
  names(betas) = c( beta_mu_names, beta_sigma_names)

  return(betas)
}
