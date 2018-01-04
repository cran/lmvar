#' @title Variance-covarience matrix of the coefficients beta for an object of class 'lmvar'
#'
#' @description Variance-covarience matrix (also simply called the 'covariance matrix') for the
#' maximum-likelihood estimators of \eqn{\beta_\mu} and \eqn{\beta_\sigma}.
#' The matrix is calculated with the assumption of asymptotic normality of maximum likelihood estimators.
#' This assumption is only valid in the limit of a large number of observations.
#'
#' @param object Object of class 'lmvar'
#' @param mu Specifies whether or not the covariance matrix for \eqn{\beta_\mu} is included in the returned matrix
#' @param sigma Specifies whether or not the  covariance matrix for \eqn{\beta_\sigma} is included in the returned matrix
#' @param ... For compatibility with \code{\link[stats]{vcov}} generic
#'
#' @return A 'matrix' object containing the (approximate) variance-covariance matrix of the maximum-likelihood estimators
#' of \eqn{\beta_\mu} and \eqn{\beta_\sigma} in \code{object}.
#'
#' @details The variance-covariance matrix is calculated as \eqn{I^{-1} / n} where \eqn{I} is the Fisher
#' information matrix and \eqn{n} the number of observations.
#'
#' When \code{mu = TRUE} and \code{sigma = TRUE}, the full covariance matrix for the combined vector
#' \eqn{(\beta_\mu, \beta_\sigma)} is returned.
#'
#' When \code{mu = TRUE} and \code{sigma = FALSE}, only the covariance matrix for \eqn{\beta_\mu} is returned.
#'
#' When \code{mu = FALSE} and \code{sigma = TRUE}, only the covariance matrix for \eqn{\beta_\sigma} is returned.
#'
#' @export
#'
#' @seealso \code{\link{summary.lmvar}} for standard errors for \eqn{\beta_\mu} and \eqn{\beta_\mu}.
#'
#' \code{\link{nobs.lmvar_no_fit}} for the number of observations in an object of class 'lmvar'.
#'
#' \code{\link{fisher}} for the Fisher information matrix of an object of class 'lmvar'.
#'
#' See the vignette "Math" (to be viewed with \code{vignette("Math", "lmvar")}) for details.
#'

vcov.lmvar <- function(object, mu = TRUE, sigma = TRUE, ...){

  # Calculate the inverses per block-matrix for efficiency

  if (mu){
    S_mu = fisher( object, sigma = FALSE)
    S_mu = chol2inv(chol(S_mu)) / nobs(object)
  }

  if (sigma){
    S_sigma = fisher( object, mu = FALSE)
    S_sigma = chol2inv(chol(S_sigma)) / nobs(object)
  }

  # Create block-matrix of zeros
  if (mu & sigma){
    S_0 = matrix( 0, nrow = nrow(S_mu), ncol = ncol(S_sigma))
  }

  # Set row and column names
  beta_mu_names = names( stats::coef(object, sigma = FALSE))
  beta_sigma_names = names( stats::coef(object, mu = FALSE))

  # Compose matrix from block-matrices
  if (mu & !sigma){
    rownames(S_mu) = beta_mu_names
    colnames(S_mu) = beta_mu_names
    return(S_mu)
  }
  else if (!mu & sigma){
    rownames(S_sigma) = beta_sigma_names
    colnames(S_sigma) = beta_sigma_names
    return(S_sigma)
  }
  else if (mu & sigma){
    S = rbind( cbind( S_mu, S_0), cbind( Matrix::t(S_0), S_sigma))
    beta_sigma_names = beta_sigma_names( beta_mu_names, beta_sigma_names)
    beta_names = c( beta_mu_names, beta_sigma_names)
    rownames(S) = beta_names
    colnames(S) = beta_names
    return(S)
  }
  else{
    return(matrix( 0, nrow = 0, ncol = 0))
  }
}
