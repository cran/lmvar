#' @title Fisher information matrix for an object of class 'lmvar'
#'
#' @description Fisher information matrix for
#' an object of class 'lmvar'.
#'
#' @param object Object of class 'lmvar'
#' @param mu Specifies whether or not the block-matrix for \eqn{\beta_\mu} is included in the returned matrix
#' @param sigma Specifies whether or not the block-matrix for \eqn{\beta_\sigma} is included in the returned matrix
#' @param ... Additional arguments, not used in the current implementation
#'
#' @return An object of class 'matrix' containing the Fisher information matrix of \code{object}.
#'
#' @details The Fisher information matrix is calculated as minus \eqn{-E[H]/n} with \eqn{E[H]} the expected value of
#' the Hessian matrix \eqn{H} of the
#' log-likelihood and \eqn{n} the number of observations.
#'
#' The matrix is calculated using the maximum-likelihood estimators of \eqn{\mu} and \eqn{\sigma}.
#'
#' If \code{mu = TRUE} and \code{sigma = TRUE}, the full Fisher information matrix is returned.
#'
#' If \code{mu = TRUE} and \code{sigma = FALSE}, only the left-upper block-matrix is returned, corresponding to the part of
#' the Fisher information matrix pertaining to \eqn{\beta_\mu}.
#'
#' If \code{mu = FALSE} and \code{sigma = TRUE}, only the right-lower block-matrix is returned, corresponding to the part of
#' the Fisher information matrix pertaining to \eqn{\beta_\sigma}.
#'
#' @export
#'
#' @seealso \code{\link{vcov.lmvar}} calculates the covariance matrix for the maximum-likelihood estimators of
#' \eqn{\beta_\mu} and \eqn{\beta_\mu}
#'
#' \code{\link{nobs.lmvar}} for the number of observations in an object of class 'lmvar'
#'
#' \code{\link{coef.lmvar}} for the coefficients \eqn{\beta_\mu} and \eqn{\beta_\sigma}
#'
#' \code{\link{fitted.lmvar}} for the expectation values \eqn{\mu} and standard deviations \eqn{\sigma}.
#'
#' See the vignette "Math" (to be viewed with \code{vignette("Math", "lmvar")}) for details.
#'
#' @example R/examples/fisher_examples.R
#'
fisher <- function(object, mu = TRUE, sigma = TRUE, ...){

  X_sigma = object$X_sigma

  # Calculate block-matrix for mu
  if (mu){
    X = object$X_mu
    sigma_object = as.numeric(exp(X_sigma %*% stats::coef(object, mu=FALSE)))
    I_mu = X * (1/sigma_object)
    I_mu = (Matrix::t(I_mu) %*% I_mu) / nobs(object)
  }

  # Calculate block-matrix for sigma
  if (sigma){
    I_sigma = 2 * (Matrix::t(X_sigma) %*%  X_sigma) / nobs(object)
  }

  # Create block-matrix of zeros
  if (mu & sigma){
    H_0 = matrix( 0, nrow = nrow(I_mu), ncol = ncol(I_sigma))
  }

  # Set row and column names
  beta_mu_names = names( stats::coef(object, sigma = FALSE))
  beta_sigma_names = names( stats::coef(object, mu = FALSE))

  # Compose matrix from block-matrices
  if (mu & !sigma){
    rownames(I_mu) = beta_mu_names
    colnames(I_mu) = beta_mu_names
    return(I_mu)
  }
  else if (!mu & sigma){
    rownames(I_sigma) = beta_sigma_names
    colnames(I_sigma) = beta_sigma_names
    return(I_sigma)
  }
  else if (mu & sigma){
    I = rbind( cbind( I_mu, H_0), cbind( Matrix::t(H_0), I_sigma))
    beta_sigma_names = beta_sigma_names( beta_mu_names, beta_sigma_names)
    beta_names = c( beta_mu_names, beta_sigma_names)
    rownames(I) = beta_names
    colnames(I) = beta_names
    return(I)
  }
  else{
    return(matrix( 0, nrow = 0, ncol = 0))
  }
}
