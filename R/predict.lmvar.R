#' @title Predictions for model matrices
#'
#' @description Estimators and confidence intervals for the expected values and standard deviations of the response vector \eqn{Y},
#' given model matrices \code{X_mu} and \code{X_sigma}. Prediction intervals for \eqn{Y} as well.
#' The estimators are based on the maximum likelihood-estimators for \eqn{\beta_\mu} and \eqn{\beta_\sigma}
#' present in an 'lmvar' object. Alternatively, estimators and intervals can be for \eqn{e^Y}.
#'
#' @param object Object of class 'lmvar'
#' @param X_mu Model matrix for the expected values
#' @param X_sigma Model matrix for the logarithm of the standard deviations
#' @param mu Boolean, specifies whether or not to include the predictions for the expected values
#' @param sigma Boolean, specifies whether or not to include the predictions for the standard deviations
#' @param log Boolean, specifies whether expected values, standard deviations (as well as their confidence intervals) and
#' prediction intervals should be for \eqn{Y} (\code{log = FALSE}) or for \eqn{e^Y} (\code{log = TRUE}).
#' @param interval Character string, specifying the type of interval. Possible values are
#' \itemize{
#' \item "none" No interval, this is the default
#' \item "confidence" Confidence intervals for the expected values (if \code{mu = TRUE}) and the standard deviation
#' (if \code{sigma = TRUE})
#' \item "prediction" Prediction intervals for the response vector \eqn{Y} (\code{log = FALSE}) or for \eqn{e^Y} (\code{log = TRUE})
#' }
#' @param level Numeric value between 0 and 1, specifying the confidence level
#' @param ... For compatibility with \code{\link[stats]{predict}} generic
#'
#' @return In the case \code{mu = FALSE} and \code{interval = "none"}: a numeric vector containing the estimators for
#' the standard deviation.
#'
#' In the case code{sigma = FALSE} and \code{interval = "none"}: a numeric vector containing the estimators for
#' the expected values.
#'
#' In all other cases: a matrix with one column for each requested feature and one row for each observation. The column names are
#' \itemize{
#' \item \code{mu} Estimators for the expected value \eqn{\mu}
#' \item \code{sigma} Estimators for the standard deviation \eqn{\sigma}
#' \item \code{mu_lwr} Lower bound of the confidence interval for \eqn{\mu}
#' \item \code{mu_upr} Upper bound of the confidence interval for \eqn{\mu}
#' \item \code{sigma_lwr} Lower bound of the confidence interval for \eqn{\sigma}
#' \item \code{sigma_upr} Upper bound of the confidence interval for \eqn{\sigma}
#' \item \code{lwr} Lower bound of the prediction interval
#' \item \code{upr} Upper bound of the prediction interval
#' }
#'
#' @export
#'
#' @details When \code{X_mu = NULL}, the model matrix \eqn{X_\mu} is taken from \code{object}. Likewise, when
#' \code{X_sigma = NULL}, \eqn{X_\sigma} is taken from \code{object}.
#'
#' Both \code{X_mu} and \code{X_sigma} must have column names. Column names are matched with the names of the elements of
#' \eqn{\beta_\mu} and \eqn{\beta_\sigma} in \code{object}. Columns with non-matching names are ignored. In case not all
#' names in \eqn{\beta_\mu} can be matched with a column in \code{X_mu}, a warning is given. The same is true for \eqn{\beta_\sigma}
#' and \code{X_sigma}.
#'
#' \code{X_mu} can not have a column with the name  "(Intercept)". This column is added by \code{predict.lmvar} in case
#' it is present in \code{object}. Likewise,
#' \code{X_sigma} can not have a column with the name  "(Intercept_s)". It is added by \code{predict.lmvar} in case
#' it is present in \code{object}
#'
#' Both matrices must be numeric and can not contain special values like
#' \code{NULL}, \code{NaN}, etc.
#'
#' If \code{log = FALSE}, \code{predict.lmvar} returns
#' expected values and standard deviations for the observations \eqn{Y} corresponding to the model matrices \eqn{X_\mu}
#' and \eqn{X_\sigma}.
#'
#' If \code{log = TRUE}, \code{predict.lmvar} returns expected values and standard deviations for \eqn{e^Y}.
#'
#' Confidence intervals are calculated under the asumption of asymptotic normality. This asumption holds when the number
#' of observations is large.
#' Intervals must be treated cautiously in case of a small number of observations.
#'
#' \code{predict.lmvar} with \code{X_mu = NULL} and \code{X_sigma = NULL} is equivalent to the function
#' \code{\link{fitted.lmvar}}.
#'
#' @seealso \code{\link{coef.lmvar}} and \code{\link[stats]{confint}} for maximum likelihood estimators and confidence intervals
#' for \eqn{\beta_\mu} and \eqn{\beta_\sigma}.
#'
#' \code{\link{fitted.lmvar}} is equivalent to \code{predict.lmvar} with \code{X_mu = NULL} and \code{X_sigma = NULL}.
#' I.e. \code{fitted.lmvar} gives estimators and confidence intervals for the expected values and standard deviations of
#' the response vector of the fit in \code{object}.
#'
#' @example R/examples/predict_examples.R
#'

predict.lmvar <- function( object, X_mu = NULL, X_sigma = NULL, mu = TRUE, sigma = TRUE, log = FALSE,
                           interval = c("none", "confidence", "prediction"), level = 0.95, ...){

  # Check column names
  if (is.element('(Intercept)', colnames(X_mu))){
    stop("Matrix X_mu is not allowed to have a column named '(Intercept)'")
  }
  if (is.element('(Intercept_s)', colnames(X_sigma))){
    stop("Matrix X_sigma is not allowed to have a column named '(Intercept_s)'")
  }

  # Check interval
  if (missing(interval)){
    interval = interval[1]
  }
  else{
    interval = match.arg(interval)
  }
  # if (!is.element( interval, c("none", "confidence"))){
  #   stop("The value of 'interval' is invalid")
  # }

  # Check level
  if (level <= 0 | level >=1){
    stop("The value of 'level' must be larger than 0 and smaller than 1")
  }

  # Set matrices, check number of matrix rows
  if (is.null(X_mu)){
    X_mu = object$X_mu
  }
  else {
    X_mu = as.matrix(X_mu)
    if (object$intercept_mu){
      Intercept = rep.int(1, nrow(X_mu))
      X_mu = cbind( Intercept, X_mu)
    }
  }

  if (is.null(X_sigma)){
    X_sigma = object$X_sigma
  }
  else {
    X_sigma = as.matrix(X_sigma)
    if (object$intercept_sigma){
      Intercept = rep.int(1, nrow(X_sigma))
      X_sigma = cbind( Intercept, X_sigma)
    }
  }

  # Set column names
  names = matrix_column_names( X_mu, X_sigma,
                               intercept_mu = object$intercept_mu, intercept_sigma = object$intercept_sigma)
  colnames(X_mu) = names$colnames_X
  colnames(X_sigma) = names$colnames_X_sigma

  # Check number of rows
  if (nrow(X_mu) != nrow(X_sigma)){
    stop("Number of rows in X_mu and X_sigma are not equal")
  }

  # Check column names
  covnames = names(coef.lmvar( object, sigma = FALSE))
  if (!all(is.element( covnames, colnames(X_mu)))){
    warning("Model matrix for mu does not contain all covariates for mu in 'lmvar' object.
            Beware of unreliable predictions.")
  }
  covnames = names(coef.lmvar( object, mu = FALSE))
  if (!all(is.element( covnames, colnames(X_sigma)))){
    warning("Model matrix for sigma does not contain all covariates for sigma in 'lmvar' object.
            Beware of unreliable predictions.")
  }

  # Remove columns not present in the coefficients
  beta_mu = coef.lmvar( object, sigma = FALSE)
  beta_sigma = coef.lmvar( object, mu = FALSE)

  cols = which(is.element( colnames(X_mu), names(beta_mu)))
  if (length(cols) == 1){
    cnames = colnames(X_mu)[cols]
    X_mu = as.matrix(X_mu[,cols])
    colnames(X_mu) = cnames
  }
  else {
    X_mu = X_mu[,cols]
  }

  cols = which(is.element( colnames(X_sigma), names(beta_sigma)))
  if (length(cols) == 1){
    cnames = colnames(X_sigma)[cols]
    X_sigma = as.matrix(X_sigma[,cols])
    colnames(X_sigma) = cnames
  }
  else {
    X_sigma = X_sigma[,cols]
  }

  # Select beta_mu and beta_sigma elements and re-order
  beta_mu = beta_mu[colnames(X_mu)]
  beta_sigma = beta_sigma[colnames(X_sigma)]

  # Calculate predictions
  predicted_mu = numeric()
  predicted_sigma = numeric()

  if (mu | log | interval == "prediction"){

    predicted_mu = as.numeric(X_mu %*% beta_mu)

    # Store value for later use
    predicted_mu_not_log = predicted_mu
  }

  if (sigma | log | (interval != "none")){

    predicted_sigma = as.numeric(exp(X_sigma %*% beta_sigma))

    # Store value for later use
    predicted_sigma_not_log = predicted_sigma
  }

  if (log){
    predicted_mu = exp( predicted_mu + predicted_sigma_not_log^2 / 2)
    if (sigma){
      predicted_sigma = predicted_mu * sqrt(exp(predicted_sigma_not_log^2) - 1)
    }
  }

  # Calculate confidence intervals
  if (interval != "none"){

    z = stats::qnorm( 0.5 + level/2)

    if (interval == "confidence"){

      if (mu | log){

        # Calculate variances of mu
        S = vcov.lmvar( object, sigma = FALSE)
        S = S[ colnames(X_mu), colnames(X_mu)]
        S = chol(S)
        S = X_mu %*% Matrix::t(S)
        variances_mu = Matrix::rowSums(S*S)
      }

      if (sigma | log){

        # Calculate variances of log sigma
        S = vcov.lmvar( object, mu = FALSE)
        S = S[ colnames(X_sigma), colnames(X_sigma)]
        S = chol(S)
        S = X_sigma %*% Matrix::t(S)
        variances_log_sigma = Matrix::rowSums(S*S)
      }

      # Calculate intervals for mu
      if (mu){
        if (!log){
          delta =  z * sqrt(variances_mu)
        }
        else {
          delta = z * predicted_mu * sqrt(variances_mu + predicted_sigma_not_log^4 * variances_log_sigma)
        }
        mu_down = predicted_mu - delta
        mu_up = predicted_mu + delta
      }

      # Calculate intervals for sigma
      if (sigma){
        if (!log){
          delta = z * predicted_sigma * sqrt(variances_log_sigma)
        }
        else {
          delta = predicted_sigma^2 * variances_mu  +
            (2 * predicted_sigma + predicted_mu^2 / predicted_sigma)^2 * predicted_sigma_not_log^4 * variances_log_sigma
          delta = z * sqrt(delta)
        }
        sigma_down = predicted_sigma - delta
        sigma_up = predicted_sigma + delta
      }
    }
    else if (interval == "prediction"){

      lwr = predicted_mu_not_log - z * predicted_sigma_not_log
      upr = predicted_mu_not_log + z * predicted_sigma_not_log

      if (log){
        lwr = exp(lwr)
        upr = exp(upr)
      }
    }
  }

  # Construct output for expected values
  mu_out = numeric()
  if (mu){
    if (interval %in% c( "none", "prediction")){
      mu_out = matrix( predicted_mu, ncol = 1)
      colnames(mu_out) = "mu"
    }
    else if (interval == "confidence"){
      mu_out = cbind( predicted_mu, mu_down, mu_up)
      colnames(mu_out) = c("mu", "mu_lwr", "mu_upr")
    }
  }

  # Construct output for standard deviations
  sigma_out = numeric()
  if (sigma){
    if (interval %in% c( "none", "prediction")){
      sigma_out = matrix( predicted_sigma, ncol = 1)
      colnames(sigma_out) = "sigma"
    }
    else if (interval == "confidence"){
      sigma_out = cbind( predicted_sigma, sigma_down, sigma_up)
      colnames(sigma_out) = c("sigma", "sigma_lwr", "sigma_upr")
    }
  }

  # Construct full output
  out = cbind( mu_out, sigma_out)
  if (ncol(out) == 1){
    out = as.numeric(out)
  }
  else{
    colnames(out) = c( colnames(mu_out), colnames(sigma_out))
  }

  if (interval == "prediction"){
    out = cbind( out, lwr, upr)
  }

  return(out)
}
