#' @title Fitted values for an 'lmvar' object
#'
#' @description Estimators and confidence intervals for the expected values and standard deviations of the response vector \eqn{Y}.
#' Prediction intervals for \eqn{Y}. Alternatively, estimators and intervals can be for \eqn{e^Y}.
#'
#' @param object An 'lmvar' object
#' @param mu Boolean, specifies whether or not to return estimators and intervals for the expected values
#' @param sigma Boolean, specifies whether or not to return estimators and intervals for the standard deviations
#' @param log Boolean, specifies whether estimators and intervals should be for \eqn{Y} (\code{log = FALSE}) or
#' for \eqn{e^Y} (\code{log = TRUE}).
#' @param interval Character string, specifying the type of interval. Possible values are
#' \itemize{
#' \item "none" No interval, this is the default
#' \item "confidence" Confidence intervals for the estimators
#' \item "prediction" Prediction intervals
#' }
#' @param level Numeric value between 0 and 1, specifying the confidence level
#' @param ... For compatibility with \code{\link[stats]{fitted}} generic.
#'
#' @return In the case \code{mu = FALSE} and \code{interval = "none"}: a numeric vector containing the estimators for
#' the standard deviation.
#'
#' In the case \code{sigma = FALSE} and \code{interval = "none"}: a numeric vector containing the estimators for
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
#' @details If \code{log = FALSE}, \code{fitted.lmvar} returns
#' estimators and intervals for the observations \eqn{Y} stored in \code{object}.
#'
#' If \code{log = TRUE}, \code{fitted.lmvar} returns estimators and intervals for \eqn{e^Y}.
#'
#' Confidence intervals are calculated under teh assumption of asymptotic normality. This assumption holds when the number of observations is
#' large. Intervals must be treated cautiously in case of a small number of observations.
#' Intervals can also be unreliable if
#' \code{object} was created with a constraint on the minimum values of the standard deviations sigma.
#'
#' This function is identical to the function \code{\link{predict.lmvar}} in which the parameters \code{X_mu} and
#' \code{X_sigma} are left unspecified.
#'
#' @seealso \code{\link{predict.lmvar}} for expected values, standard deviations and intervals for model matrices different from
#' the ones present in \code{object}.
#'
#' \code{\link{coef.lmvar}} and \code{\link[stats]{confint}} for maximum likelihood estimators and confidence intervals for
#' \eqn{\beta_\mu} and \eqn{\beta_\sigma}.
#'
#' @example R/examples/fitted_examples.R
#'

fitted.lmvar <- function( object, mu = TRUE, sigma = TRUE, log = FALSE, interval = c("none", "confidence", "prediction"), level = 0.95, ...){

  # Check interval
  if (missing(interval)){
    interval = interval[1]
  }
  else{
    interval = match.arg(interval)
  }

  return(predict.lmvar( object, mu = mu, sigma = sigma, log = log, interval = interval, level = level, ...))
}
