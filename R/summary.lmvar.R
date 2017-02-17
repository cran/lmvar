#' @title Summary overview for an object of class 'lmvar'
#'
#' @description Summary overview for an object of class 'lmvar'.
#'
#' @param object Object of class 'lmvar'
#' @param ... For compatibility with \code{\link[base]{summary}} generic
#'
#' @return An object of class 'summary_lmvar'. This is a list with the following members:
#' \itemize{
#' \item \code{call} Call that created \code{object}
#'
#' \item \code{coefficients} Data frame
#' with one row for each element of \eqn{\beta_\mu} and \eqn{\beta_\sigma} and the following variables.
#' \itemize{
#' \item \code{Estimate} maximum-likelihood estimate
#' \item \code{Std. Error} standard error, defined as \eqn{\sqrt(var(\beta))} with \eqn{var(\beta)} the estimated variance
#' of \eqn{\beta}.
#' \item \code{t value} t-statistic, defined as \eqn{\beta / \sqrt(var(\beta))}
#' \item \code{Pr(>|t|)} p-value of the t-statistic, calculated from the normal distribution.
#' }
#' \item \code{residuals} A numeric vector with the minimum, the 25\% quartile, the median, the 75\% quartile and the maximum
#' standardized residual. The standardized residual of an observation is defined as \eqn{(y - \mu) / \sigma} where \eqn{y} is the value
#' of the observation, \eqn{\mu} the expectation value and \eqn{\sigma}
#' the standard deviation of the observation.
#' \item \code{sigma} A numeric vector with the minimum, the 25\% quartile, the median, the 75\% quartile and the maximum
#' standard deviation \eqn{\sigma} of all observations.
#' \item \code{aliased_mu} A named logical vector. The names are the column names of the user-supplied model matrix
#' \eqn{X_\mu}. The values (\code{TRUE} or \code{FALSE}) tell whether or not the column
#' has removed by \code{lmvar} to make the matrix full-rank.
#' \item \code{aliased_sigma} As \code{aliased_mu} but for the user-supplied model matrix \eqn{X_\sigma}.
#' \item \code{logLik_ratio} The difference in log-likelihood between the model in \code{object} and a classical linear
#' model with model matrix \eqn{X_\mu} and a constant variance for all observations.
#' \item \code{df} The difference in degrees in freedom between the model in \code{object} and a classical linear
#' model with model matrix \eqn{X_\mu} and a constant variance for all observations.
#' \item \code{p_value} The p-value of \code{2 loglik_ratio}, calculated from a chi-squared distribution with \code{df}
#' degrees of freedom.
#' }
#'
#' @details Standard errors and t-statistics are calculated from the estimated covariance matrix and may not
#' be reliable when the number of observations in \code{object} is small.
#'
#' @seealso \code{\link[stats]{coef}} to extract the matrix with estimates, standard-errors, t-statistics and
#' p-values for \eqn{\beta_\mu} and \eqn{\beta_\sigma} from a 'summary_lmvar' object.
#'
#' \code{\link{vcov.lmvar}} for the covariance matrix of the \eqn{\beta_\mu} and \eqn{\beta_\sigma} in an object of class
#' 'lmvar'.
#'
#' \code{\link{print.summary_lmvar}} for a print method for a 'summary_lmvar' object.
#'
#' \code{\link{fitted.lmvar}} for the expected values and standard deviations
#' of the observations in an object of class 'lmvar'.
#'
#' \code{\link{logLik.lmvar}} for the log-likelihood of a fit in an object of class 'lmvar'.
#'
#' \code{\link{alias.lmvar}} to obtain the aliased columns of the user-supplied model matrices in the call of \code{\link{lmvar}}.
#'
#' @method summary lmvar
#'
#' @export
#'
#' @example R/examples/summary_examples.R
#'
summary.lmvar <- function(object, ...){

  beta_mu = coef.lmvar( object, sigma = FALSE)
  beta_sigma = coef.lmvar( object, mu = FALSE)

  # get the names of the betas
  beta_sigma_names = beta_sigma_names( names(beta_mu), names(beta_sigma))
  beta_names = c( names(beta_mu), beta_sigma_names)

  # Calculate varianve-covariance matrices
  I_mu = vcov.lmvar( object, sigma = FALSE)
  I_sigma = vcov.lmvar( object, mu = FALSE)

  # Calculate variance-related statistics
  variances = c( Matrix::diag(I_mu), Matrix::diag(I_sigma))
  sterr = sqrt(variances)
  z_values = c(beta_mu, beta_sigma) / sterr
  pr_z = 2 * stats::pnorm( abs(z_values), lower.tail = FALSE)

  # calculate data-vectors for dataframe 'coefficients'
  estimate = c( beta_mu, beta_sigma)

  coeff = data.frame( "Estimate"= estimate, "Std. Error"=sterr, "z value"=z_values, "Pr(>|z|)"=pr_z,
                      row.names = beta_names, check.names = FALSE)

  # Retrieve sigmas
  sigma = fitted.lmvar( object, mu=FALSE)

  # Calculate standardized residuals
  res = residuals.lmvar(object) / sigma
  res = stats::quantile( res, c( 0, 0.25, 0.5, 0.75, 1))
  names(res) = c( "Min", "1Q", "Median", "3Q", "Max")

  # Summarize sigmas in quantiles
  sigma = stats::quantile( sigma, c( 0, 0.25, 0.5, 0.75, 1))
  names(sigma) = c( "Min", "1Q", "Median", "3Q", "Max")


  rlist = list( call = object$call,
                residuals = res,
                coefficients = coeff,
                sigma = sigma,
                aliased_mu = object$aliased_mu,
                aliased_sigma = object$aliased_sigma,
                logLik_ratio = object$logLik - object$logLik_lm,
                df = dfree( object, mu = FALSE) - 1,
                p_value = stats::pchisq(2 * (object$logLik - object$logLik_lm), dfree( object, mu = FALSE) - 1, lower.tail = FALSE))

  class(rlist) = "summary_lmvar"

  return(rlist)
}
