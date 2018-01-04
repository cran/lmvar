#' @title Summary overview for an object of class 'lmvar'
#'
#' @description Summary overview for an object of class 'lmvar'.
#'
#' @param object Object of class 'lmvar'
#' @param mu Boolean, specifies whether or not to include the coefficients \eqn{\beta_\mu} in the table of coefficients
#' @param sigma Boolean, specifies whether or not to include the coefficients \eqn{\beta_\sigma} in the table of coefficients
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
#' \item \code{z value} z-statistic, defined as \eqn{\beta / \sqrt(var(\beta))}
#' \item \code{Pr(>|z|)} p-value of the z-statistic, calculated from the standard normal distribution.
#' }
#' \item \code{residuals} A numeric vector with the minimum, the 25\% quartile, the median, the 75\% quartile and the maximum
#' standardized residual. The standardized residual of an observation is defined as \eqn{(y - \mu) / \sigma} where \eqn{y} is the value
#' of the observation, \eqn{\mu} the expectation value and \eqn{\sigma}
#' the standard deviation of the observation.
#' \item \code{sigma} A numeric vector with the minimum, the 25\% quartile, the median, the 75\% quartile and the maximum
#' standard deviation \eqn{\sigma} of all observations.
#' \item \code{aliased_mu} A named logical vector. The names are the column names of the user-supplied model matrix
#' \eqn{X_\mu}. The values (\code{TRUE} or \code{FALSE}) tell whether or not the column
#' has been removed by \code{lmvar} to make the matrix full-rank.
#' \item \code{aliased_sigma} As \code{aliased_mu} but for the user-supplied model matrix \eqn{X_\sigma}.
#' \item \code{logLik_ratio} The difference in log-likelihood between the model in \code{object} and a classical linear
#' model with model matrix \eqn{X_\mu} and a constant variance for all observations.
#' \item \code{df_additional} The difference in degrees in freedom between the model in \code{object} and a classical linear
#' model with model matrix \eqn{X_\mu} and a constant variance for all observations. Is equal to \code{NULL} if
#' \eqn{X_\sigma} does not contain an intercept term.
#' \item \code{p_value} The p-value of \code{2 loglik_ratio}, calculated from a chi-squared distribution with \code{df}
#' degrees of freedom. Is equal to \code{NULL} if there are no additional degrees of freedom.
#' \item \code{nobs} The number of observations in \code{object}.
#' \item \code{df} The degrees of freedom of the fit in \code{object}.
#' \item \code{options} A list of argument-values of the function call.
#' }
#'
#' @details Standard errors and z-statistics are calculated under the assumption of asymptotic normality for maximum
#' likelihood estimators. They may not
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
#' \code{\link{alias.lmvar_no_fit}} to obtain the aliased columns of the user-supplied model matrices in the call of \code{\link{lmvar}}.
#'
#' @method summary lmvar
#'
#' @export
#'
#' @example R/examples/summary_examples.R
#'
summary.lmvar <- function(object, mu = TRUE, sigma = TRUE, ...){

  # Calculate betas
  beta_mu = numeric()
  beta_sigma = numeric()
  if (mu){
    beta_mu = coef.lmvar( object, sigma = FALSE)
  }
  if (sigma){
    beta_sigma = coef.lmvar( object, mu = FALSE)
  }

  # get the names of the betas
  beta_sigma_names = beta_sigma_names( names(beta_mu), names(beta_sigma))
  beta_names = c( names(beta_mu), beta_sigma_names)

  # Calculate variances
  variances_mu = numeric()
  variances_sigma = numeric()
  if (mu){
    M = vcov.lmvar( object, sigma = FALSE)
    variances_mu = Matrix::diag(M)
  }
  if (sigma){
    M = vcov.lmvar( object, mu = FALSE)
    variances_sigma = Matrix::diag(M)
  }

  # Calculate variance-related statistics
  variances = c( variances_mu, variances_sigma)
  sterr = sqrt(variances)
  z_values = c(beta_mu, beta_sigma) / sterr
  pr_z = 2 * stats::pnorm( abs(z_values), lower.tail = FALSE)

  # calculate data-vectors for dataframe 'coefficients'
  estimate = c( beta_mu, beta_sigma)

  coeff = cbind( estimate, sterr, z_values, pr_z)
  colnames(coeff) = c( "Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(coeff) = beta_names

  # Retrieve sigmas
  sigma_fit = fitted.lmvar( object, mu=FALSE)

  # Calculate standardized residuals
  res = residuals.lmvar(object) / sigma_fit
  res = stats::quantile( res, c( 0, 0.25, 0.5, 0.75, 1))
  names(res) = c( "Min", "1Q", "Median", "3Q", "Max")

  # Summarize sigmas in quantiles
  sigma_fit = stats::quantile( sigma_fit, c( 0, 0.25, 0.5, 0.75, 1))
  names(sigma_fit) = c( "Min", "1Q", "Median", "3Q", "Max")

  # Calculate additional degrees of freedom
  if (object$intercept_sigma){
    df_additional = dfree( object, mu = FALSE) - 1
  }
  else {
    df_additional = NULL
  }

  # Calculate p-value for difference in deviance
  if (is.null(df_additional)){
    p_value = NULL
  }
  else {
    p_value = stats::pchisq(2 * (object$logLik - object$logLik_lm), df_additional, lower.tail = FALSE)
  }

  rlist = list( call = object$call,
                residuals = res,
                coefficients = coeff,
                sigma = sigma_fit,
                aliased_mu = object$aliased_mu,
                aliased_sigma = object$aliased_sigma,
                logLik_ratio = object$logLik - object$logLik_lm,
                df_additional = df_additional,
                p_value = p_value,
                nobs = nobs(object),
                df = dfree(object),
                options = list(mu = mu, sigma = sigma))

  class(rlist) = "summary_lmvar"

  return(rlist)
}
