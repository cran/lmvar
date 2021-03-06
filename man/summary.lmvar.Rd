% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.lmvar.R
\name{summary.lmvar}
\alias{summary.lmvar}
\title{Summary overview for an object of class 'lmvar'}
\usage{
\method{summary}{lmvar}(object, mu = TRUE, sigma = TRUE, ...)
}
\arguments{
\item{object}{Object of class 'lmvar'}

\item{mu}{Boolean, specifies whether or not to include the coefficients \eqn{\beta_\mu} in the table of coefficients}

\item{sigma}{Boolean, specifies whether or not to include the coefficients \eqn{\beta_\sigma} in the table of coefficients}

\item{...}{For compatibility with \code{\link[base]{summary}} generic}
}
\value{
An object of class 'summary_lmvar'. This is a list with the following members:
\itemize{
\item \code{call} Call that created \code{object}

\item \code{coefficients} Data frame
with one row for each element of \eqn{\beta_\mu} and \eqn{\beta_\sigma} and the following variables.
\itemize{
\item \code{Estimate} maximum-likelihood estimate
\item \code{Std. Error} standard error, defined as \eqn{\sqrt(var(\beta))} with \eqn{var(\beta)} the estimated variance
of \eqn{\beta}.
\item \code{z value} z-statistic, defined as \eqn{\beta / \sqrt(var(\beta))}
\item \code{Pr(>|z|)} p-value of the z-statistic, calculated from the standard normal distribution.
}
\item \code{residuals} A numeric vector with the minimum, the 25\% quartile, the median, the 75\% quartile and the maximum
standardized residual. The standardized residual of an observation is defined as \eqn{(y - \mu) / \sigma} where \eqn{y} is the value
of the observation, \eqn{\mu} the expectation value and \eqn{\sigma}
the standard deviation of the observation.
\item \code{sigma} A numeric vector with the minimum, the 25\% quartile, the median, the 75\% quartile and the maximum
standard deviation \eqn{\sigma} of all observations.
\item \code{aliased_mu} A named logical vector. The names are the column names of the user-supplied model matrix
\eqn{X_\mu}. The values (\code{TRUE} or \code{FALSE}) tell whether or not the column
has been removed by \code{lmvar} to make the matrix full-rank.
\item \code{aliased_sigma} As \code{aliased_mu} but for the user-supplied model matrix \eqn{X_\sigma}.
\item \code{logLik_ratio} The difference in log-likelihood between the model in \code{object} and a classical linear
model with model matrix \eqn{X_\mu} and a constant variance for all observations.
\item \code{df_additional} The difference in degrees in freedom between the model in \code{object} and a classical linear
model with model matrix \eqn{X_\mu} and a constant variance for all observations. Is equal to \code{NULL} if
\eqn{X_\sigma} does not contain an intercept term.
\item \code{p_value} The p-value of \code{2 loglik_ratio}, calculated from a chi-squared distribution with \code{df}
degrees of freedom. Is equal to \code{NULL} if there are no additional degrees of freedom.
\item \code{nobs} The number of observations in \code{object}.
\item \code{df} The degrees of freedom of the fit in \code{object}.
\item \code{options} A list of argument-values of the function call.
}
}
\description{
Summary overview for an object of class 'lmvar'.
}
\details{
Standard errors and z-statistics are calculated under the assumption of asymptotic normality for maximum
likelihood estimators. They may not
be reliable when the number of observations in \code{object} is small.
}
\examples{
# As example we use the dataset 'attenu' from the library 'datasets'. The dataset contains
# the response variable 'accel' and two explanatory variables 'mag'  and 'dist'.
library(datasets)

# Create the model matrix for the expected values
X = cbind(attenu$mag, attenu$dist)
colnames(X) = c("mag", "dist")

# Create the model matrix for the standard deviations.
X_s = cbind(attenu$mag, 1 / attenu$dist)
colnames(X_s) = c("mag", "dist_inv")

# Carry out the fit
fit = lmvar(attenu$accel, X, X_s)

# Print a summary of the fit
summary(fit)

# Include only the coefficients beta for the expected values
summary(fit, sigma = FALSE)

# Include only the coefficients beta for the standard deviations
summary(fit, mu = FALSE)

# Extract the matrix of coefficients from the summary
coef(summary(fit))
}
\seealso{
\code{\link[stats]{coef}} to extract the matrix with estimates, standard-errors, t-statistics and
p-values for \eqn{\beta_\mu} and \eqn{\beta_\sigma} from a 'summary_lmvar' object.

\code{\link{vcov.lmvar}} for the covariance matrix of the \eqn{\beta_\mu} and \eqn{\beta_\sigma} in an object of class
'lmvar'.

\code{\link{print.summary_lmvar}} for a print method for a 'summary_lmvar' object.

\code{\link{fitted.lmvar}} for the expected values and standard deviations
of the observations in an object of class 'lmvar'.

\code{\link{logLik.lmvar}} for the log-likelihood of a fit in an object of class 'lmvar'.

\code{\link{alias.lmvar_no_fit}} to obtain the aliased columns of the user-supplied model matrices in the call of \code{\link{lmvar}}.
}
