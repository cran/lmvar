#' @title Linear regression with non-constant variances
#'
#' @description Performs a Gaussian regression with non-constant variances and outputs an 'lmvar' object.
#'
#' @param y Vector of observations
#' @param X_mu Model matrix for the expected values \eqn{\mu}
#' @param X_sigma Model matrix for the logarithms of the standard deviations \eqn{\sigma}
#' @param check_hessian Boolean, if TRUE it is checked whether the Hessian is negative-definite
#' @param ... Additional arguments, not used in the current implementation
#'
#' @return An object of class 'lmvar'.
#'
#' @details If the matrix \code{X_mu} is not specified, the model for the expected values \eqn{\mu} will consist of an intercept term only.
#' Likewise, the model for \eqn{\log \sigma} will consist of an intercept term only if \code{X_sigma} is not specified.
#' In the latter case, the model reduces to a standard linear model.
#'
#' The input matrices can be of class \code{\link{matrix}} or of class \link[Matrix:Matrix-class]{Matrix}
#'
#' When \code{check_hessian = TRUE}, it is checked whether the fitted log-likelihood is at a maximum. A warning will be issued if
#' that is not the case.
#'
#' The vector \code{y} must be a numeric vector. It can not contain special values like \code{NULL}, \code{NaN}, etc.
#'
#' Both model matrices must be numeric matrices. They can not contain special values like \code{NULL}, \code{NaN}, etc.
#'
#' All columns in the matrix \code{X_mu} must either have a name, or no column has a name at all. It is not allowed that some
#' colums have a
#' name while others don't. The same is true for \code{X_sigma}.
#'
#' If supplied, the column names
#' must be unique for \code{X_mu}. The same is true for \code{X_sigma}. A column name can be used in both \code{X_mu} and
#' \code{X_sigma} though.
#'
#' In case the matrix \code{X_mu} has no columns with a column name, \code{lmvar}  gives the names \code{v1},
#' \code{v2} etc. to the columns. Likewise, if the matrix \code{X_sigma} has no columns with a column name, \code{lmvar}
#' gives the names \code{v1_s}, \code{v2_s} etc. to the columns.
#'
#' Matrix \code{X_mu} can not have a column with the name '(Intercept)'.
#' Matrix \code{X_sigma} can not have a column with the name '(Intercept_s)'. Both names are reserved names.
#'
#' The function \code{nleqslv} from the \code{nleqslv} package, is used to solve the non-linear equations that comprise the fit.
#' \code{nleqslv} returns a
#' termination code which reports whether a solution was found, see \code{\link[nleqslv]{nleqslv}}.
#' A potential problem is reported by a termination code different from \code{1}. Termination codes different from \code{1}
#' as passed on by \code{lmvar}. Termination codes \code{2} and \code{3} are passed on as a warning.
#' Other termination codes are passed on as an error.
#'
#' An \code{lmvar} object is a list. Users are discouraged to access list-members directly.
#' Instead, list-members are to be accessed with the various accessor and utility functions in the package.
#' Exceptions are the following list members for which no accessor functions exist:
#' \itemize{
#' \item \code{object$y} the vector of observations
#' \item \code{object$X_mu} the  model matrix for \eqn{\mu}. In general, it differs from the user-supplied \code{X_mu} because
#' \code{lmvar} adds an intercept-column and makes the matrix full-rank.
#' \item \code{object$X_sigma} the  model matrix for \eqn{\log \sigma}. In general, it differs from the user-supplied
#' \code{X_sigma} because
#' \code{lmvar} adds an intercept-column and makes the matrix full-rank.
#' }
#'
#' See the vignettes
#' that come with the \code{lmvar} package for more info. Run \code{vignette(package="lmvar")} to list the available vignettes.
#'
#' @export
#'
#' @example R/examples/lmvar_examples.R
#'
lmvar <- function( y, X_mu = NULL, X_sigma = NULL, check_hessian = TRUE, ...){

  dlogL <- function(beta_sigma){

    # Functions returns derivative (gradient) of log likelihood wrt beta_sigma

    sigma = as.numeric(exp( X_sigma %*% beta_sigma))

    M = X_mu * (1/sigma)
    M = tryCatch( chol2inv(chol(Matrix::t(M) %*% M)),          # Exploit that t(M) %*% M is symmetric and positive-definite

                  error = function(e){
                    i = which.min(beta_sigma)
                    name = colnames(X_sigma)[i]
                    value = beta_sigma[i]
                    message( "Smallest element of beta_sigma is '", name, "' with value ", signif( value, 3))
                    i = which.max(beta_sigma)
                    name = colnames(X_sigma)[i]
                    value = beta_sigma[i]
                    message( "Largest element of beta_sigma is '", name, "' with value ", signif( value, 3))
                    stop(e)
                  })

    beta_mu = M %*% Matrix::t(X_mu) %*% (y / (sigma * sigma))
    mu = as.numeric(X_mu %*% beta_mu)

    res = (y - mu) / sigma

    dBeta_sigma = as.numeric(Matrix::t(X_sigma) %*% (res*res - 1))  # gradient of log-likelihood w.r.t. beta_sigma

    return(dBeta_sigma)
  }

  hesslogL <- function(beta_sigma){

    # Functions returns hessian matrix of log likelihood

    sigma = as.numeric(exp( X_sigma %*% beta_sigma))

    M = X_mu * (1/sigma)
    M = chol2inv(chol(Matrix::t(M) %*% M))          # Exploit that t(M) %*% M is symmetric and positive-definite

    beta_mu = M %*% Matrix::t(X_mu) %*% (y / (sigma * sigma))
    mu = as.numeric(X_mu %*% beta_mu)
    res = (y - mu) / sigma

    XT = chol(M) %*% Matrix::t(X_mu) %*% (X_sigma * (res / sigma)) # Use Choleski decomposition to have H_ss_1 symmetric
    H_ss_1 = 4 * as.matrix(Matrix::t(XT) %*% XT)

    XT = X_sigma * abs(res)
    H_ss_2 = -2 * as.matrix(Matrix::t(XT) %*% XT)

    return(H_ss_1 + H_ss_2)
  }

  # Check inputs
  if (!is.null(X_mu)){
    if (nrow(as.matrix(X_mu)) != length(y)){
      stop("Number of rows of matrix X_mu not equal to length of response vector y")
    }
    if (is.element( "(Intercept)", colnames(X_mu))){
      stop("Matrix X_mu can not have a column with the name '(Intercept)'")
    }
  }
  if
  (!is.null(X_sigma)){
    if (nrow(as.matrix(X_sigma)) != length(y)){
      stop("Number of rows of matrix X_sigma not equal to length of response vector y")
    }
    if (is.element( "(Intercept_s)", colnames(X_sigma))){
      stop("Matrix X_sigma can not have a column with the name '(Intercept_s)'")
    }
  }

  n = length(y)

  # Add intercept column to X_mu and X_sigma
  Intercept = rep.int(1, n)
  if (is.null(X_mu)){
    X_mu = matrix(Intercept, nrow=n, ncol=1)
  }
  else {
    X_mu = cbind( Intercept, X_mu)
  }

  if (is.null(X_sigma)){
    X_sigma = matrix(Intercept, nrow=n, ncol=1)
  }
  else {
    X_sigma = cbind( Intercept, X_sigma)
  }

  # Initialize aliased columns
  aliased_mu = logical(length = ncol(X_mu))
  names(aliased_mu) = matrix_column_names( X_mu, X_sigma)$colnames_X
  aliased_sigma = logical(length = ncol(X_sigma))
  names(aliased_sigma) = matrix_column_names( X_mu, X_sigma)$colnames_X_sigma

  # Turn matrices into a full-rank matrices
  X_mu = make_matrix_full_rank(X_mu)
  if (class(X_mu) == 'numeric'){
    X_mu = as.matrix(X_mu)
  }
  X_sigma = make_matrix_full_rank(X_sigma)
  if (class(X_sigma) == 'numeric'){
    X_sigma = as.matrix(X_sigma)
  }

  # Give matrices column names if there are none
  colnames(X_mu) = matrix_column_names( X_mu, X_sigma)$colnames_X
  colnames(X_sigma) = matrix_column_names( X_mu, X_sigma)$colnames_X_sigma

  # set aliased columns
  a_names = names(aliased_mu)
  aliased_mu = !is.element( a_names, colnames(X_mu))
  names(aliased_mu) = a_names
  a_names = names(aliased_sigma)
  aliased_sigma = !is.element( a_names, colnames(X_sigma))
  names(aliased_sigma) = a_names

  if (nrow(X_mu) < (ncol(X_mu) + ncol(X_sigma))){
    stop("Too few observations. There must be at least ", ncol(X_mu) + ncol(X_sigma), " observations.")
  }

  # Initial estimate by ordinary linear regression
  fit = stats::lm( y ~ .+0, as.data.frame(as.matrix(X_mu)), model = FALSE)
  beta = stats::coef(fit)
  sigma = summary(fit)$sigma
  beta_sigma = log(sigma) * solve( Matrix::t(X_sigma) %*% X_sigma) %*% Matrix::colSums(X_sigma)  # Minimize sum over all elements of (X_sigma %*% beta_sigma - log(sigma))^2
  beta_sigma = as.numeric(beta_sigma)

  l = length(beta)

  # Likelihood of model with constant variance
  sigma2 = fit$residuals %*% fit$residuals / n
  logLik_lm = -n * (log(2*pi*sigma2) + 1) / 2
  logLik_lm = as.numeric(logLik_lm)

  # Calculate betas
  solve_result = nleqslv::nleqslv( beta_sigma, dlogL, jac = hesslogL)

  if (solve_result$termcd == 2 | solve_result$termcd == 3){
    max_value = max(abs(dlogL(solve_result$x)))
    warning( c( solve_result$message, ". Maximum absolute value in gradient should be within 1e-08 but is ", signif( max_value, 4)))
  }
  else if (solve_result$termcd > 3 | solve_result$termcd < 1){
    stop( solve_result$message)
  }
  beta_sigma = solve_result$x

  sigma = as.numeric(exp( X_sigma %*% beta_sigma))
  M = X_mu * (1/sigma)
  M = solve(Matrix::t(M) %*% M)
  beta = as.numeric(M %*% Matrix::t(X_mu) %*% (y / (sigma * sigma)))


  # Check Hessian
  if (check_hessian){
    if(!matrixcalc::is.negative.definite(hesslogL(beta_sigma))){
      warning("Log-likelihood appears not to be at a maximum!")
    }
  }

  # Calculate log-likelihood
  mu = as.numeric(X_mu %*% beta)
  sigma = as.numeric(exp(X_sigma %*% beta_sigma))
  res = (y-mu) / sigma
  logLik= -0.5*n*log(2*pi) - sum(log(sigma)) - 0.5*sum(res*res)

  # Associate betas and other result vectors with names of covariates
  names(beta) = colnames(X_mu)
  names(beta_sigma) = colnames(X_sigma)

  rlist = list( call = match.call(),
                coefficients_mu = beta,
                coefficients_sigma = beta_sigma,
                logLik = logLik,
                logLik_lm = logLik_lm,
                aliased_mu = aliased_mu,
                aliased_sigma = aliased_sigma,
                y = y,
                X_mu = X_mu,
                X_sigma = X_sigma)

  class(rlist) = "lmvar"

  return(rlist)

}
