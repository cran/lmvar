#' @title Linear regression with non-constant variances
#'
#' @description Performs a Gaussian regression with non-constant variances and outputs an 'lmvar' object.
#'
#' @param y Vector of observations
#' @param X_mu Model matrix for the expected values \eqn{\mu}
#' @param X_sigma Model matrix for the logarithms of the standard deviations \eqn{\sigma}
#' @param check_hessian Boolean, if TRUE it is checked whether the Hessian is negative-definite
#' @param slvr_options A list with options to be passed on to the function
#' \code{maxNR} from the package \code{maxLik} which maximizes the log-liklihood
#' @param slvr_log Boolean, if TRUE the output of \code{maxNR} is added to the 'lmvar' object
#' @param intercept_mu Boolean, if TRUE a column with the name '(Intercept)' is added to the matrix \code{X_mu}. This column
#' represents the intercept term in the model for \eqn{\mu}.
#' @param intercept_sigma Boolean, if TRUE a column with the name '(Intercept_s)' is added to the matrix \code{X_sigma}. This column
#' represents the intercept term in the model for \eqn{\log \sigma}.
#' @param ... Additional arguments, not used in the current implementation
#'
#' @return An object of class 'lmvar'.
#'
#' @details If the matrix \code{X_mu} is not specified, the model for the expected values \eqn{\mu} will consist of an intercept term only.
#' Likewise, the model for \eqn{\log \sigma} will consist of an intercept term only if \code{X_sigma} is not specified.
#' In the latter case, the model reduces to a standard linear model.
#'
#' The input matrices can be of class \code{\link{matrix}} or of class \link[Matrix:Matrix-class]{Matrix}. They can also be of
#' class \code{\link{numeric}} in case they are matrices with one column only.
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
#' The function \code{maxNR} from the \code{maxLik} package, is used to maximize the (profile) log-likelihood.
#' \code{maxNR} returns a
#' termination code which reports whether a maximum was found, see \code{\link[maxLik]{maxNR}}.
#' A potential problem is reported by a termination code different from \code{1}, \code{2} or \code{8}. In case of a different
#' termination code, the message from \code{maxNR} is output as a warning.
#'
#' All list elements in \code{slvr_options} are passed on as arguments to \code{maxNR}. The name of the list element is
#' the argument name, the value of the list element is the argument value.
#' It is not allowed to pass on the arguments \code{fn}, \code{grad}
#' or \code{hess}.
#'
#' In case one wants to supply an initial estimate for the coefficients, one has to supply an initial estimate
#' for \eqn{\beta_\sigma}. If \code{beta_sigma_initial} is the initial estimate, one pass on the argument
#' \code{slvr_options = list(start = beta_sigma_initial)}.
#' The inital estimate \code{beta_sigma_initial} must be a numeric vector with a length equal to the number
#' of columns of \code{X_sigma} plus one. The first element of the vector is the initial estimate of the intercept term for
#' \eqn{sigma}, the second element is the inital estimate corresponding to the first column in \code{X_sigma}, etc.
#'
#' It is not possible (and there is no need) to supply an inital estimate for \eqn{\beta_\mu}.
#'
#' A \code{lmvar} object is a list. Users are discouraged to access list-members directly.
#' Instead, list-members are to be accessed with the various accessor and utility functions in the package.
#' Exceptions are the following list members for which no accessor functions exist:
#' \itemize{
#' \item \code{object$y} the vector of observations
#' \item \code{object$X_mu} the  model matrix for \eqn{\mu}. In general, it differs from the user-supplied \code{X_mu} because
#' \code{lmvar} adds an intercept-column and makes the matrix full-rank.
#' \item \code{object$X_sigma} the  model matrix for \eqn{\log \sigma}. In general, it differs from the user-supplied
#' \code{X_sigma} because
#' \code{lmvar} adds an intercept-column and makes the matrix full-rank.
#' \item \code{slvr_log} the output of \code{maxNR} (the solver routine used to maximize the likelihood). Included only
#' if the argument \code{slvr_log} has the value \code{TRUE}.
#' See \code{\link[maxLik]{maxNR}} for details about this output.
#' }
#'
#' See the vignettes
#' that come with the \code{lmvar} package for more info. Run \code{vignette(package="lmvar")} to list the available vignettes.
#'
#' @export
#'
#' @example R/examples/lmvar_examples.R
#'
lmvar <- function( y, X_mu = NULL, X_sigma = NULL, check_hessian = TRUE, slvr_options = list(), slvr_log = FALSE,
                   intercept_mu = TRUE, intercept_sigma = TRUE, ...){

  logLHood <- function(beta_sigma){

    sigma_inv = as.numeric(exp(- X_sigma %*% beta_sigma))
    M = X_mu * sigma_inv
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

                    i = which.min(sigma_inv)
                    value = sigma_inv[i]
                    message( "Smallest element of sigma inverse is ", signif( value, 3), " for observation ", i)

                    i = which.max(sigma_inv)
                    value = sigma_inv[i]
                    message( "Largest element of sigma inverse is ", signif( value, 3), " for observation ", i)
                    terms = X_sigma[i,] * beta_sigma
                    i = which.max(terms)
                    message("   largest contribution to predictor is ", signif( terms[i], 3), " for covariate ", colnames(X_sigma)[i],
                            " with beta is ", signif( beta_sigma[i], 3))
                    stop(e)
                  })

    beta_mu = M %*% Matrix::t(X_mu) %*% (y * sigma_inv * sigma_inv)
    mu = as.numeric(X_mu %*% beta_mu)
    res = (y - mu) * sigma_inv

    # Calculate log-likelihood
    logLik= -0.5*n*log(2*pi) + sum(log(sigma_inv)) - 0.5*sum(res*res)

    # Calculate gradient
    dBeta_sigma = as.numeric(Matrix::t(X_sigma) %*% (res*res - 1))  # gradient of log-likelihood w.r.t. beta_sigma

    # Calculate Hessian
    XT = chol(M) %*% Matrix::t(X_mu) %*% (X_sigma * (res * sigma_inv)) # Use Choleski decomposition to have H_ss_1 symmetric
    H_ss_1 = 4 * Matrix::t(XT) %*% XT
    XT = X_sigma * abs(res)
    H_ss_2 = -2 * Matrix::t(XT) %*% XT

    # Set attributes
    attr( logLik, "gradient") = dBeta_sigma
    attr( logLik, "hessian") = H_ss_1 + H_ss_2

    return(logLik)
  }

  # Set whether or not matrices have been specified
  isnull_X = FALSE
  isnull_X_sigma = FALSE
  if(is.null(X_mu)){
    isnull_X = TRUE
    intercept_mu = TRUE
  }
  if(is.null(X_sigma)){
    isnull_X_sigma = TRUE
    intercept_sigma = TRUE
  }

  # Turn vectors into matrices
  if (class(X_mu) == 'numeric'){
    X_mu = as.matrix(X_mu)
  }
  if (class(X_sigma) == 'numeric'){
    X_sigma = as.matrix(X_sigma)
  }

  # Check inputs
  if (!isnull_X){
    if (nrow(as.matrix(X_mu)) != length(y)){
      stop("Number of rows of matrix X_mu not equal to length of response vector y")
    }
    if (is.element( "(Intercept)", colnames(X_mu))){
      stop("Matrix X_mu can not have a column with the name '(Intercept)'")
    }
  }
  if
  (!isnull_X_sigma){
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
  if (isnull_X){
    X_mu = matrix(Intercept, nrow=n, ncol=1)
  }
  else if (intercept_mu){
    X_mu = cbind( Intercept, X_mu)
  }

  if (isnull_X_sigma){
    X_sigma = matrix(Intercept, nrow=n, ncol=1)
  }
  else if (intercept_sigma){
    X_sigma = cbind( Intercept, X_sigma)
  }

  # Initialize aliased columns
  aliased_mu = logical(length = ncol(X_mu))
  aliased_sigma = logical(length = ncol(X_sigma))
  names = matrix_column_names( X_mu, X_sigma, intercept_mu, intercept_sigma)
  names(aliased_mu) = names$colnames_X
  names(aliased_sigma) = names$colnames_X_sigma

  # Turn matrices into a full-rank matrices
  X_mu = make_matrix_full_rank(X_mu)
  X_sigma = make_matrix_full_rank(X_sigma)

  # Give matrices column names if there are none
  names = matrix_column_names( X_mu, X_sigma, intercept_mu, intercept_sigma)
  colnames(X_mu) = names$colnames_X
  colnames(X_sigma) = names$colnames_X_sigma

  # set aliased columns
  names = names(aliased_mu)
  aliased_mu = !is.element( names, colnames(X_mu))
  names(aliased_mu) = names
  names = names(aliased_sigma)
  aliased_sigma = !is.element( names, colnames(X_sigma))
  names(aliased_sigma) = names

  if (nrow(X_mu) < (ncol(X_mu) + ncol(X_sigma))){
    stop("Too few observations. There must be at least ", ncol(X_mu) + ncol(X_sigma), " observations.")
  }

  # Likelihood of model with ordinary linear regression
  fit = stats::lm( y ~ .+0, as.data.frame(as.matrix(X_mu)), model = FALSE)
  sigma2 = fit$residuals %*% fit$residuals / n
  logLik_lm = -n * (log(2*pi*sigma2) + 1) / 2
  logLik_lm = as.numeric(logLik_lm)

  # Set inital estimate of betas if none has been specified
  if ("start" %in% names(slvr_options)){
    slvr_options$start = slvr_options$start[!aliased_sigma]
    slvr_default = list(fn = logLHood)
  }
  else {
    beta = stats::coef(fit)
    sigma = summary(fit)$sigma
    beta_sigma = log(sigma) * solve( Matrix::t(X_sigma) %*% X_sigma) %*% Matrix::colSums(X_sigma)  # Minimize sum over all elements of (X_sigma %*% beta_sigma - log(sigma))^2
    beta_sigma = as.numeric(beta_sigma)
    slvr_default = list( fn = logLHood, start = beta_sigma )
  }

  # Calculate beta for standard deviation sigma
  solve_result = do.call(maxLik::maxNR, c( slvr_default, slvr_options))

  # Check on errors and warnings from solver
  if (!(solve_result$code %in% c( 1, 2, 8))){
    warning(solve_result$message)
  }

  beta_sigma = solve_result$estimate

  # Calculate beta for expectation value mu
  sigma = as.numeric(exp( X_sigma %*% beta_sigma))
  M = X_mu * (1/sigma)
  M = solve(Matrix::t(M) %*% M)
  beta = as.numeric(M %*% Matrix::t(X_mu) %*% (y / (sigma * sigma)))

  # Check Hessian
  if (check_hessian){
    if(!matrixcalc::is.negative.definite(as.matrix(solve_result$hessian))){
      warning("Log-likelihood appears not to be at a maximum!")
    }
  }

  # Calculate log-likelihood
  logLik = solve_result$maximum

  # Associate betas and other result vectors with names of covariates
  names(beta) = colnames(X_mu)
  names(beta_sigma) = colnames(X_sigma)

  # Create ouput list
  rlist = list( call = match.call(),
                coefficients_mu = beta,
                coefficients_sigma = beta_sigma,
                logLik = logLik,
                logLik_lm = logLik_lm,
                aliased_mu = aliased_mu,
                aliased_sigma = aliased_sigma,
                intercept_mu = intercept_mu,
                intercept_sigma = intercept_sigma,
                y = y,
                X_mu = X_mu,
                X_sigma = X_sigma)

  # Add ouput of solver routine
  if (slvr_log){
    rlist$slvr_log = solve_result
  }

  class(rlist) = "lmvar"

  return(rlist)
}
