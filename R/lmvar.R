#' @title Linear regression with non-constant variances
#'
#' @description Performs a Gaussian regression with non-constant variances and outputs an 'lmvar' object.
#'
#' @param y Vector of observations
#' @param X_mu Model matrix for the expected values \eqn{\mu}
#' @param X_sigma Model matrix for the logarithms of the standard deviations \eqn{\sigma}
#' @param slvr_options A list with options to be passed on to the function
#' \code{\link[maxLik]{maxLik}} which maximizes the log-liklihood
#' @param intercept_mu Boolean, if TRUE a column with the name '(Intercept)' is added to the matrix \code{X_mu}. This column
#' represents the intercept term in the model for \eqn{\mu}.
#' @param intercept_sigma Boolean, if TRUE a column with the name '(Intercept_s)' is added to the matrix \code{X_sigma}. This column
#' represents the intercept term in the model for \eqn{\log \sigma}.
#' @param sigma_min Numeric, the minimum value for the standard deviations sigma. Can be a single number which applies to all observations or
#' a vector containing a separate minimum for each observation.
#' @param control A list with further options. The following options can be set.
#' \itemize{
#' \item \code{check_hessian} Boolean, if TRUE it is checked whether the Hessian is negative-definite, i.e., whether the
#' log-likelihood is at a maximum. The default value is TRUE.
#' \item \code{slvr_log} Boolean, if TRUE the output of \code{maxLik} is added to the 'lmvar' object.
#' The default value is FALSE.
#' \item \code{running_diagnostics} Boolean, if TRUE diagnostic messages about errors that occur will be
#' printed during the run. The default value is FALSE.
#' }
#' @param ... Additional arguments, not used in the current implementation
#'
#' @details
#' \subsection{Response vector}{
#' The vector \code{y} must be a numeric vector. It can not contain special values like \code{NULL}, \code{NaN}, etc.
#' }
#' \subsection{Model matrices}{
#' If the matrix \code{X_mu} is not specified, the model for the expected values \eqn{\mu} will consist of an intercept term only.
#' The argument \code{intercept_mu} is ignored in this case. Likewise, the model for \eqn{\log \sigma} will consist of an
#' intercept term only if \code{X_sigma} is not specified.
#' In the latter case, the model reduces to a standard linear model.
#'
#' Both model matrices must be numeric matrices. They can not contain special values like \code{NULL}, \code{NaN}, etc.
#'
#' The model matrices can be of class \code{\link{matrix}} or of a class from the package \code{\link[Matrix]{Matrix}}.
#' They can also be of class \code{\link{numeric}} in case they are matrices with one column only.
#'
#' All columns in the matrix \code{X_mu} must either have a name, or no column has a name at all. It is not allowed that some
#' colums have a
#' name while others don't. The same is true for \code{X_sigma}.
#'
#' If supplied, the column names
#' must be unique within \code{X_mu}. The same is true for \code{X_sigma}. A column name can be used in both \code{X_mu} and
#' \code{X_sigma} though.
#'
#' In case the matrix \code{X_mu} has no columns with a column name, \code{lmvar}  gives the names \code{v1},
#' \code{v2} etc. to the columns. Likewise, if the matrix \code{X_sigma} has no columns with a column name, \code{lmvar}
#' gives the names \code{v1_s}, \code{v2_s} etc. to the columns.
#'
#' Matrix \code{X_mu} can not have a column with the name '(Intercept)'.
#' Matrix \code{X_sigma} can not have a column with the name '(Intercept_s)'. Both names are reserved names.
#' }
#'
#' \subsection{Minimum sigma}{
#' The argument \code{sigma_min} allows one run a regression under the constraint of a minimum standard deviation for each
#' observation. The argument can be a single number, which applies to all observations, or a vector which contains
#' a separate
#' minimum for each observation. In the latter case, all vector elements must be zero
#' (in which case no constraint is applied) or all vector elements must be larger than zero.
#'
#' The option of a minimum sigma is intended for situations in which an unconstrained regression attempts to make
#' sigma equal to zero for one or more observations. This will cause extreme values for the \eqn{\beta_\sigma} and
#' numerical instabilities.Such a situation can be remedied by bounding sigma from below.
#'
#' In case \code{sigma_min} has a value (or a vector of values) larger than zero and the solve-method is "NR", the
#' solve method is set to "BFGS". If the argument \code{constraints} is passed on to \code{maxlik} (as a list member of
#' \code{slvr_options}), it is ignored.
#'
#' Error estimates and confidence intervals (e.g. such as given by \code{\link{summary.lmvar}} and
#' \code{\link{predict.lmvar}}) can be unreliable
#' if minimum sigmas are specified.
#' }
#'
#' \subsection{Likelihood maximization}{
#' The function \code{\link[maxLik]{maxLik}} from the \code{maxLik} package, is used to maximize the (profile) log-likelihood.
#' \code{maxLik} returns a
#' termination code which reports whether a maximum was found, see \code{maxLik}.
#' For the method "NR", a potential problem is reported by a \code{maxLik} return code different from \code{1}, \code{2} or \code{8}.
#' For other methods, a code different from \code{0} flags a potential problem.
#' In case the return code flags a potential problem, the message from \code{maxLik} is output as a warning.
#'
#' All list elements in \code{slvr_options} are passed on as arguments to \code{maxLik}. The name of the list element is
#' the argument name, the value of the list element is the argument value.
#' It is not allowed to pass on the arguments \code{fn}, \code{grad} or \code{hess}. In case the list does not contain an element
#' \code{method}, it is set to "NR". In case the list does not contain an element \code{start}, an initial estimate for
#' \eqn{\beta_\sigma} is set by \code{lmvar}.
#'
#' In case one wants to supply an initial estimate for the coefficients, one has to supply an initial estimate
#' for \eqn{\beta_\sigma}. If \code{beta_sigma_initial} is the initial estimate, one passes on the argument
#' \code{slvr_options = list(start = beta_sigma_initial)}.
#' The inital estimate \code{beta_sigma_initial} must be a numeric vector. Its length must be as follows.
#' \itemize{
#' \item In case \code{X_sigma} is not specified or has the value \code{NULL}, the initial estimate must be
#' a single value.
#' \item In case \code{X_sigma} is specified and \code{intercept_sigma = TRUE}: the length must be equal to the number
#' of columns of \code{X_sigma} plus one. The first element of the vector is the initial estimate of the
#' intercept term for \eqn{\log \sigma}, the second element is the inital estimate corresponding to the first
#' column in \code{X_sigma}, the third element is the inital estimate corresponding to the second
#' column in \code{X_sigma}, etc.
#' \item  In case \code{X_sigma} is specified and \code{intercept_sigma = FALSE}: the length must be equal to the
#' number of columns of \code{X_sigma}. The first element of the vector is the initial estimate corresponding to the
#' first column in \code{X_sigma}, the second element is the inital estimate corresponding to the second
#' column in \code{X_sigma}, etc.
#' }
#'
#' There is no need to supply an inital estimate for \eqn{\beta_\mu}, see the vignette 'Math' for details.
#' }
#' \subsection{Other}{
#' When \code{check_hessian = TRUE}, it is checked whether the fitted log-likelihood is at a maximum. A warning will be issued if
#' that is not the case.
#'
#' See the vignettes
#' that come with the \code{lmvar} package for more info. Run \code{vignette(package="lmvar")} to list the available vignettes.
#' }
#'
#' @return An object of class 'lmvar', which is a list. Users are discouraged to access list-members directly.
#' Instead, list-members are to be accessed with the various accessor and utility functions in the package.
#' Exceptions are the following list members for which no accessor functions exist:
#' \itemize{
#' \item \code{object$y} the vector of observations
#' \item \code{object$X_mu} the  model matrix for \eqn{\mu}. In general, it differs from the user-supplied \code{X_mu} because
#' \code{lmvar} adds an intercept-column and makes the matrix full-rank.
#' \item \code{object$X_sigma} the  model matrix for \eqn{\log \sigma}. In general, it differs from the user-supplied
#' \code{X_sigma} because \code{lmvar} adds an intercept-column and makes the matrix full-rank.
#' \item \code{intercept_mu} boolean which tells whether or not an intercept column \code{(Intercept)} has been added to the
#' model matrix \eqn{X_\mu}
#' \item \code{intercept_sigma} boolean which tells whether or not an intercept column \code{(Intercept_s)} has been added to the
#' model matrix \eqn{X_\sigma}
#' \item \code{sigma_min} the value of the argument \code{sigma_min} in the call of \code{lmvar}
#' \item \code{slvr_options} the value of the argument \code{slvr_options} in the call of \code{lmvar}
#' \item \code{control} the value of the argument \code{control} in the call of \code{lmvar}
#' \item \code{slvr_log} the output of \code{maxLik} (the solver routine used to maximize the likelihood). Included only
#' if the argument \code{slvr_log} has the value \code{TRUE}.
#' See \code{\link[maxLik]{maxLik}} for details about this output.
#' }
#'
#' @export
#'
#' @example R/examples/lmvar_examples.R
#'
lmvar <- function( y, X_mu = NULL, X_sigma = NULL,
                   intercept_mu = TRUE, intercept_sigma = TRUE, sigma_min = 0,
                   slvr_options = list(method = "NR"), control = list(), ...){

  print_diagnostics <- function( beta_sigma, sigma_inv){

    i = which.min(beta_sigma)
    name = colnames(X_sigma)[i]
    value = beta_sigma[i]
    message( "    beta_sigma smallest: '", name, "' = ", signif( value, 3))

    i = which.max(beta_sigma)
    name = colnames(X_sigma)[i]
    value = beta_sigma[i]
    message( "                largest: '", name, "' = ", signif( value, 3))
    cat("\n")

    i = which.min(sigma_inv)
    value = sigma_inv[i]
    message( "    sigma inverse smallest: ", signif( value, 3), " for observation ", i)

    i = which.max(sigma_inv)
    value = sigma_inv[i]
    message( "                   largest: ", signif( value, 3), " for observation ", i)
    cat("\n")
  }

  logLHood <- function(beta_sigma, return_hessian = FALSE){

    sigma_inv = as.numeric(exp(- X_sigma %*% beta_sigma))
    M = X_mu * sigma_inv
    M = tryCatch( chol2inv(chol(Matrix::t(M) %*% M)),          # Exploit that t(M) %*% M is symmetric and positive-definite

                  error = function(e){

                    if (control$running_diagnostics){
                      message (e)
                      cat("\n")
                      print_diagnostics( beta_sigma, sigma_inv)
                    }
                    return(NA)
                  })

    if(class(M) == 'logical'){
      return (NA)
    }

    beta_mu = M %*% Matrix::t(X_mu) %*% (y * sigma_inv^2)
    mu = as.numeric(X_mu %*% beta_mu)
    res = (y - mu) * sigma_inv

    # Calculate log-likelihood
    logLik= -0.5*n*log(2*pi) + sum(log(sigma_inv)) - 0.5*sum(res^2)

    # Calculate gradient
    if(!(slvr_options$method %in% c( "NM", "SANN"))){

      dBeta_sigma = as.numeric(Matrix::t(X_sigma) %*% (res^2 - 1))  # gradient of log-likelihood w.r.t. beta_sigma

      # Set attribute
      attr( logLik, "gradient") = dBeta_sigma
    }

    # Calculate Hessian
    if(slvr_options$method == "NR" | return_hessian){

      M = tryCatch( chol(M),

                    error = function(e){

                      if (control$running_diagnostics){
                        message(e)
                        cat("\n")
                        print_diagnostics( beta_sigma, sigma_inv)
                      }
                      return(NA)
                   })
      if(class(M) == 'logical'){
        return (NA)
      }

      XT = M %*% Matrix::t(X_mu) %*% (X_sigma * (res * sigma_inv)) # Use Choleski decomposition to have H_ss_1 symmetric
      H_ss_1 = 4 * Matrix::t(XT) %*% XT
      XT = X_sigma * abs(res)
      H_ss_2 = -2 * Matrix::t(XT) %*% XT

      # Set attribute
      attr( logLik, "hessian") = as.matrix(H_ss_1 + H_ss_2)
    }

    return(logLik)
  }

  # store argument values
  sigma_min_call = sigma_min
  slvr_options_call = slvr_options
  control_call = control

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

  # set default solve method
  if (!("method" %in% names(slvr_options))){
    slvr_options$method = "NR"
  }

  # Make sure method is upper-case always
  slvr_options$method = toupper(slvr_options$method)

  # set default control options
  if (!("check_hessian" %in% names(control))){
    control$check_hessian = TRUE
  }
  if (!("slvr_log" %in% names(control))){
    control$slvr_log = FALSE
  }
  if (!("running_diagnostics" %in% names(control))){
    control$running_diagnostics = FALSE
  }

  # Check inputs
  if (!isnull_X){
    if (nrow(X_mu) != length(y)){
      stop("Number of rows of matrix X_mu not equal to length of response vector y")
    }
    if (is.element( "(Intercept)", colnames(X_mu))){
      stop("Matrix X_mu can not have a column with the name '(Intercept)'")
    }
  }
  if
  (!isnull_X_sigma){
    if (nrow(X_sigma) != length(y)){
      stop("Number of rows of matrix X_sigma not equal to length of response vector y")
    }
    if (is.element( "(Intercept_s)", colnames(X_sigma))){
      stop("Matrix X_sigma can not have a column with the name '(Intercept_s)'")
    }
  }

  if (length(sigma_min) > 1 & length(sigma_min) != length(y)){
    stop("The length of 'sigma_min' must be 1, or equal to the length of the response vector y")
  }

  if(any(sigma_min > 0) & any(sigma_min <= 0)){
    stop("All values in 'sigma_min' must be strictly positive, or all values must be zero")
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

  # Give matrices column names if there are none
  names = matrix_column_names( X_mu, X_sigma, intercept_mu, intercept_sigma)
  colnames(X_mu) = names$colnames_X
  colnames(X_sigma) = names$colnames_X_sigma

  # Turn matrices into a full-rank matrices
  qX = Matrix::qr(X_mu)
  X_mu = make_matrix_full_rank( qX, X_mu)
  X_sigma = make_matrix_full_rank(X_sigma)

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
  sigma2 = sum(Matrix::qr.resid( qX, y)^2) / n
  logLik_lm = -n * (log(2*pi*sigma2) + 1) / 2

  # Set inital estimate of betas
  if ("start" %in% names(slvr_options)){

    slvr_options$start = slvr_options$start[!aliased_sigma]
    slvr_default_start = list()
  }
  else {

    sigma = sqrt(sigma2)
    if (intercept_sigma){
      beta_sigma = c( log(sigma), rep.int( 0, ncol(X_sigma) - 1))
    }
    else{
      beta_sigma = log(sigma) * solve( Matrix::t(X_sigma) %*% X_sigma) %*% Matrix::colSums(X_sigma)  # Minimize sum over all elements of (X_sigma %*% beta_sigma - log(sigma))^2
    }
    slvr_default_start = list(start = as.numeric(beta_sigma))
  }

  # set function
  slvr_default_fn = list(logLik = logLHood)

  # Set options to solve under boundary condition of minimum value for sigma
  if (all(sigma_min > 0)){
    if (length(sigma_min) == 1){
      sigma_min = rep( sigma_min, times = length(y))
    }
    slvr_options$constraints = list( ineqA = as.matrix(X_sigma), ineqB = - log(sigma_min))

    if (slvr_options$method == "NR"){
      slvr_options$method = "BFGS"
    }
  }

  # Calculate beta for standard deviation sigma
  solve_result = do.call(maxLik::maxLik, c( slvr_default_fn, slvr_default_start, slvr_options))

  # Check on errors and warnings from solver
  if (slvr_options$method == "NR")
    ok_codes = c( 1, 2, 8)
  else {
    ok_codes = 0
    }
  if (!(solve_result$code %in% ok_codes)){
    warning(solve_result$message)
  }

  # Extract beta_sigma from solve results
  beta_sigma = solve_result$estimate

  # Calculate beta for expectation value mu
  sigma = as.numeric(exp( X_sigma %*% beta_sigma))
  M = X_mu * (1/sigma)
  M = solve(Matrix::t(M) %*% M)
  beta = as.numeric(M %*% Matrix::t(X_mu) %*% (y / (sigma^2)))

  # Check Hessian
  if (control$check_hessian){
    if (slvr_options$method == "NR"){
      hessian = maxLik::hessian(solve_result)
    }
    else {
      hessian = attr( logLHood( beta_sigma, return_hessian = TRUE), "hessian")
    }
    if(!matrixcalc::is.negative.definite(hessian)){
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
                y = y,
                X_mu = X_mu,
                X_sigma = X_sigma,
                intercept_mu = intercept_mu,
                intercept_sigma = intercept_sigma,
                sigma_min = sigma_min_call,
                slvr_options = slvr_options_call,
                control = control_call)

  # Add ouput of solver routine
  if (control$slvr_log){
    rlist$slvr_log = solve_result
  }

  class(rlist) = "lmvar"

  return(rlist)
}
