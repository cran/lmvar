#' @title Create an 'lmvar'-like object without a model fit
#'
#' @description Creates an 'lmvar'-like object without carrying out a model fit. This object is a 'lmvar' object from which all members
#' have been left out
#' that are the result of the fit. Such an object can be used in functions which typically use an 'lmvar' object as input but do not
#' need the fit results. Since no fit is performed, \code{lmvar_no_fit} is convenient when the fit is time-consuming or, e.g.,
#' does not converge.
#'
#' @param y Vector of observations
#' @param X_mu Model matrix for the expected values \eqn{\mu}
#' @param X_sigma Model matrix for the logarithms of the standard deviations \eqn{\sigma}
#' @param slvr_options see the function \code{lmvar}
#' @param intercept_mu see the function \code{\link{lmvar}}
#' @param intercept_sigma see the function \code{lmvar}.
#' @param sigma_min see the function \code{lmvar}.
#' @param control see the function \code{lmvar}.
#' @param ... Additional arguments, not used in the current implementation
#'
#' @details
#' See \code{\link{lmvar}} for the requirements and a further explanation of all the arguments.
#'
#' The class 'lmvar' is an extension of the class 'lmvar_no_fit'. This means that each object which is of class 'lmvar', is  of class
#' 'lmvar_no_fit' as well. Wherever an object of class 'lmvar_no_fit' is required, an object of class 'lmvar' can be used as well.
#'
#' Accessor and utility functions for a 'lmvar_no_fit' object are: \code{\link{nobs.lmvar_no_fit}}, \code{\link{alias.lmvar_no_fit}}
#' and \code{\link{dfree}}
#'
#' The function \code{lmvar_no_fit} is especially useful in combination with \code{\link{fwbw.lmvar_no_fit}}. In case a
#' model with many degrees of freedom does not converge with \code{\link{lmvar}}, one can create an 'lmvar_no_fit' object. This
#' is used as
#' input for \code{fwbw} with the argument \code{fw = TRUE}. The \code{fwbw} algorithm will try to select an optimal subset of all
#' degrees of freedom, starting with the smallest subset possible.
#'
#' Although \code{lmvar_no_fit} does not carry out a model fit, it will do the following:
#' \itemize{
#' \item add an intercept term to the model matrices (unless \code{intercept_mu} is FALSE and/or
#' \code{intercept_sigma} is FALSE), and
#' \item make the model matrices full rank.
#' }
#'
#' @return An object of class 'lmvar_no_fit', which is a list. The list-members are the same as for an object of call 'lmvar', but with
#' members that are the result of the model fit left out.
#'
#' Users are discouraged to access list-members directly.
#' Instead, list-members are to be accessed with the various accessor and utility functions in the package.
#' Exceptions are the following list members for which no accessor functions exist:
#' \itemize{
#' \item \code{call} the function call
#' \item \code{y} the vector of observations
#' \item \code{X_mu} the  model matrix for \eqn{\mu}. In general, it differs from the user-supplied \code{X_mu} because
#' \code{lmvar_no_fit} adds an intercept-column (unless \code{intercept_mu} is FALSE) and makes the matrix full-rank.
#' \item \code{X_sigma} the  model matrix for \eqn{\log \sigma}. In general, it differs from the user-supplied
#' \code{X_sigma} because \code{lmvar_no_fit} adds an intercept-column (unless \code{intercept_sigma} is FALSE) and makes
#' the matrix full-rank.
#' \item \code{intercept_mu} boolean which tells whether or not an intercept column \code{(Intercept)} has been added to the
#' model matrix \eqn{X_\mu}
#' \item \code{intercept_sigma} boolean which tells whether or not an intercept column \code{(Intercept_s)} has been added to the
#' model matrix \eqn{X_\sigma}
#' \item \code{sigma_min} the value of the argument \code{sigma_min} in the call of \code{lmvar_no_fit}
#' \item \code{slvr_options} the value of the argument \code{slvr_options} in the call of \code{lmvar_no_fit}
#' \item \code{control} the value of the argument \code{control} in the call of \code{lmvar_no_fit}
#' }
#'
#' @export
#'
#' @example R/examples/lmvar_no_fit_examples.R
#'
lmvar_no_fit <- function( y, X_mu = NULL, X_sigma = NULL,
                   intercept_mu = TRUE, intercept_sigma = TRUE, sigma_min = 0,
                   slvr_options = list(method = "NR"), control = list(), ...){

  # store argument values
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
  if (inherits( X_mu, 'numeric')){
    X_mu = as.matrix(X_mu)
  }
  if (inherits( X_sigma, 'numeric')){
    X_sigma = as.matrix(X_sigma)
  }

  # set default control options
  if (!("mu_full_rank" %in% names(control))){
    control$mu_full_rank = FALSE
  }
  if (!("sigma_full_rank" %in% names(control))){
    control$sigma_full_rank = FALSE
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

  # Give matrices column names if there are none
  names = matrix_column_names( X_mu, X_sigma, intercept_mu, intercept_sigma)
  colnames(X_mu) = names$colnames_X
  colnames(X_sigma) = names$colnames_X_sigma

  # Initialize aliased columns
  aliased_mu = logical(length = ncol(X_mu))
  aliased_sigma = logical(length = ncol(X_sigma))
  names(aliased_mu) = colnames(X_mu)
  names(aliased_sigma) = colnames(X_sigma)

  # Turn matrices into a full-rank matrices
  equal_matrices = FALSE
  if (!control$sigma_full_rank & ncol(X_mu) == ncol(X_sigma)){
    if (all(X_mu == X_sigma)){
      equal_matrices = TRUE
    }
  }

  qX = Matrix::qr(X_mu)
  if (!control$mu_full_rank){
    X_mu = make_matrix_full_rank( qX, X_mu)
  }
  if (!control$sigma_full_rank){
    if (equal_matrices){
      i = names(aliased_mu) %in% colnames(X_mu)
      X_sigma = X_mu
      colnames(X_sigma) = names(aliased_sigma)[i]
    }
    else X_sigma = make_matrix_full_rank(X_sigma)
  }

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

  # Create ouput list
  rlist = list( call = match.call(),
                aliased_mu = aliased_mu,
                aliased_sigma = aliased_sigma,
                y = y,
                X_mu = X_mu,
                X_sigma = X_sigma,
                intercept_mu = intercept_mu,
                intercept_sigma = intercept_sigma,
                sigma_min = sigma_min,
                slvr_options = slvr_options,
                control = control_call)

  class(rlist) = "lmvar_no_fit"

  return(rlist)
}
