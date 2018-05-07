#' @title Cross-validation for an object of class 'lmvar'
#'
#' @description k-fold cross-validation for an object of class 'lmvar'
#'
#' @param object Object of class 'lmvar'
#' @param k Integer, number of folds
#' @param ks_test Boolean, if \code{TRUE}, a Kolmogorov-Smirnov test is carried out. See details.
#' @param fun User-specified function for which cross-validation results are to be obtained. See details.
#' @param log Boolean, specifies whether \code{object} contains a fit to the response vector \eqn{Y} or its logarithm \eqn{\log Y}
#' @param seed Integer, seed for the random number generator. The seed is not set when \code{seed} equals \code{NULL}.
#' @param sigma_min Minimum value for the standard deviations. Can be a single number which applies to all
#' observations, or a vector giving a minimum per observation. In case of the the default value \code{NULL}, the
#' value is the same as the value in \code{object}.
#' @param exclude Numeric vector with observations that must be excluded for error statistics. The default
#' \code{NULL} means no observations are excluded. See 'Details' for more information.
#' @param slvr_options List of options passed on to the function \code{\link[maxLik]{maxLik}} which carries out the
#' fits for the \eqn{k} folds. See 'Details' for more information.
#' @param max_cores Integer, maximum number of CPU-cores that can be used. For the default value \code{NULL},
#' the number is set to the number of available cores minus one.
#' @param ... Other parameters, not used in the current implementation.
#'
#'
#' @return In case none of the fits in the cross-validations returns an error or a warning, a 'cvlmvar'
#' object is returned. This is a list with the following items:
#' \itemize{
#' \item \code{MAE} a list with two items
#' \itemize{
#' \item \code{mean} the sample mean of the absolute prediction error over the k folds
#' \item \code{sd} the sample standard deviation of the absolute prediction error over the k folds
#' }
#' \item \code{MSE} a list with two items
#' \itemize{
#' \item \code{mean} the sample mean of the mean squared prediction error over the k folds
#' \item \code{sd} the sample standard deviation of the mean squared prediction error over the k folds
#' }
#' \item \code{MSE_sqrt} a list with two items
#' \itemize{
#' \item \code{mean} the sample mean of the root mean squared prediction error over the k folds
#' \item \code{sd} the sample standard deviation of the root mean squared prediction error
#' over the k folds
#' }
#' \item \code{KS_distance} a list with two items
#' \itemize{
#' \item \code{mean} the sample mean of the Kolmogorov-Smirnov distance over the k folds
#' \item \code{sd} the sample standard deviation of the Kolmogorov-Smirnov distance over the k folds
#' }
#' \item \code{KS_p.value} a list with two items
#' \itemize{
#' \item \code{mean} the sample mean of the p-value of Kolmogorov-Smirnov distance over the k folds
#' \item \code{sd} the sample standard deviation of the p-value of the Kolmogorov-Smirnov distance over the k folds
#' }
#' \item \code{fun} a list with two items
#' \itemize{
#' \item \code{mean} the sample mean of the user-specified function \code{fun}
#' \item \code{sd} the sample standard deviation of the of the user-specified function over the k folds
#' }
#' }
#' The items \code{KS_distance} and \code{KS_p.value} are added only in case \code{ks_test = TRUE}.
#'
#' In case a fit returns an error or a warning, the return value of \code{cv.lmvar}
#' lists the arguments of the first call to \code{\link{lmvar}} which failed. In addition, it lists the row
#' number of the observations in \code{object} that formed the training set for which the fit
#' returned an error or warning. These items are returned as a list:
#' \itemize{
#' \item \code{y} the argument \code{y} of the failing call
#' \item \code{X_mu} the argument \code{X_mu} of the failing call
#' \item \code{X_sigma} the argument \code{X_sigma} of the failing call
#' \item \code{intercept_mu} the argument \code{intercept_mu} of the failing call
#' \item \code{intercept_sigma} the argument \code{intercept_sigma} of the failing call
#' \item \code{sigma_min} the argument \code{sigma_min} of the failing call
#' \item \code{slvr_options} the argument \code{slvr_options} of the failing call
#' \item \code{control} the argument \code{control} of the failing call
#' \item \code{training_rows} numeric vector containing the rows of the observations in \code{object} that were
#' used in the failing fit
#' }
#'
#' @details
#' \subsection{Cross-validations}{
#' The function \code{cv.lmvar} carries out a k-fold cross-validation for an 'lmvar' model. For each fold, an 'lmvar'
#' model is fit to all observations that are not in the fold (the 'training set') and prediction errors are calculated
#' for the  observations in the fold (the 'test set'). The prediction errors are the  absolute error \eqn{|y - \mu|}
#' and its square \eqn{(y - \mu)^2}. The average prediction errors over the observations in the fold are calculated,
#' and the square root of the average of the squared errors is taken. Optionally, one can calculate a user-specified
#' function \code{fun} for the test set and the 'lmvar' model resulting from the
#' training set. Optionally, one can also calculate the Kolmogorov-Smirnov (KS) distance for the test set and its p-value.
#'
#' The results for the k folds are averaged over the folds and standard deviations are calculated from the k results.
#' }
#'
#' \subsection{User defined function}{
#' The argument \code{fun} allows a user to specify a function for which cross-validation results
#' must be obtained. This function must meet the following requirements.
#' \itemize{
#' \item Its arguments are:
#' \itemize{
#' \item \code{object_t} an object of class 'lmvar',
#' \item \code{y} a numerical vector of response values and
#' \item \code{X_mu} the model matrix for the expected values of the response vector \code{y}.
#' \item \code{X_sigma} the model matrix for the standard deviations of the response vector \code{y}.
#' }
#' \item It returns a single numerical value.
#' }
#' Carrying out a k-fold cross-validation, the function is called k times with \code{object_t} equal to the fit
#' to the training set, \code{y} equal
#' to the response vector of the test set, and
#' \code{X_mu} and \code{X_sigma} the design matrices of the test set.
#'
#' If the evaluation of \code{fun} gives an error, \code{cv.lmvar} will give a warning and exclude that
#' evaluation from the mean and the standard deviation of \code{fun} over the k folds. If the evaluation
#' of \code{fun} gives a warning, it will be ignored.
#'
#' In the cross-validations, \code{object_t} contains the design matrices of the training set as
#' \code{object_t$X_mu} and \code{object_t$X_sigma}. \code{object_t$X_mu} was formed by taking
#' \code{object$X_mu} and removing the fold-rows. In addition, columns may have been removed to make the matrix
#' full-rank. Therefore, \code{object_t$X_mu} may have fewer columns than \code{object$X_mu}. The same is true
#' for \code{object_t$X_sigma} compared to \code{object$X_sigma}.
#'}
#'
#' \subsection{Kolmogorov-Smirnov test}{
#' When \code{ks_test = TRUE}, a Kolmogorov-Smirnov (KS) test is carried out for each fold. The test checks whether the
#' standardized residuals \eqn{(y - \mu) / \sigma} in a fold are distributed as a standard normal distribution. The
#' KS-distance and the corresponding p-value are calculated for each fold. The test uses the
#' function \code{\link[stats]{ks.test}}. The expectation values \eqn{\mu} and standard deviations \eqn{\sigma} are
#' calculated from the model matrices for the test set (the fold) and the 'lmvar' fit to the training set.
#' }
#'
#' \subsection{Excluding observations}{
#' The observations specified in the argument \code{exclude} are not used to calculate the error statistics MAE
#' (mean absolute error), MSE (mean squared error) and the square root of MSE. They are also not used to calculate
#' the statistics for the user-defined function \code{fun}. This is useful when there are a few observations
#' with such large residuals that they dominate the error estimates. Note that the excluded observations are not
#' excluded from the training sets. It is only in the calculation of the statistics of the test sets that the
#' observations are
#' excluded. They are not excluded from the KS-test: when observations have large residuals, they should have large
#' standard deviations as well,
#' to give the standardized residuals normal values.
#' }
#'
#' \subsection{Minimum sigma}{
#' The argument \code{sigma_min} gives the option to enforce a minimum standard deviation. This is
#' useful when, in a cross-validation, a fit fails because the maximum likelihood occurs when the standard
#' deviation of one or more observations becomes zero.
#' When a minimum  standard deviation is specified, all fits are carried out under the
#' boundary condition that the standard deviation is larger than the minimum. If \code{sigma_min = NULL} the same value
#' is used as was used to create \code{object}.
#' }
#'
#' \subsection{Other}{
#' The fits are carried out with the options \code{slvr_options} stored in the 'lmvar' object \code{object}.
#' However, these options can be overwritten with an explicit argument \code{slvr_options} in the call of
#' \code{cv.lmvar}. Some of the options are affected by a \code{sigma_min} larger than zero, see \code{\link{lmvar}} for
#' details.
#'
#' The argument \code{slvr_options} is a list, members of which can be a list themselves.
#' If  members of a sublist are overwritten, the other members of the sublist remain unchanged. E.g., the
#' argument \code{slvr_options = list(control = list(iterlim = 600))} will set \code{control$iterlim} to 600 while
#' leaving other members of the list \code{control} unchanged.
#'
#' The number of available CPU cores is detected with \code{\link[parallel]{detectCores}}.
#' }
#'
#' @seealso See \code{\link{lmvar}} for the options \code{slvr_options} stored in an 'lmvar' object.
#'
#' \code{\link{cv.lm}} is the equivalent function for an object of class 'lm'. It is supplied in case one wants to
#' compare an 'lmvar' fit with an 'lm' fit.
#'
#' \code{\link{print.cvlmvar}} provides a print-method for an object of class 'cvlmvar'.
#'
#' @export
#'
#' @example R/examples/cv.lmvar_examples.R
#'
cv.lmvar <- function( object, k = 10, ks_test = FALSE, fun = NULL, log = FALSE, seed = NULL, sigma_min = NULL, exclude = NULL,
                      slvr_options = list(), max_cores = NULL, ...){

  # make matrix a matrix again
  matrix_as_matrix <- function( X, col_names){
    if (inherits(X, "numeric")){
      X = as.matrix(X)
      colnames(X) = col_names
    }
    return(X)
  }

  # Check input
  if (!inherits( object, "lmvar")){
    stop("Object must be of type 'lmvar'")
  }

  size = trunc(nobs(object) / k)
  if (size == 0){
    stop("Number of folds is too large for the number of observations in 'object'")
  }

  # set minimum sigma
  if(is.null(sigma_min)){
    sigma_min = object$sigma_min
  }

  if (length(sigma_min) > 1 & length(sigma_min) != nobs(object)){
    stop("The length of 'sigma_min' must be 1, or equal to the number of observations in 'object'")
  }

  if(any(sigma_min > 0) & any(sigma_min <= 0)){
    stop("All values in 'sigma_min' must be strictly positive, or all values must be zero")
  }

  # Set RNG seed
  if (!is.null(seed)){
    set.seed(seed)
  }

  # Set whether or not function is present
  isFunc = !is.null(fun)

  # Retrieve info from object
  y = object$y
  intercept_mu = object$intercept_mu
  intercept_sigma = object$intercept_sigma
  slvr_options_object = object$slvr_options
  control = object$control

  # Set excluded rows
  excluded_rows = is.element( 1:length(y), exclude)

  # Set solver options
  for (option in names(slvr_options_object)){
    if (!(option %in% names(slvr_options))){
      slvr_options[option] = slvr_options_object[option]
    }
    else if(inherits( slvr_options_object[option], "list")){
      for (suboption in names(slvr_options_object[[option]])){
        if (!(suboption %in% names(slvr_options[[option]]))){
          slvr_options[[option]][suboption] = slvr_options_object[[option]][suboption]
        }
      }
    }
  }

  # Set beta_sigma
  beta_sigma = coef.lmvar( object, mu = FALSE)

  # set the folds
  remaining = seq.int( from = 1, to = length(y))
  selected_obs = matrix(nrow = k, ncol = size)
  for (i in 1:k){
    selected_obs[i,] = sample( remaining, size)
    remaining = setdiff( remaining, selected_obs[i,])
  }

  # set up cluster of cores
  no_cores = parallel::detectCores() - 1
  no_cores = max( no_cores, 1)
  if (!is.null(max_cores)){
    no_cores = min( max_cores, no_cores, k)
  }
  cl = parallel::makeCluster(no_cores)

  # prepare cores
  parallel::clusterEvalQ( cl, library(Matrix))

  # loop over cross-validations
  cv_results = parallel::parLapply( cl, 1:k, function(i){

    # select elements from response vector
    foldrows = is.element(1:length(y), selected_obs[i,])

    # create corresponding model matrices
    X_mu = object$X_mu[!foldrows,]
    if(inherits( X_mu, "numeric")){
      X_mu = matrix_as_matrix( X_mu, colnames(object$X_mu))
    }

    if(intercept_mu){
      if(ncol(X_mu) > 1){
        X_mu = X_mu[,-1]
        if(inherits( X_mu, "numeric")){
          X_mu = matrix_as_matrix( X_mu, colnames(object$X_mu)[-1])
        }
      }
      else {
        X_mu = NULL
      }
    }

    X_sigma = object$X_sigma[!foldrows,]
    if(inherits( X_sigma, "numeric")){
      X_sigma = matrix_as_matrix( X_sigma, colnames(object$X_sigma))
    }

    if (intercept_sigma){
      if(ncol(X_sigma) > 1){
        X_sigma = X_sigma[,-1]
        if(inherits( X_sigma, "numeric")){
          X_sigma = matrix_as_matrix( X_sigma, colnames(object$X_sigma)[-1])
        }
      }
      else {
        X_sigma = NULL
      }
    }

    # model matrices are not guaranteed to be full rank
    control$mu_full_rank = FALSE
    control$sigma_full_rank = FALSE

    # perform fit on rows not in fold
    fit = tryCatch( lmvar( y[!foldrows], X_mu, X_sigma,
                           intercept_mu = intercept_mu, intercept_sigma = intercept_sigma,
                           sigma_min = sigma_min, slvr_options = slvr_options,
                           control = control),

                    error = function(e){
                      outlist = list( error = e, warning = NULL,
                                      y = y[!foldrows], X_mu = X_mu, X_sigma = X_sigma,
                                      intercept_mu = intercept_mu, intercept_sigma = intercept_sigma,
                                      sigma_min = sigma_min, slvr_options = slvr_options, control = control,
                                      training_rows = which(!foldrows))
                      class(outlist) = "error_fit"
                      return(outlist)
                    },
                    warning = function(w){
                      outlist = list( error = NULL, warning = w,
                                      y = y[!foldrows], X_mu = X_mu, X_sigma = X_sigma,
                                      intercept_mu = intercept_mu, intercept_sigma = intercept_sigma,
                                      sigma_min = sigma_min, slvr_options = slvr_options, control = control,
                                      training_rows = which(!foldrows))
                      class(outlist) = "error_fit"
                      return(outlist)
                    })
    if (inherits( fit, "error_fit")){
      return(fit)
    }

    # create prediction matrices
    X_mu = object$X_mu[foldrows & !excluded_rows,]
    if (inherits( X_mu, "numeric")){
      X_mu = matrix_as_matrix( X_mu, colnames(object$X_mu))
    }

    if (intercept_mu){
      if(ncol(X_mu) > 1){
        X_mu = X_mu[,-1]
        if (inherits( X_mu, "numeric")){
          X_mu = matrix_as_matrix( X_mu, colnames(object$X_mu)[-1])
        }
      }
      else {
        X_mu = matrix( 1, nrow = sum(foldrows), ncol = 1)
      }
    }

    X_sigma = object$X_sigma[foldrows & !excluded_rows,]
    if (inherits( X_sigma, "numeric")){
      X_sigma = matrix_as_matrix( X_sigma, colnames(object$X_sigma))
    }

    if (intercept_sigma){
      if(ncol(X_sigma) > 1){
        X_sigma = X_sigma[,-1]
        if (inherits( X_sigma, "numeric")){
          X_sigma = matrix_as_matrix( X_sigma, colnames(object$X_sigma)[-1])
        }
      }
      else {
        X_sigma = matrix( 1, nrow = sum(foldrows), ncol = 1)
      }
    }

    # predict values for rows in fold
    predictions = predict.lmvar( fit, X_mu, X_sigma, log = FALSE)
    mu_not_log = predictions[,"mu"]
    sigma_not_log = predictions[,"sigma"]

    if (log){
      predictions = predict.lmvar( fit, X_mu, X_sigma, log = log)
      mu = predictions[,"mu"]
      sigma = predictions[,"sigma"]
    }
    else {
      mu = mu_not_log
      sigma = sigma_not_log
    }

    # Calculate MAE and MSE for rows in fold
    if (log){
      mae = abs(exp(y[foldrows & !excluded_rows]) - mu)
    }
    else {
      mae = abs(y[foldrows & !excluded_rows] - mu)
    }
    mse = mae^2

    # Calculate Kolmogorov-Smirnov distance
    if (ks_test){
      z = (y[foldrows] - mu_not_log) / sigma_not_log
      ks = tryCatch(stats::ks.test( z, "pnorm"),
                    warning = function(w){
                      warning(w)

                      # Determine duplicates
                      dup = duplicated(z)
                      if (sum(dup) > 0){
                        dup = which(dup)[1]
                        dup = which(z == z[dup])

                        # Convert duplicate indices to observation indices
                        dup = which(foldrows)[dup]

                        # Re-run KS-test
                        ks = suppressWarnings(stats::ks.test( z, "pnorm"))

                        # Return warning with duplicates
                        outlist = list( warning = w, duplicates = dup,
                                        statistic = ks$statistic, p.value = ks$p.value)
                        class(outlist) = "error_ks"
                        return(outlist)
                      }
                    })
    }

    # Calculate user-specified function
    if (isFunc){
      f_user = tryCatch(fun( fit, y[foldrows & !excluded_rows],
                             object$X_mu[foldrows & !excluded_rows,],
                             object$X_sigma[foldrows & !excluded_rows,]),
                        error = function(e){
                          outlist = list(error = e)
                          class(outlist) = "error_func"
                          return(outlist)
                        })
    }

    # Return results
    outlist = list( mae = mean(mae), mse = mean(mse), mse_sqrt = sqrt(mean(mse)))
    if (ks_test){
      outlist$ks = ks
    }
    if (isFunc){
      outlist$fun = f_user
    }
    return(outlist)
  })

  # release cluster of cores
  parallel::stopCluster(cl)

  # manipulate list with cross-validation results and print warnings
  cv_results = lapply( cv_results, function(cv_result){

    # Handle failing fit
    if (inherits( cv_result, "error_fit")){
      if (is.null(cv_result$warning)){
        warning(cv_result$error)
      }
      else {
        warning(cv_result$warning)
      }
    }
    else {

      # unpack the list 'ks'
      if ("ks" %in% names(cv_result)){

        # handle warnings
        if (inherits( cv_result$ks, "error_ks")){
          warning(cv_result$ks$warning)
          dup = cv_result$ks$duplicates
          warning("  Observations ", dup[1], " and ", dup[2],
                  " (and maybe others) have identical standardized residuals", call. = FALSE)
        }
        cv_result$ks_distance = cv_result$ks$statistic
        cv_result$ks_p.value = cv_result$ks$p.value
        cv_result$ks = NULL
      }

      if ("fun" %in% names(cv_result)){
        if (inherits( cv_result$fun, "error_func")){
          warning( "Error in user-specified function: ", cv_result$fun$error, call. = FALSE)
          cv_result$fun= NA
        }
      }
    }
    return(cv_result)
  })

  # Check whether any fit failed
  classes = sapply( cv_results, function(cv_result){
    return(class(cv_result))
  })

  if (any(classes == "error_fit")){

    # Return info about failed fit
    index = which(classes == "error_fit")[1]
    outlist = cv_results[[index]]
    outlist$error = NULL
    outlist$warning = NULL
    return(outlist)
  }
  else {

    # Create list with outputs
    outlist = list()
    for (statistic in names(cv_results[[1]])){

      # Store cross-validation results for statistic in vector
      vec = sapply( cv_results, function(cv_result){
        return(cv_result[[statistic]])
      })

      statistic_value = list( mean = mean( vec, na.rm = TRUE), sd = stats::sd( vec, na.rm = TRUE))

      # set the name of the statistic in the 'cvlmvar' object
      name = statistic
      if (statistic %in% c( "mae", "mse")){
        name = toupper(name)
      }
      else if (name == "mse_sqrt"){
        name = "MSE_sqrt"
      }
      else if (name == "ks_distance"){
        name = "KS_distance"
      }
      else if (name == "ks_p.value"){
        name = "KS_p.value"
      }
      outlist[[name]] = statistic_value
    }

    class(outlist) = "cvlmvar"

    return( outlist)
  }
}
