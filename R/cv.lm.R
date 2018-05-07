#' @title Cross-validation for an object of class 'lm'
#'
#' @description k-fold cross-validation for an object of class 'lm'
#'
#' @param object Object of class 'lm'
#' @param k Integer, number of folds
#' @param ks_test Boolean, if \code{TRUE}, a Kolmogorov-Smirnov test is carried out. See details.
#' @param fun User-specified function for which cross-validation results are to be obtained. See details.
#' @param log Boolean, specifies whether \code{object} contains a fit to the response vector \eqn{Y} or its logarithm \eqn{\log Y}
#' @param seed Integer, seed for the random number generator. The seed is not set when \code{seed} equals \code{NULL}.
#' @param max_cores Integer, maximum number of CPU-cores that can be used. For the default value \code{NULL},
#' the number is set to the number of available cores minus one.
#' @param ... Other parameters, not used in the current implementation.
#'
#'
#' @return An object of class 'cvlmvar', which is a list with the following items:
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
#' The items \code{KS_distance} and \code{KS_p.value} are added only in case \code{ks_test = TRUE}. The item
#' \code{fun} is added only in case a function \code{fun} has been specified.
#'
#' @details
#' \subsection{Cross-validations}{
#' The function \code{cv.lm} carries out a k-fold cross-validation for a linear model (i.e. a 'lm' model).
#' For each fold, an 'lm'
#' model is fit to all observations that are not in the fold (the 'training set') and prediction errors are calculated
#' for the  observations in the fold (the 'test set'). The prediction errors are the  absolute error \eqn{|y - \mu|}
#' and its square \eqn{(y - \mu)^2}. The average prediction errors over the observations in the fold are calculated,
#' and the square root of the average of the squared errors is taken. Optionally, one can calculate a user-specified
#' function \code{fun} for the test set and the 'lmvar' model resulting from the
#' training set. Optionally, one can also calculate the Kolmogorov-Smirnov (KS) distance for the test set and its p-value.
#'
#' The results for the k folds are averaged over the folds and standard deviations are calculated from the k results.
#' }
#' \subsection{Requirements on the 'lm' object}{
#' \code{object} must contain the list-members \code{x} and \code{y}. I.e., it must be created by running
#' \code{\link[stats]{lm}} with the options \code{x = TRUE} and \code{y = TRUE}.
#' }
#'
#' \subsection{User defined function}{
#' The argument \code{fun} allows a user to specify a function for which cross-validation results
#' must be obtained. This function must meet the following requirements.
#' \itemize{
#' \item Its arguments are:
#' \itemize{
#' \item \code{object_t} an object of class 'lm',
#' \item \code{y} a numerical vector of response values and
#' \item \code{X} the model matrix for the response vector \code{y}.
#' }
#' \item It returns a single numerical value.
#' }
#' Carrying out a k-fold cross-validation, the function is called k times with \code{object_t} equal to the fit
#' to the training set, \code{y} equal
#' to the response vector of the test set, and
#' \code{X_mu} the design matrix of the test set.
#'
#' If the evaluation of \code{fun} gives an error, \code{cv.lm} will give a warning and exclude that
#' evaluation from the mean and the standard deviation of \code{fun} over the k folds. If the evaluation
#' of \code{fun} gives a warning, it will be ignored.
#'
#' In the cross-validations, \code{object_t} contains the design matrix used in the fit to the training set as
#' \code{object_t$x}.
#'}
#'
#' \subsection{Kolmogorov-Smirnov test}{
#' When \code{ks_test = TRUE}, a Kolmogorov-Smirnov (KS) test is carried out for each fold. The test checks whether the
#' standardized residuals \eqn{(y - \mu) / \sigma} in a fold are distributed as a standard normal distribution. The
#' KS-distance and the corresponding p-value are calculated for each fold. The test uses the
#' function \code{\link[stats]{ks.test}}. The expectation values \eqn{\mu} and standard deviation \eqn{\sigma} are
#' calculated from the model matrices for the test set (the fold) and the 'lm' fit to the training set.
#' }
#'
#' \subsection{Other}{
#' The number of available CPU cores is detected with \code{\link[parallel]{detectCores}}.
#' }
#'
#' @seealso \code{\link{cv.lmvar}} is the equivalent function for an object of class 'lmvar'. It is supplied in
#' case one wants to compare an 'lmvar' fit with an 'lm' fit.
#'
#' \code{\link{print.cvlmvar}} provides a print-method for an object of class 'cvlmvar'.
#'
#' @export
#'
#' @example R/examples/cv.lm_examples.R
#'
cv.lm <- function( object, k = 10, ks_test = FALSE, fun = NULL, log = FALSE, seed = NULL, max_cores = NULL, ...){

  # Check input
  if (!inherits( object, "lm")){
    stop("Object must be of type 'lm'")
  }

  if (!('y' %in% names(object))){
    stop( "Object must contain a list member 'y' containing the response vector. Please run 'lm' with 'y = TRUE'")
  }
  if (!('x' %in% names(object))){
    stop( "Object must contain a list member 'x' containing the model matrix. Please run 'lm' with 'x = TRUE'")
  }

  size = trunc(nobs(object) / k)
  if (size == 0){
    stop("Number of folds is too large for the number of observations in 'object'")
  }

  # Set RNG seed
  if (!is.null(seed)){
    set.seed(seed)
  }

  # Set whether or not function is present
  isFunc = !is.null(fun)

  # Retrieve info from object
  y = object$y

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

  # loop over cross-validations
  cv_results = parallel::parLapply( cl, 1:k, function(i){

    # select elements from response vector
    foldrows = is.element(1:length(y), selected_obs[i,])

    # create model matrix
    XX = object$x[!foldrows,]
    if (inherits( XX, "numeric")){
      XX = as.matrix(XX)
      colnames(XX) = colnames(object$x)
    }
    XX = make_matrix_full_rank(XX)

    # perform fit on rows not in fold
    fit = stats::lm( y[!foldrows] ~ . - 1, as.data.frame(XX), x = TRUE)

    # predict values for rows in fold
    cols = colnames(XX)
    XX = object$x[ foldrows, cols]
    if (inherits( XX, "numeric")){
      XX = as.matrix(XX)
      colnames(XX) = cols
    }

    mu_not_log = stats::predict( fit, as.data.frame(XX))
    sigma_not_log = summary(fit)$sigma

    if (log){
      mu = exp(mu_not_log + sigma_not_log^2 / 2)
      sigma = mu * sqrt(exp(sigma_not_log^2) - 1)
    }
    else {
      mu = mu_not_log
      sigma = sigma_not_log
    }

    # Calculate MAE and MSE for rows in fold
    if (log){
      mae = abs(exp(y[foldrows]) - mu)
    }
    else {
      mae = abs(y[foldrows] - mu)
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
      f_user = tryCatch(fun( fit, y[foldrows], XX),
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
    return(cv_result)
  })

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
