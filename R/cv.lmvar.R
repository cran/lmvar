#' @title Cross-validation for an object of class 'lmvar'
#'
#' @description k-fold cross-validation for an object of class 'lmvar'
#'
#' @param object Object of class 'lmvar'
#' @param k Integer, number of folds
#' @param ks_test Boolean, if \code{TRUE}, a Kolmogorov-Smirnov test is carried out. See details.
#' @param log Boolean, specifies whether \code{object} contains a fit to the response vector \eqn{Y} or its logarithm \eqn{\log Y}
#' @param seed Integer, seed for the random number generator. The seed is not set when \code{seed} equals \code{NULL}.
#' @param sigma_min Minimum value for the standard deviations. Can be a single number which applies to all
#' observations, or a vector giving a minimum per observation.
#' @param exclude Numeric vector with observations that must be excluded for error statistics. The default
#' \code{NULL} means no observations are excluded. See 'Details' for more information.
#' @param slvr_options List of options passed on to the function \code{\link[maxLik]{maxLik}} which carries out the
#' fits for the \eqn{k} folds. See 'Details' for more information.
#' @param ... Other parameters, not used in the current implementation.
#'
#'
#' @return In case none of the fits in the cross-validations returns an error, a 'cvlmvar' object is returned.
#' This is a list with the following items:
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
#' \item \code{mean} the sample mean of the square root of the mean squared prediction error over the k folds
#' \item \code{sd} the sample standard deviation of the square root of the mean squared prediction error
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
#' }
#' The items \code{KS_distance} and \code{KS_p.value} are added only in case \code{ks_test = TRUE}.
#'
#' In case a fit returns an error, the return value of \code{cv.lmvar}
#' lists the arguments of the first call to \code{\link{lmvar}} which failed. In addition, it lists the rows of the
#' full set of observations in \code{object} that were used in the fit. These items are returned as a list:
#' \itemize{
#' \item \code{y} the argument \code{y} of the failing call
#' \item \code{X_mu} the argument \code{X_mu} of the failing call
#' \item \code{X_sigma} the argument \code{X_sigma} of the failing call
#' \item \code{intercept_mu} the argument \code{intercept_mu} of the failing call
#' \item \code{intercept_sigma} the argument \code{intercept_sigma} of the failing call
#' \item \code{slvr_options} the argument \code{slvr_options} of the failing call
#' \item \code{control} the argument \code{control} of the failing call
#' \item \code{training_rows} numeric vector containing the rows of all observations in \code{object} that were
#' used in the failing fit
#' }
#'
#' @details When \code{ks_test = TRUE}, a Kolmogorov-Smirnov (KS) test is carried out. The test checks whether the
#' standardized residuals \eqn{(y - \mu) / \sigma} in a fold are distributed as a standard normal distribution. The
#' KS-distance and the corresponding p-value are calculated for each fold. The test uses the function
#' \code{\link[stats]{ks.test}}.
#'
#' The argument \code{sigma_min} gives the option to enforce a minimum standard deviation. This is
#' useful when, in a cross-validation, a fit fails because the maximum likelihood occurs when the standard
#' deviation of one or more observations becomes zero.
#' When a minimum  standard deviation is specified, all fits are carried out under the
#' boundary condition that the standard deviation is larger than the minimum.
#'
#' The fits are carried out with the options \code{slvr_options} stored in the 'lmvar' object \code{object}.
#' However, these options can be overwritten with an explicit argument \code{slvr_options} in the call of
#' \code{cv.lmvar}.
#'
#' The argument \code{slvr_options} is a list, members of which can be a list themselves.
#' If  members of a sublist are overwritten, the other members of the sublist remain unchanged. E.g., the
#' argument \code{slvr_options = list(control = list(iterlim = 600))} will set \code{control$iterlim} to 600 while
#' leaving other members of the list \code{control} unchanged.
#'
#' If a minimum value larger than zero for the standard deviation is specified and the
#' solver-option \code{method} is \code{"NR"}, the solve method used is \code{"BFGS"} and not \code{"NR"}. Also, the
#' solver-option \code{constraints} (if present) is overridden when a minimum value larger than zero for
#' the standard deviation is specified.
#'
#' The observations specified in the argument \code{exclude} are not used to calculate the error statistics MAE
#' (mean absolute error),
#' MSE (mean squared error) and the square root of MSE. This is useful when there are a few observations
#' with such large residuals that they dominate these error estimates (these observations can have
#' very large estimates for the standard deviation as well,
#' in which case the standardized residuals have normal values). Note that the excluded observations are used
#' in the fits carried out in the cross-validation. It is only in the calculation of the statistics that they are
#' excluded. They are not excluded from the KS-test.
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
cv.lmvar <- function( object, k = 10, ks_test = FALSE, log = FALSE, seed = NULL, sigma_min = 0, exclude = NULL,
                      slvr_options = list(), ...){

  # Check input
  if (!any(class(object) == "lmvar")){
    stop("Object must be of type 'lmvar'")
  }

  size = trunc(nobs(object) / k)
  if (size == 0){
    stop("Number of folds is too large for the number of observations in 'object'")
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

  # Set excluded rows
  excluded_rows = is.element( 1:length(y), exclude)

  # Retrieve info from object
  y = object$y
  X_mu = object$X_mu
  X_sigma = object$X_sigma
  intercept_mu = object$intercept_mu
  intercept_sigma = object$intercept_sigma
  slvr_options_object = object$slvr_options
  control = object$control

  # Set solver options
  for (option in names(slvr_options_object)){
    if (!(option %in% names(slvr_options))){
      slvr_options[option] = slvr_options_object[option]
    }
    else if(class(slvr_options_object[option]) == 'list'){
      for (suboption in names(slvr_options_object[[option]])){
        if (!(suboption %in% names(slvr_options[[option]]))){
          slvr_options[[option]][suboption] = slvr_options_object[[option]][suboption]
        }
      }
    }
  }

  # Convert minimim sigma to vector
  if (length(sigma_min) == 1){
    sigma_min = rep( sigma_min, times = length(y))
  }

  # Set options to solve under boundary condition of minimum value for sigma
  constraint = FALSE
  if (all(sigma_min > 0)){

    constraint = TRUE
    if (slvr_options$method == "NR"){
      slvr_options$method = "BFGS"
    }
    ineqA = as.matrix(X_sigma)
    ineqB = - log(sigma_min)
  }

  # Remove intercept terms from model matrices
  if (intercept_mu){
    if (ncol(X_mu) > 1){
      X_mu = X_mu[,-1]
      if(class(X_mu) == "numeric"){
        X_mu = as.matrix(X_mu)
      }
    }
    else {
      X_mu = NULL
    }
  }
  if (intercept_sigma){
    if (ncol(X_sigma) > 1){
      X_sigma = X_sigma[,-1]
      if(class(X_sigma) == "numeric"){
        X_sigma = as.matrix(X_sigma)
      }
    }
    else {
      X_sigma = NULL
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

  # loop over cross-validations
  cv_results = lapply( 1:k, function(i){

    # select elements from response vector
    foldrows = is.element(1:length(y), selected_obs[i,])

    # create corresponding model matrices
    if (is.null(X_mu)){
      XX = NULL
      XX_predict = NULL
    }
    else {
      XX = X_mu[!foldrows,]
      XX_predict = X_mu[foldrows,]
    }

    if (is.null(X_sigma)){
      XX_sigma = NULL
      XX_sigma_predict = NULL
    }
    else {
      XX_sigma = X_sigma[ !foldrows,]
      XX_sigma_predict = X_sigma[foldrows,]
    }

    # create constrains
    if (constraint){
      slvr_options$constraints = list( ineqA = ineqA[!foldrows,], ineqB = ineqB[!foldrows])
    }

    # perform fit on rows not in fold
    fit = tryCatch( lmvar( y[!foldrows], XX, XX_sigma,
                           intercept_mu = intercept_mu, intercept_sigma = intercept_sigma,
                           slvr_options = slvr_options,
                           control = control),

                    error = function(e){
                      warning(e)
                      outlist = list( y = y[!foldrows], X_mu = XX, X_sigma = XX_sigma,
                                      intercept_mu = intercept_mu, intercept_sigma = intercept_sigma,
                                      slvr_options = slvr_options, control = control,
                                      training_rows = which(!foldrows))
                      return(outlist)
                    })
    if (class(fit) == "list"){
      return(fit)
    }

    # predict values for rows in fold
    predictions = predict.lmvar( fit, XX_predict, XX_sigma_predict, log = FALSE)
    mu_not_log = predictions[,"mu"]
    sigma_not_log = predictions[,"sigma"]

    if (log){
      predictions = predict.lmvar( fit, XX_predict, XX_sigma_predict, log = log)
      mu = predictions[,"mu"]
      sigma = predictions[,"sigma"]
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

    # Exclude rows
    excluded_rows = excluded_rows[foldrows]
    mae = mae[!excluded_rows]
    mse = mse[!excluded_rows]

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
                        message( "  Observations ", dup[1], " and ", dup[2], " (and maybe others) have identical standardized residuals")
                      }
                    })
    }

    return_vec = c( mean(mae), mean(mse), sqrt(mean(mse)))
    if (ks_test){
      return_vec = c( return_vec, ks$statistic, ks$p.value)
    }
    return(return_vec)
  })

  # retrieve results from folds
  mae = numeric()
  mse = numeric()
  ks_distance = numeric()
  ks_p.value = numeric()
  mse_sqrt = numeric()
  for (i in 1:k){
    if (class(cv_results[[i]]) == "list"){
      return(cv_results[[i]])
    }
    else {
      mae[i] = cv_results[[i]][1]
      mse[i] = cv_results[[i]][2]
      mse_sqrt[i] = cv_results[[i]][3]
      if (ks_test){
        ks_distance[i] = cv_results[[i]][4]
        ks_p.value[i] = cv_results[[i]][5]
      }
    }
  }

  outlist = list(MAE = list( mean = mean(mae), sd = stats::sd(mae)),
                 MSE = list( mean = mean(mse), sd = stats::sd(mse)),
                 MSE_sqrt = list( mean = mean(mse_sqrt), sd = stats::sd(mse_sqrt))
                )
  if (ks_test){
    outlist$KS_distance = list( mean = mean(ks_distance), sd = stats::sd(ks_distance))
    outlist$KS_p.value = list( mean = mean(ks_p.value), sd = stats::sd(ks_p.value))
  }
  class(outlist) = "cvlmvar"

  return( outlist)
}
