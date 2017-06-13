#' @title Cross-validation for an object of class 'lm'
#'
#' @description k-fold cross-validation for an object of class 'lm'
#'
#' @param object Object of class 'lm'
#' @param k Integer, number of folds
#' @param ks_test Boolean, if \code{TRUE}, a Kolmogorov-Smirnov test is carried out. See details.
#' @param log Boolean, specifies whether \code{object} contains a fit to the response vector \eqn{Y} or its logarithm \eqn{\log Y}
#' @param seed Integer, seed for the random number generator. The seed is not set when \code{seed} equals \code{NULL}.
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
#' @details \code{object} must contain the list-members \code{x} and \code{y}. I.e., it must be created by running
#' \code{\link[stats]{lm}} with the options \code{x = TRUE} and \code{y = TRUE}.
#'
#' When \code{ks_test = TRUE}, a Kolmogorov-Smirnov (KS) test is carried out. The test checks whether the
#' standardized residuals \eqn{(y - \mu) / \sigma} in a fold are distributed as a standard normal distribution. The
#' KS-distance and the corresponding p-value are calculated for each fold. The test uses the function
#' \code{\link[stats]{ks.test}}.
#'
#' @seealso \code{\link{cv.lmvar}} is the equivalent function for an object of class 'lmvar'. It is supplied in
#' case one wants to compare an  'lmvar' fit with an 'lm' fit.
#'
#' \code{\link{print.cvlmvar}} provides a print-method for an object of class 'cvlmvar'.
#'
#' @export
#'
#' @example R/examples/cv.lm_examples.R
#'
cv.lm <- function( object, k = 10, ks_test = FALSE, log = FALSE, seed = NULL,  ...){

  # Check input
  if (!any(class(object) == "lm")){
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

  # Retrieve info from object
  y = object$y
  X = object$x

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

    # create model matrix
    XX = make_matrix_full_rank(X[!foldrows,])
    XX_predict = X[ foldrows, colnames(XX)]

    # perform fit on rows not in fold
    fit = stats::lm( y[!foldrows] ~ . - 1, as.data.frame(XX))

    # predict values for rows in fold
    mu_not_log = stats::predict( fit, as.data.frame(XX_predict))
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
