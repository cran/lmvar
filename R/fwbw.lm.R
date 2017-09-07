#' @title Forward / backward-step model selection for an object of class 'lm'
#'
#' @description Model selection by a forward / backward-stepping algorithm. The algorithm reduces the degrees of freedom of an existing
#' 'lm' object. It searches for the subset of degrees of freedom that results in an optimal goodness-of-fit. The optimal subset is the
#' subset for which a user-specified function reaches its minimum.
#'
#' @param object Object of class 'lm'
#' @param fun User-specified function which measures the goodness-of-fit. See 'Details'.
#' @param fw Boolean, if \code{TRUE} the search will start with a minimum degrees of freedom ('forward search'). If \code{FALSE}
#' the search will start with the full model ('backward search').
#' @param counter Boolean, if \code{TRUE} and \code{fw = TRUE}, the algorithm will carry out backward steps (attempts to
#' remove degrees of freedom) while searching for the optimal subset. If \code{FALSE} and \code{fw = TRUE}, the algorithm will only carry out
#' forward steps (attempts to insert degrees if freedom). The effect of \code{counter} is opposite if \code{fw = FALSE}.
#' @param removal_percentage Percentage of degrees of freedom that the algorithm attempts to remove at a backward-step.
#' Must be a number between 0 and 1.
#' @param control List of control options. The following options can be set
#' \itemize{
#' \item \code{monitor} Boolean, if \code{TRUE} information about the attempted removals and insertions will be printed during the run.
#' Default is \code{FALSE}.
#' \item \code{plot} Boolean, if \code{TRUE} a plot will be shown at the end of the run. It shows how the value of \code{fun}
#' decreases during the run. Default is \code{FALSE}.
#' }
#' @param ... for compatibility with \code{\link{fwbw}} generic
#'
#' @return A list with the following members.
#' \itemize{
#' \item \code{object} An object of class 'lm' which contains the model for which \code{fun} is minimized.
#' \item \code{fun} The minimum value of the user-specified function \code{fun}.
#' }
#'
#' @details
#' \subsection{Description of the algorithm}{
#' The function \code{fwbw.lm} selects the subset of all the degrees of freedom present in \code{object} for which the user-specified function
#' \code{fun} is minimized. This function is supposed to be a measure for the foodness-of-fit. Typical examples would be
#' \code{fun=AIC} or \code{fun=BIC}. The function \code{fun} can also be a measure of the prediction error,
#' determined by cross-validation.
#'
#' This function is intended for situations in which the degrees of freedom in \code{object} is so large that it is not feasible to go
#' through all possible subsets systematically to find the smallest value of \code{fun}. Instead, the algorithm generates subsets by removing
#' degrees of freedom from the current-best subset (a 'backward' step) and reinserting degrees of freedom that were previously removed
#' (a 'forward' step). Whenever a backward or forward step results in a subset for which \code{fun} is smaller than for the current-best
#' subset, the new subset becomes current-best.
#'
#' The start set depends on the argument \code{fw}. If \code{fw = TRUE}, the algorithm starts with only one degree of freedom
#' for the expected values \eqn{\mu}. This degree is the intercept term, if the model in \code{object} contains an intercept term.
#' If \code{fw = FALSE} (the default), the algorithm starts with all degrees of freedom present in \code{object}.
#'
#' At a backward step, the model removes \code{removal_percentage} of the degrees of freedom of the current-best subset (with a minimum
#' of 1 degree of freedom). The degrees
#' that are removed are the ones with the largest p-value (p-values can be seen with the function \code{\link[stats]{summary.lm}}). If the removal
#' results in a larger value of \code{fun}, the algorithm will try again by halving the degrees of freedom it removes.
#'
#' At a forward step, the algorithm goes through all degrees of freedom that have
#' been removed from \code{object} so far. For each degree of freedom, the algorithm estimates the increase of the log-likelihood if the
#' degree were to be inserted. It inserts the one which is estimated to increase the log-likelihood the most, followed by the one with
#' the second-largest estimated increase, etc. After each insertion, the value of \code{fun} is calculated. As long as \code{fun} decreases,
#' the insertion process continues.
#'
#' If \code{counter = FALSE}, the algorithm is 'greedy': it will only carry out forward-steps in case \code{fw = TRUE} or backward-steps
#' in case \code{fw = FALSE}.
#'
#' The algorithm stops if neither the backward nor the forward step resulted in a lower value of \code{fun}. It returns the current-best model
#' and the minimum value of \code{fun}.
#' }
#'
#' \subsection{The user-defined function}{
#' The function \code{fun} must be a function which is a measure for the goodness-of-fit. It must take one argument: an object of class
#' 'lm'. Its return value must be a single number. A smaller number (more negative) must represent a better fit.
#' During the run, a fit to the data is
#' carried out for each new subset of degrees of freedom. The result of the fit is an object of class 'lm'. This object is passed on to
#' \code{fun} to evaluate the goodness-of-fit. Typical examples for \code{fun} are \code{\link[stats]{AIC}} and
#' \code{\link[stats]{BIC}}.
#'}
#'
#' \subsection{Monitor information}{
#' When the \code{control}-option \code{monitor} is equal to \code{TRUE}, information is displayed about the progress of the run.
#' The following information is displayed:
#' \itemize{
#' \item \code{Iteration} A counter which first value is always \code{0}, followed by \code{1}. From then on, the counter is increased
#' whenever the addition or removal of degrees of freedom results in a smaller function value than the smallest so far.
#' \item \code{attempted removals/insertions} The number of degrees of freedoms that one attempts to remove or insert
#' \item \code{function value} The value of the user-specified function \code{fun} after the removal or insertion of the degrees of freedom
#' \item The last column shows the word \code{insert} when the attempt regards the insertion of degrees of freedom. When nothing is shown,
#' the algorithm attempted to remove degrees of freedom.
#' }
#' }
#'
#' \subsection{Other}{
#' If the model matrix present in \code{object} conatains a column with the name "(Intercept)", the intercept term for
#' the expected values \eqn{\mu} will not be removed by
#' \code{fwbw.lm}.
#'
#' When a new subset of degrees of freedom is generated by either a backward or a forward step, the response vector in \code{object}
#' is fitted to the new model. The fit is carried out by \code{\link[stats]{lm}}.
#' }
#'
#' @seealso
#' \code{\link{fwbw}} for the generic method
#'
#' \code{\link{fwbw.lmvar}} for the corresponding function for an 'lmvar' object
#'
#' @export
#'
#' @example R/examples/fwbw.lm_examples.R
#'

fwbw.lm <- function( object, fun, fw = FALSE, counter = TRUE, removal_percentage = 0.05, control = list(), ...){

  hessian <- function(iterobject){

    # Calculate Hessian

    res = stats::residuals( iterobject)
    sigma = sqrt(sum(res^2) / stats::nobs(iterobject))
    lambda = 1 / sigma

    X = iterobject$x * lambda
    H = - Matrix::t(X) %*% X

    outlist = list( residuals = res, sigma = sigma, H_mu_mu = H, H_sigma_sigma = -2 * stats::nobs(iterobject) * lambda^2)

    return(outlist)
  }

  sort_logl <- function(iterobject){

    # sort degrees of freedom on estomated effect on likelihood

    # determine which columns can be inserted
    pool_mu = !(colnames(object$x) %in% itercols)
    pool_mu = colnames(object$x)[pool_mu]

    if (length(pool_mu) == 0){
      return(numeric())
    }
    else {

      # calculate Hessian
      hessian = hessian(iterobject)
      lambda_1 = 1 / hessian$sigma
      lambda_2 = lambda_1 * lambda_1
      lambda_3 = hessian$residuals * lambda_2
      lambda_4 = lambda_1 * lambda_3

      # estimate effect on likelihood for degrees of freedom for mu
      effect_mu = sapply( pool_mu, function(col){

        # augment H_mu_mu
        vec = object$x[,col] * lambda_2
        row = as.numeric(Matrix::t(iterobject$x) %*% vec)
        H_mu_mu = cbind( hessian$H_mu_mu, -row)

        vec = object$x[,col] * lambda_1
        row = c( row, as.numeric(vec %*% vec))
        H_mu_mu = rbind( H_mu_mu, -row)

        # augment H_mu_sigma
        n = nrow(H_mu_mu)
        vec = numeric(n)
        vec[n] = -2 * object$x[,col] %*% lambda_4
        H_mu_sigma = matrix(vec)

        # augment Hessian
        H = rbind( cbind( H_mu_mu, H_mu_sigma), cbind( Matrix::t(H_mu_sigma), hessian$H_sigma_sigma))


        # invert Hessian
        H = tryCatch( solve(H),
                      error = function(e){
                        return("error")
                      })

        if (class(H) == "character"){
          return(NA)
        }
        else{
          dLogLdB = as.numeric(object$x[,col] %*% lambda_3)

          n = ncol(H_mu_mu)
          return( -0.5 * dLogLdB^2 * H[ n, n])
        }
      })

      # order effects
      if (length(effect_mu) == 0){
        effect_mu = numeric()
      }
      effect_sorted = sort( effect_mu, decreasing = TRUE)

      return(effect_sorted)
    }
  }

  remove_intercept <- function(X){

    # Function removes intercept terms from model matrix, if intercept term is present in model

    intercept = "(Intercept)" %in% colnames(X)

    # remove intercept terms from model matrices
    if(intercept){
      cols = colnames(X)
      i = which(colnames(X) == "(Intercept)")
      X = X[,-i]
      if (class(X) == "numeric"){
        X = as.matrix(X)
        colnames(X) = cols[-i]
      }
    }
    return (X)
  }

  monitor <- function( iter, n, fun, insert = FALSE){
    s = "\n"
    if (insert){
      s = paste(" insert", s)
    }
    cat( format( iter, width = 9), format( n, width = 32), format( fun, width = 17), s)
  }

  # check inputs
  if (class(object) != "lm"){
    stop("object must be of class 'lm'")
  }
  if (!("x" %in% names(object)) | !("y" %in% names(object))){
    stop("object must contain response vector and model matrix, please run 'lm' with 'x = TRUE' and 'y = TRUE'")
  }
  if (missing(fun)){
    stop("a function must be specified to measure the goodness-of-fit")
  }
  if (removal_percentage<=0 | removal_percentage>=1){
    stop("the removal percentage must be between 0 and 1")
  }

  # set defaults
  if (!("monitor" %in% names(control))){
    control$monitor = FALSE
  }
  if (!("plot" %in% names(control))){
    control$plot = FALSE
  }

  # set whether inserts and removals are done
  if (fw){
    insert = TRUE
    remove = counter
  }
  else {
    remove = TRUE
    insert = counter
  }

  # check if intercept column is present
  intercept_mu = "(Intercept)" %in% colnames(object$x)

  # initialize current model
  if (fw){
    if (intercept_mu){
      iterobject = stats::lm( object$y ~ 1, x = TRUE, y = TRUE)
      itercols = "(Intercept)"          # keep track of column names because 'lm' modifies them
    }
    else{
      cols = colnames(object$x)[1]
      X_mu = as.matrix(object$x[,1])
      colnames(X_mu) = cols
      iterobject = stats::lm( object$y ~ . + 0, as.data.frame(X_mu), x = TRUE, y = TRUE)
      itercols = cols
    }
  }
  else {
    iterobject = object
    itercols = colnames(object$x)
  }

  # set-up list with iteration results
  iterlist = list()
  iterlist[[1]] = list(fun = fun(iterobject))

  # print monitor header
  if (control$monitor){
    cat( "Iteration    attempted removals/insertions    function value\n")
    monitor( 0, 0, iterlist[[1]]$fun)
  }

  # iterate over backward-forward steps
  proceed = TRUE
  while(proceed){

    # degrees of freedom that can be removed
    df = ncol(iterobject$x) + 1
    df_remove = df - 2

    # Remove degrees of freedom
    success_remove = FALSE
    if (remove & df_remove > 0){

      # iterate over attempts to remove degrees of freedom
      removal_percentage_iter = removal_percentage
      proceed_remove = TRUE
      while(proceed_remove){

        # calculate z-values of coefficients
        z_values = stats::summary.lm(iterobject)$coefficients[, "t value"]

        # make sure at least one degree fo freedom remains
        if (intercept_mu){
          i = which(names(z_values) == "(Intercept)")
        }
        else {
          i = which.max(abs(z_values))
        }
        col_mu = names(z_values)[i]
        z_values = z_values[-i]

        # number of degrees of freedom to remove
        n_remove = max( 1, trunc(removal_percentage_iter * df_remove))

        # calculate columns to keep
        z_sorted = sort( abs(z_values), index.return = TRUE)

        if (n_remove < df_remove){
          z_sorted = z_sorted$x[(n_remove + 1):df_remove]
        }
        else {
          z_sorted = numeric()
        }

        # remove columns from model matrix for mu
        cols = names(z_sorted)
        cols = c( cols, col_mu)
        cols = colnames(iterobject$x) %in% cols
        X_mu = iterobject$x[, cols]
        if (class(X_mu) == "numeric"){
          X_mu = as.matrix(X_mu)
        }
        colnames(X_mu) = itercols[cols]

        # remove intercept terms from model matrices
        X_mu = remove_intercept(X_mu)

        # fit
        if (intercept_mu){
          if (ncol(X_mu) > 0){
            fit = stats::lm(object$y ~ ., as.data.frame(as.matrix(X_mu)), x = TRUE, y = TRUE)
            cols_fit = c( "(Intercept)", colnames(X_mu))
          }
          else {
            fit = stats::lm(object$y ~ 1, x = TRUE, y = TRUE)
            cols_fit = "(Intercept)"
          }
        }
        else {
          fit = stats::lm(object$y ~ . + 0, as.data.frame(as.matrix(X_mu)), x = TRUE, y = TRUE)
          cols_fit = colnames(X_mu)
        }
        fun_iter = fun(fit)

        # set iteration counter
        iter = length(iterlist)

        # monitor
        if (control$monitor){
          monitor( iter, n_remove, fun_iter)
        }

        if (fun_iter <= iterlist[[iter]]$fun){
          success_remove = TRUE
          proceed_remove = FALSE
          iterlist[[iter + 1]] = list(fun = fun_iter)
          iterobject = fit
          itercols = cols_fit
        }
        else {
          if (n_remove > 1){
            removal_percentage_iter = removal_percentage_iter / 2
          }
          else {
            proceed_remove = FALSE
          }
        }
      }
    }
    # Insert degrees of freedom

    success_insert = FALSE
    if (insert){

      # Get the most promising degrees of freedom
      effect_sorted = sort_logl(iterobject)

      if (length(effect_sorted) > 0){
        proceed_insert = TRUE
      }
      else {
        proceed_insert = FALSE
      }
      while (proceed_insert){

        # add new degree of freedom to model matrix
        col = names(effect_sorted[1])

        # add column to model matrix for mu
        cols = colnames(object$x) %in% c( itercols, col)
        X_mu = object$x[, cols]

        # remove intercept terms from model matrices
        X_mu = remove_intercept(X_mu)

        # fit
        if (intercept_mu){
          fit = stats::lm(object$y ~ ., as.data.frame(as.matrix(X_mu)), x = TRUE, y = TRUE)
          cols_fit = c( "(Intercept)", colnames(X_mu))
        }
        else {
          fit = stats::lm(object$y ~ . + 0, as.data.frame(as.matrix(X_mu)), x = TRUE, y = TRUE)
          cols_fit = colnames(X_mu)
        }
        fun_iter = fun(fit)

        # set iteration counter
        iter = length(iterlist)

        # monitor
        if (control$monitor){
          monitor( iter, 1, fun_iter, insert = TRUE)
        }

        if (fun_iter < iterlist[[iter]]$fun){
          success_insert = TRUE
          iterlist[[iter + 1]] = list(fun = fun_iter)
          iterobject = fit
          itercols = cols_fit
        }
        else {
          proceed_insert = FALSE
        }

        # proceed with next degree of freedom
        effect_sorted = effect_sorted[-1]
        if (length(effect_sorted) == 0){
          proceed_insert = FALSE
        }
      }
    }

    if (!success_remove & !success_insert){
      proceed = FALSE
    }

  }
  iter = length(iterlist)

  # plot sequence of fun-values
  if (control$plot){
    funs = sapply( iterlist, function(iter){
      return(iter$fun)
    })
    graphics::plot( 0:(iter - 1), funs, xlab = "iteration", ylab = "function value")
  }

  outlist = list( object = iterobject, fun = iterlist[[iter]]$fun)

  return(outlist)
}
