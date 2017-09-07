#' @title Forward / backward-step model selection for an 'lmvar' object
#'
#' @description Model selection by a forward / backward-stepping algorithm. The algorithm reduces the degrees of freedom of an existing
#' 'lmvar' object. It searches for the subset of degrees of freedom that results in an optimal goodness-of-fit. This is the
#' subset for which a user-specified function reaches its minimum.
#'
#' @param object Object of class 'lmvar'
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
#' \item \code{object} An object of class 'lmvar' which contains the model for which \code{fun} is minimized.
#' \item \code{fun} the minimum value of the user-specified function \code{fun}.
#' }
#'
#' @details
#' \subsection{Description of the algorithm}{
#' The function \code{fwbw.lmvar} selects the subset of all the degrees of freedom present in \code{object} for which the user-specified function
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
#' The start set depends on the argument \code{fw}. If \code{fw = TRUE}, the algorithm starts with only two degrees of freedom: one
#' for the expected values \eqn{\mu} and one for the standard deviations \eqn{\sigma}.
#' These degrees are the intercept terms, if the model in \code{object} contains them.
#' If \code{fw = FALSE} (the default), the algorithm starts with all degrees of freedom present in \code{object}.
#'
#' At a backward step, the model removes \code{removal_percentage} of the degrees of freedom of the current-best subset (with a minimum
#' of 1 degree of freedom). The degrees
#' that are removed are the ones with the largest p-value (p-values can be seen with the function \code{\link{summary.lmvar}}). If the removal
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
#' 'lmvar'. Its return value must be a single number. A smaller (more negative) number must represent a better fit.
#' During the run, a fit to the data is
#' carried out for each new subset of degrees of freedom. The result of the fit is an object of class 'lmvar'. This object is passed on to
#' \code{fun} to evaluate the goodness-of-fit. Typical examples for \code{fun} are \code{\link{AIC.lmvar}} and \code{\link[stats]{BIC}}.
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
#' If \code{object} was created with \code{intercept_mu = TRUE}, the intercept term for the expected values \eqn{\mu} will not be removed by
#' \code{fwbw.lmvar}. Likewise for \code{intercept_sigma}.
#'
#' When a new subset of degrees of freedom is generated by either a backward or a forward step, the response vector in \code{object}
#' is fitted to the new model. The fit is carried out by \code{\link{lmvar}}. The arguments used in the call to
#' \code{lmvar} (other than \code{X_mu}
#' and \code{X_sigma}) are the same as used to create \code{object}.
#' }
#'
#' @seealso
#' \code{\link{fwbw}} for the S3 generic method
#'
#' \code{\link{fwbw.lm}} for the corresponding function for an 'lm' object
#'
#' The degrees of freedom of an 'lmvar' model are given by \code{\link{dfree}}.
#'
#' @export
#'
#' @example R/examples/fwbw.lmvar_examples.R

fwbw.lmvar <- function( object, fun, fw = FALSE, counter = TRUE, removal_percentage = 0.05, control = list(), ...){

  hessian <- function(iterobject){

    # Calculate Hessian

    res = residuals.lmvar( iterobject)
    sigma = fitted.lmvar( iterobject, mu = FALSE)
    lambda = 1 / sigma

    X = iterobject$X_mu * lambda
    H_mu_mu = - Matrix::t(X) %*% X

    lambda_2 = -2 * res * lambda * lambda
    X = iterobject$X_sigma * lambda_2
    H_mu_sigma = Matrix::t(iterobject$X_mu) %*% X

    lambda = res * lambda
    X = iterobject$X_sigma * lambda
    H_sigma_sigma = -2 * Matrix::t(X) %*% X

    outlist = list( residuals = res, sigma = sigma, H_mu_mu = H_mu_mu, H_mu_sigma = H_mu_sigma,
                    H_sigma_sigma = H_sigma_sigma)

    return(outlist)
  }

  sort_logl <- function(iterobject){

    # sort degrees of freedom on estomated effect on likelihood

    # determine which columns can be inserted
    pool_mu = !(colnames(object$X_mu) %in% colnames(iterobject$X_mu))
    pool_mu = colnames(object$X_mu)[pool_mu]
    pool_sigma = !(colnames(object$X_sigma) %in% colnames(iterobject$X_sigma))
    pool_sigma = colnames(object$X_sigma)[pool_sigma]

    if (length(pool_mu) == 0 & length(pool_sigma) == 0){
      return(numeric())
    }
    else {

      # calculate Hessian
      hessian = hessian(iterobject)
      lambda_1 = 1 / hessian$sigma
      lambda_2 = lambda_1 * lambda_1
      lambda_3 = hessian$residuals * lambda_2
      lambda_4 = hessian$residuals * lambda_3
      lambda_5 = hessian$residuals * lambda_1

      # estimate effect on likelihood for degrees of freedom for mu
      effect_mu = sapply( pool_mu, function(col){

        # augment H_mu_mu
        vec = object$X_mu[,col] * lambda_2
        row = as.numeric(Matrix::t(iterobject$X_mu) %*% vec)
        H_mu_mu = cbind( hessian$H_mu_mu, -row)

        vec = object$X_mu[,col] * lambda_1
        row = c( row, as.numeric(vec %*% vec))
        H_mu_mu = rbind( H_mu_mu, -row)

        # augment H_mu_sigma
        vec = -2 * object$X_mu[,col] * lambda_3
        row = vec %*% iterobject$X_sigma
        H_mu_sigma = rbind( hessian$H_mu_sigma, row)

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
          dLogLdB = as.numeric(object$X_mu[,col] %*% lambda_3)

          n = ncol(H_mu_mu)
          return( -0.5 * dLogLdB^2 * H[ n, n])
        }
      })

      # estimate effect on likelihood for degrees of freedom for sigma
      effect_sigma = sapply( pool_sigma, function(col){

        # augment H_mu_sigma
        vec = 2 * object$X_sigma[,col] * lambda_3
        row = Matrix::t(iterobject$X_mu) %*% vec
        H_mu_sigma = cbind( hessian$H_mu_sigma, -row)

        # augment H_sigma_sigma
        vec = 2 * object$X_sigma[,col] * lambda_4
        row = as.numeric(Matrix::t(iterobject$X_sigma) %*% vec)
        H_sigma_sigma = cbind( hessian$H_sigma_sigma, -row)

        vec = object$X_sigma[,col] * lambda_5
        row = c( row, -2 * as.numeric(vec %*% vec))
        H_sigma_sigma = rbind( H_sigma_sigma, row)

        # augment Hessian
        H = rbind( cbind( hessian$H_mu_mu, H_mu_sigma), cbind( Matrix::t(H_mu_sigma), H_sigma_sigma))

        # invert Hessian
        H = tryCatch( solve(H),
                      error = function(e){
                        return("error")
                      })

        if (class(H) == "character"){
          return(NA)
        }
        else{
          dLogLdB = as.numeric(object$X_sigma[,col] %*% (lambda_4 - 1))

          n = ncol(H)
          return( -0.5 * dLogLdB^2 * H[ n, n])
        }
      })

      # set names of estimats for sigma
      names_sigma = beta_sigma_names( colnames(object$X_mu), colnames(object$X_sigma))
      names_sigma = names_sigma[names(effect_sigma)]
      names(effect_sigma) = names_sigma

      # order effects
      if (length(effect_mu) == 0){
        effect_mu = numeric()
      }
      if (length(effect_sigma) == 0){
        effect_sigma = numeric()
      }
      effect_sorted = sort( c( effect_mu, effect_sigma), decreasing = TRUE)

      return(effect_sorted)
    }
  }

  remove_intercept <- function( X, mu){

    # Function removes intercept terms from model matrix, if intercept term is present in model

    if (mu){
      intercept = object$intercept_mu
    }
    else {
      intercept = object$intercept_sigma
    }

    # remove intercept terms from model matrices
    if(intercept){
      cols = colnames(X)
      X = X[,-1]
      if (class(X) == "numeric"){
        X = as.matrix(X)
        colnames(X) = cols[2]
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
  if (class(object) != "lmvar"){
    stop("object must be of class 'lmvar'")
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

  # set whether inserts and removals are done
  if (fw){
    insert = TRUE
    remove = counter
  }
  else {
    remove = TRUE
    insert = counter
  }

  # initialize current model
  if (fw){
    col = colnames(object$X_mu)[1]
    X_mu = Matrix::Matrix(object$X_mu[,1])
    colnames(X_mu) = col

    col = colnames(object$X_sigma)[1]
    X_sigma = Matrix::Matrix(object$X_sigma[,1])
    colnames(X_sigma) = col

    if (object$intercept_mu){
      X_mu = NULL
    }
    if (object$intercept_sigma){
      X_sigma = NULL
    }
    iterobject = lmvar( object$y, X_mu, X_sigma, intercept_mu = object$intercept_mu,
                        intercept_sigma = object$intercept_sigma, sigma_min = object$sigma_min,
                        slvr_options = object$slvr_options, control = object$control)
  }
  else {
    iterobject = object
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
    df = dfree(iterobject)
    df_remove = df - 2

    # Remove degrees of freedom
    success_remove = FALSE
    if (remove & df_remove > 0){

      # iterate over attempts to remove degrees of freedom
      removal_percentage_iter = removal_percentage
      proceed_remove = TRUE
      while(proceed_remove){

        # calculate z-values of coefficients
        z_values = summary.lmvar(iterobject)$coefficients[, "z value"]
        df_mu = ncol(iterobject$X_mu)
        z_in_mu = c( !logical(df_mu), logical(df - df_mu))

        # make sure at least one degree fo freedom remains in both model matrices
        if (object$intercept_sigma){
          i = which(names(z_values) == "(Intercept_s)")
        }
        else {
          i = which.max(abs(z_values[(df_mu + 1):df]))
          i = df_mu + i
        }
        col_sigma = names(z_values)[i]
        z_values = z_values[-i]
        z_in_mu = z_in_mu[-i]

        if (object$intercept_mu){
          i = which(names(z_values) == "(Intercept)")
        }
        else {
#          df_mu = ncol(iterobject$X_mu)
          i = which.max(abs(z_values[1:df_mu]))
        }
        col_mu = names(z_values)[i]
        z_values = z_values[-i]
        z_in_mu = z_in_mu[-i]

        # number of degrees of freedom to remove
        n_remove = max( 1, trunc(removal_percentage_iter * df_remove))

        # calculate columns to keep
        z_sorted = sort( abs(z_values), index.return = TRUE)
        z_in_mu = z_in_mu[z_sorted$ix]

        if (n_remove < df_remove){
          z_sorted = z_sorted$x[(n_remove + 1):df_remove]
          z_in_mu = z_in_mu[(n_remove + 1):df_remove]
        }
        else {
          z_sorted = numeric()
          z_in_mu = logical()
        }

        # remove columns from model matrix for mu
        cols = names(z_sorted)[z_in_mu]
        cols = c( cols, col_mu)
        cols = colnames(iterobject$X_mu) %in% cols
        X_mu = iterobject$X_mu[, cols]
        if (class(X_mu) == "numeric"){
          X_mu = as.matrix(X_mu)
          colnames(X_mu) = colnames(iterobject$X_mu)[cols]
        }

        # remove columns from model matrix for sigma
        cols = names(z_sorted)[!z_in_mu]
        cols = c( cols, col_sigma)
        colnames_adapted = beta_sigma_names( colnames(iterobject$X_mu), colnames(iterobject$X_sigma))
        cols = colnames_adapted %in% cols
        X_sigma = iterobject$X_sigma[, cols]
        if (class(X_sigma) == "numeric"){
          X_sigma = as.matrix(X_sigma)
          colnames(X_sigma) = colnames(iterobject$X_sigma)[cols]
        }

        # remove intercept terms from model matrices
        X_mu = remove_intercept( X_mu, mu = TRUE)
        X_sigma = remove_intercept( X_sigma, mu = FALSE)

        # fit
        fit = lmvar( object$y, X_mu, X_sigma, intercept_mu = object$intercept_mu,
                     intercept_sigma = object$intercept_sigma, sigma_min = object$sigma_min,
                     slvr_options = object$slvr_options, control = object$control)
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
        if (col %in% colnames(object$X_mu)){

          # add column to model matrix for mu
          cols = colnames(object$X_mu) %in% c( colnames(iterobject$X_mu), col)
          X_mu = object$X_mu[, cols]

          X_sigma = iterobject$X_sigma
        }
        else {

          # restore the name of column
          names_sigma = beta_sigma_names( colnames(object$X_mu), colnames(object$X_sigma))
          i = which(names_sigma == col)
          col = names(names_sigma)[i]

          # add column to model matrix for sigma
          cols = colnames(object$X_sigma) %in% c( colnames(iterobject$X_sigma), col)
          X_sigma = object$X_sigma[, cols]

          X_mu = iterobject$X_mu
        }
        # remove intercept terms from model matrices
        X_mu = remove_intercept( X_mu, mu = TRUE)
        X_sigma = remove_intercept( X_sigma, mu = FALSE)

        # fit
        fit = lmvar( object$y, X_mu, X_sigma, intercept_mu = object$intercept_mu,
                     intercept_sigma = object$intercept_sigma, sigma_min = object$sigma_min,
                     slvr_options = object$slvr_options, control = object$control)
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
