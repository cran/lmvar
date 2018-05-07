#' @title Forward / backward-step model selection for an 'lmvar' object
#'
#' @description Model selection by a forward / backward-stepping algorithm. The algorithm reduces the degrees of freedom of an existing
#' 'lmvar' object. It searches for the subset of degrees of freedom that results in an optimal goodness-of-fit. This is the
#' subset for which a user-specified function reaches its minimum.
#'
#' @param object Object of class 'lmvar_no_fit' (hence it can also be of class 'lmvar')
#' @param fun User-specified function which measures the goodness-of-fit. See 'Details'.
#' @param fw Boolean, if \code{TRUE} the search will start with a minimum degrees of freedom ('forward search'). If \code{FALSE}
#' the search will start with the full model ('backward search').
#' @param counter Boolean, if \code{TRUE} and \code{fw = TRUE}, the algorithm will carry out backward steps (attempts to
#' remove degrees of freedom) while searching for the optimal subset. If \code{FALSE} and \code{fw = TRUE}, the algorithm will only carry out
#' forward steps (attempts to insert degrees if freedom). The effect of \code{counter} is opposite if \code{fw = FALSE}.
#' @param df_percentage Percentage of degrees of freedom that the algorithm attempts to remove at a backward-step,
#' or insert at a forward-step.
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
#' The function \code{fwbw} selects the subset of all the degrees of freedom present in \code{object} for which the user-specified function
#' \code{fun} is minimized. This function is supposed to be a measure for the goodness-of-fit. Typical examples would be
#' \code{fun=AIC} or \code{fun=BIC}. Another example is where \code{fun} is a measure of the prediction error,
#' determined by cross-validation or otherwise.
#'
#' The function \code{fwbw} is intended for situations in which the degrees of freedom in \code{object} is so large that it is not
#' feasible to go
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
#' At a backward step, the model removes degrees of freedom of the current-best subset. It removes at least 1 degree of freeedom
#' and at most \code{df_percentage} of the degrees in the current-best subset. The degrees
#' that are removed are the ones with the largest p-value (p-values can be seen with the function
#' \code{\link{summary.lmvar}}). If the removal
#' results in a larger value of \code{fun}, the algorithm will try again by halving the degrees of freedom it removes.
#'
#' At a forward step, the algorithm inserts degrees of freedom that are present in
#' \code{object} but left out in the current-best subset. It inserts at least 1 degree of freedom and at most \code{df_percentage} of
#' the  current-best subset. It inserts those
#' degees of freedom which are estimated to
#' increase the likelihood most. If the insertion
#' results in a larger value of \code{fun}, the algorithm will try again by halving the degrees of freedom it inserts.
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
#' and \code{X_sigma}) are the same as used to create \code{object}, except that the control options
#' \code{mu_full_rank} and \code{sigma_full_rank} are both set to TRUE. Setting them to TRUE can be done safely
#' because the model matrices \code{object$X_mu} and \code{object$X_sigma} are guaranteed to be full-rank.
#' }
#'
#' @seealso
#' \code{\link{fwbw}} for the S3 generic method
#'
#' \code{\link{fwbw.lm}} for the corresponding function for an 'lm' object
#'
#' \code{\link{lmvar}} for the constructor of a 'lmvar' object
#'
#' \code{\link{lmvar_no_fit}} for the constructor of a 'lmvar_no_fit' object
#'
#' The number of degrees of freedom is given by \code{\link{dfree}}.
#'
#' @export
#'
#' @example R/examples/fwbw.lmvar_no_fit_examples.R

fwbw.lmvar_no_fit <- function( object, fun, fw = FALSE, counter = TRUE, df_percentage = 0.05, control = list(), ...){

  prep_for_mu <- function(iterobject){

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

  prep_for_sigma <- function(iterobject){

    res = residuals.lmvar( iterobject)
    sigma = fitted.lmvar( iterobject, mu = FALSE)
    lambda = 1 / sigma

    X = iterobject$X_mu * lambda
    M = Matrix::t(X) %*% X
    M = chol2inv(chol(M))
    M = chol(M)

    vec = 2 * res * lambda^2
    M_1 = M %*% Matrix::t(iterobject$X_mu * vec)
    H_base_1 = M_1 %*% iterobject$X_sigma

    vec = res * lambda
    M_2 = iterobject$X_sigma * vec
    H_base_2 = -2 * (Matrix::t(M_2) %*% M_2)
    M_2 = M_2 * (-2 * vec)

    outlist = list( H_base_1 = H_base_1, M_1 = M_1,
                    H_base_2 = H_base_2, M_2 = M_2,
                    lambda = vec^2 - 1, lambda_2 = vec)
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

      # estimate effect on likelihood for degrees of freedom for mu
      effect_mu = numeric()
      if (length(pool_mu) > 0){

        prepared = prep_for_mu(iterobject)
        lambda_1 = 1 / prepared$sigma
        lambda_2 = lambda_1 * lambda_1
        lambda_3 = prepared$residuals * lambda_2
        lambda_4 = prepared$residuals * lambda_3
        lambda_5 = prepared$residuals * lambda_1

        effect_mu = sapply( pool_mu, function(col){

          # augment H_mu_mu
          vec = object$X_mu[,col] * lambda_2
          row = as.numeric(Matrix::crossprod( iterobject$X_mu, vec))
          H_mu_mu = cbind( prepared$H_mu_mu, -row)

          vec = object$X_mu[,col] * lambda_1
          row = c( row, as.numeric(vec %*% vec))
          H_mu_mu = rbind( H_mu_mu, -row)

          # augment H_mu_sigma
          vec = -2 * object$X_mu[,col] * lambda_3
          row = vec %*% iterobject$X_sigma
          H_mu_sigma = rbind( prepared$H_mu_sigma, row)

          # augment Hessian
          H = rbind( cbind( H_mu_mu, H_mu_sigma), cbind( Matrix::t(H_mu_sigma), prepared$H_sigma_sigma))

          # invert Hessian
          H = tryCatch( solve(H),
                        error = function(e){
                          return(NA)
                        })

          if (inherits( H, "logical")){
            return(NA)
          }
          else{

            # Calculate estmated change in log-likelihood
            dLogLdB = as.numeric(object$X_mu[,col] %*% lambda_3)
            n = ncol(H_mu_mu)
            est_logl = -0.5 * dLogLdB^2 * H[ n, n]
            if (est_logl < 0){
              est_logl = -3 * est_logl
            }
            return(est_logl)
          }
        })
      }

      # estimate effect on likelihood for degrees of freedom for sigma
      effect_sigma = numeric()
      if (length(pool_sigma) > 0){

        prepared = prep_for_sigma(iterobject)

        effect_sigma = sapply( pool_sigma, function(col){

          # calculate profile Hessian
          vec = object$X_sigma[,col]

          X = cbind( prepared$H_base_1, prepared$M_1 %*% vec)
          H = Matrix::crossprod(X)

          X = vec %*% prepared$M_2
          x = vec * prepared$lambda_2
          x = -2 * (x %*% x)
          M = rbind( cbind( prepared$H_base_2, Matrix::t(X)),
                     cbind( X, x))
          H = H + M

          # invert Hessian
          H = tryCatch( solve(H),
                        error = function(e){
                          return(NA)
                        })

          if (inherits( H, "logical")){
            return(NA)
          }
          else{

            # Calculate estmated change in log-likelihood
            dLogLdB = as.numeric(vec %*% prepared$lambda)
            n = ncol(H)
            est_logl = -0.5 * H[n,n] * dLogLdB^2
            if (est_logl < 0){
              est_logl = -3 * est_logl
            }
            return(est_logl)
          }
        })

        # set names of estimates for sigma
        names_sigma = beta_sigma_names( colnames(object$X_mu), colnames(object$X_sigma))
        names_sigma = names_sigma[names(effect_sigma)]
        names(effect_sigma) = names_sigma
      }

      # order effects
      if (length(effect_mu) == 0){
        effect_mu = numeric()
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
      if (inherits( X, "numeric")){
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

  slvr_options <- function(){

    slvr_options = object$slvr_options

    if("start" %in% names(slvr_options)){

      start = slvr_options$start
      names(start) = colnames(object$X_sigma)

      if (is.null(X_sigma)){
        start = start[1]
      }
      else {
        if (object$intercept_sigma){
          start = c( start[1], start[colnames(X_sigma)])
          }
        else {
          start = start[colnames(X_sigma)]
        }
      }
      slvr_options$start = start
    }
    return(slvr_options)
  }

  # check inputs
  if (!inherits( object, "lmvar_no_fit")){
    stop("object must be of class 'lmvar' or 'lmvar_no_fit'")
  }
  if (missing(fun)){
    stop("a function must be specified to measure the goodness-of-fit")
  }
  if (df_percentage<=0 | df_percentage>=1){
    stop("the removal percentage must be between 0 and 1")
  }

  # set defaults
  if (!("monitor" %in% names(control))){
    control$monitor = FALSE
  }
  if (!("plot" %in% names(control))){
    control$plot = FALSE
  }
  object_control = object$control
  object_control$mu_full_rank = TRUE
  object_control$sigma_full_rank = TRUE

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

    # set matrices
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
                        slvr_options = slvr_options(), control = object_control)
  }
  else {
    if (inherits( object, "lmvar")){
      iterobject = object
      iterobject$control = object_control
    }
    else {

      X_mu = remove_intercept( object$X_mu, mu = TRUE)
      X_sigma = remove_intercept( object$X_sigma, mu = FALSE)

      iterobject = lmvar( object$y, X_mu, X_sigma, intercept_mu = object$intercept_mu,
                          intercept_sigma = object$intercept_sigma, sigma_min = object$sigma_min,
                          slvr_options = object$slvr_options, control = object_control)
    }
  }

  # set-up list with iteration results
  iterlist = list()
  iterlist[[1]] = list(fun = fun(iterobject))

  # print monitor header
  if (control$monitor){
    cat( "Iteration    attempted removals/insertions    function value\n")
    monitor( 0, 0, iterlist[[1]]$fun)
  }

  # initialize with large value
  n_success_remove = dfree(object)
  n_success_insert = dfree(object)

  # initialize counter
  failed_steps = 0

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
      df_percentage_iter = df_percentage
      proceed_remove = TRUE

      # number of degrees of freedom to remove
      n_remove = min( 2 * n_success_remove, trunc(df_percentage * df_remove))
      n_remove = 2 * n_remove

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
          i = which.max(abs(z_values[1:df_mu]))
        }
        col_mu = names(z_values)[i]
        z_values = z_values[-i]
        z_in_mu = z_in_mu[-i]

        # number of degrees of freedom to remove
        n_remove = max( 1, trunc(n_remove / 2))

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
        if (inherits( X_mu, "numeric")){
          X_mu = as.matrix(X_mu)
          colnames(X_mu) = colnames(iterobject$X_mu)[cols]
        }

        # remove columns from model matrix for sigma
        cols = names(z_sorted)[!z_in_mu]
        cols = c( cols, col_sigma)
        colnames_adapted = beta_sigma_names( colnames(iterobject$X_mu), colnames(iterobject$X_sigma))
        cols = colnames_adapted %in% cols
        X_sigma = iterobject$X_sigma[, cols]
        if (inherits( X_sigma, "numeric")){
          X_sigma = as.matrix(X_sigma)
          colnames(X_sigma) = colnames(iterobject$X_sigma)[cols]
        }

        # remove intercept terms from model matrices
        X_mu = remove_intercept( X_mu, mu = TRUE)
        X_sigma = remove_intercept( X_sigma, mu = FALSE)

        # fit
        fit = lmvar( object$y, X_mu, X_sigma, intercept_mu = object$intercept_mu,
                     intercept_sigma = object$intercept_sigma, sigma_min = object$sigma_min,
                     slvr_options = slvr_options(), control = object_control)
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
          n_success_remove = n_remove
        }
        else {
          if (n_remove == 1){
            proceed_remove = FALSE
            n_success_remove = 0
          }
        }
      }
    }

    # check whether iteration must continue
    if (success_remove){
      failed_steps = 0
    }
    else {
      failed_steps = failed_steps + 1
      if (failed_steps == 2){
        proceed = FALSE
      }
    }

    # Insert degrees of freedom

    success_insert = FALSE
    if (insert & proceed){

      # Get the most promising degrees of freedom
      effect_sorted = sort_logl(iterobject)

      if (length(effect_sorted) > 0){
        proceed_insert = TRUE
      }
      else {
        proceed_insert = FALSE
      }

      # number of degrees of freedom to insert
      n_insert = min( 2 * n_success_insert, trunc(df_percentage * length(effect_sorted)))
      n_insert = 2 * n_insert

      while (proceed_insert){

        # number of degrees of freedom to insert
        n_insert = max( 1, trunc(n_insert / 2))

        # add new degrees of freedom to model matrices
        col = names(effect_sorted[1:n_insert])
        cols = colnames(object$X_mu) %in% c( colnames(iterobject$X_mu), col)
        X_mu = object$X_mu[, cols]
        if (inherits( X_mu, "numeric")){
          X_mu = as.matrix(X_mu)
          colnames(X_mu) = colnames(object$X_mu)[cols]
        }

        names_sigma = beta_sigma_names( colnames(object$X_mu), colnames(object$X_sigma))
        col = names(names_sigma[names_sigma %in% col])
        cols = colnames(object$X_sigma) %in% c( colnames(iterobject$X_sigma), col)
        X_sigma = object$X_sigma[, cols]
        if (inherits( X_sigma, "numeric")){
          X_sigma = as.matrix(X_sigma)
          colnames(X_sigma) = colnames(object$X_sigma)[cols]
        }

        # remove intercept terms from model matrices
        X_mu = remove_intercept( X_mu, mu = TRUE)
        X_sigma = remove_intercept( X_sigma, mu = FALSE)

        # fit
        fit = lmvar( object$y, X_mu, X_sigma, intercept_mu = object$intercept_mu,
                     intercept_sigma = object$intercept_sigma, sigma_min = object$sigma_min,
                     slvr_options = slvr_options(), control = object_control)
        fun_iter = fun(fit)

        # set iteration counter
        iter = length(iterlist)

        # monitor
        if (control$monitor){
          monitor( iter, n_insert, fun_iter, insert = TRUE)
        }

        if (fun_iter < iterlist[[iter]]$fun){
          success_insert = TRUE
          proceed_insert = FALSE
          iterlist[[iter + 1]] = list(fun = fun_iter)
          iterobject = fit
          n_success_insert = n_insert
        }
        else {
          if (n_insert > 1){

            # check if insertion of degrees of freedom increases log-likelihood
            if (logLik.lmvar(fit) < logLik.lmvar(iterobject)){

              # enforce that only 1 degree of freedom is inserted in next attempt
              n_insert = 1
            }
          }
          else {
          proceed_insert = FALSE
          n_success_insert = 0
          }
        }
      }
    }

    # check whether iteration must continue
    if (success_insert){
      failed_steps = 0
    }
    else {
      failed_steps = failed_steps + 1
      if (failed_steps == 2){
        proceed = FALSE
      }
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
