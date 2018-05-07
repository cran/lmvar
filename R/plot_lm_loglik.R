#' @title Plot of the log-likelihood surface of a lineair model
#'
#' @description Creates a 3-d plot of the maximum log-likelihood of a lineair model. The maximum log-likelihood is plotted as a function
#' of two elements of the parameter vector \eqn{\beta}. Optionally, the maximum of a quadratic approximation to the log-likelihood
#' surface is plotted.
#'
#' This function is intended for development purposes only.
#'
#' @param y Vector of observations
#' @param X Model matrix
#' @param beta_or Vector of beta values around which the plot is centered. Also the origin of the quadratic appriximation.
#' @param beta_x Component of beta to be plotted at the x-axis. Can be an index of \code{beta_or} or a name in case \code{beta_or}
#' is a named vector
#' @param beta_y Component of beta to be plotted at the y-axis. Can be an index of \code{beta_or} or a name in case \code{beta_or}
#' is a named vector
#' @param add_qa Boolean, specifies whether or not a quadratic approximation of the maximum log-likelihood
#' surface is plotted
#' @param plot_width Single numeric value, half the range of x- and y-values
#' @param plot_points Integer, number of points in the x- and y-dimension at which the maximum-likelihood surface is calculated.
#'
#' @details The function plots the maximum log-likelihood of linear model as a function of two components of the vector \eqn{\beta}.
#' Optionally, it also plots a quadratic approximation to the log-likelihood and plot its maximum.
#'
#' The quadratic approximation is defined as
#'
#' \eqn{log L_q (\beta) = log L(\beta_o) + g (\beta - \beta_o)  + 0.5 (\beta - \beta_o) H (\beta - \beta_o)}
#'
#' where \eqn{log L_q} is the quadratic approximation of the log-likelihood, \eqn{\beta_o} the value of \eqn{\beta} at the origin, \eqn{g}
#' the gradient and \eqn{H} the Hessian of the log-likelihood at \eqn{\beta_o}.
#'
#' For each point \eqn{(\beta_x, \beta_y)}, the other components of the vector \eqn{\beta} are chosen such that the log-likelihood
#' is at its maximum. This maximum is plotted.
#'
#' The same is true for a quadratic approximation of the log-likelihood, except when it has no maximum (i.e. the Hessian
#' \eqn{H} is not negative-definite at \eqn{\beta_o}). In that case, a waning is issued and
#' the other components of \eqn{\beta} are set to the value they have in \code{beta_or}. The resulting quadratic approximation is plotted.
#' This does not affect the plot of the maximum of the true log-likelihood.
#'
#' To create the plot, the package \code{plotly} needs to be installed. A warning is issued if that is not the case.
#'
#' @example R/examples/plot_lm_loglik_examples.R

#'

plot_lm_loglik <- function( y, X, beta_or, beta_x, beta_y, add_qa = FALSE, plot_width = 3, plot_points = 20){

  logl <- function( y, X, beta){

    # Function calculates the log-likelihood of a linear model
    #
    # Input: y: response vector
    #        X: model matrix
    #        beta: vector of betas

    n = length(y)
    res = (y - X %*% beta)
    var = mean(res^2)

    logl = - n * (log(2 * pi * var) + 1) / 2

    return(logl)
  }

  logl_exp <- function( y, X, beta){

    # Function calculates the gradient and Hessian of the profile log-likelihood
    #
    # Input: y: response vector
    #        X: design matrix
    #        beta: vector of betas

    # Set number of observations
    n = length(y)

    # Calculate mu and variance
    res = as.vector((y - X %*% beta))
    var = mean(res^2)

    # Calculate gradient
    grad = as.vector((Matrix::t(X) %*% res) / var)

    # Calculate  Hessian
    H = sapply( grad, function(gx){
      sapply( grad, function(gy){
        gx * gy
      })
    })
    H = (2 / n) * H - (Matrix::t(X) %*% X) / var

    return(list(grad = grad, hess = H))
  }

  logl_qapprox <- function( y, X, beta, beta_or){

    # Function calculates the log-likelihood in a quadratic approximation
    #
    # Input: y: response vector
    #        X: design matrix
    #        beta: vector of betas
    #        beta_or: vector of betas around which you expand


    # Calculate log-likelihood at origin
    logl_or = logl( y, X, beta_or)

    # Calculate approximate log-likelihood
    delta_omega = beta - beta_or
    coeff = logl_exp( y, X, beta_or)

    logl = logl_or + coeff$grad %*% delta_omega + .5 * delta_omega %*% coeff$hess %*% delta_omega

    return(logl)
  }

  # Convert names to indices
  if (inherits( beta_x, "character")){
    beta_x = which(names(beta_or) == beta_x)
  }
  if (inherits( beta_y, "character")){
    beta_y = which(names(beta_or) == beta_y)
  }

  x_plot = seq( from = beta_or[beta_x] - plot_width, to = beta_or[beta_x] + plot_width, length.out = plot_points)
  y_plot = seq( from = beta_or[beta_y] - plot_width, to = beta_or[beta_y] + plot_width, length.out = plot_points)

  # Calculate exact log-likelihood
  z = sapply( x_plot, function(xx){
    sapply( y_plot, function(yy){

      beta_plot = numeric(length(beta_or))
      beta_plot[beta_x] = xx
      beta_plot[beta_y] = yy

      # Calculate betas not in plot
      if (ncol(X) > 2){
        X_r = X[, -c( beta_x, beta_y)]
        beta_r = chol2inv(chol(Matrix::t(X_r) %*% X_r)) %*% Matrix::t(X_r) %*% (y - X[, c( beta_x, beta_y)] %*% c( xx, yy))
        beta_r = as.vector(beta_r)
        beta_plot[-c( beta_x, beta_y)] = beta_r
      }

      logl( y, X, beta_plot)
    })
  })
  colnames(z) = as.character(round( x_plot, 2))
  rownames(z) = as.character(round( y_plot, 2))

  # Calculate log-likelihood in quadratic approximation
  if (add_qa){

    logl_exp_coef = logl_exp( y, X, beta_or)
    grad = logl_exp_coef$grad
    H = logl_exp_coef$hess

    exists_max = TRUE
    if (!matrixcalc::is.negative.definite(as.matrix(H))){
      exists_max = FALSE
      warning("Quadratic approximation has no maximum likelihood")
    }
    else {
      H_r = H[-c( beta_x, beta_y),-c( beta_x, beta_y)]
      H_r = Matrix::Matrix(H_r)
      grad_r = grad[-c( beta_x, beta_y)]
      H_rc = H[-c( beta_x, beta_y), c( beta_x, beta_y)]

      if (ncol(X) > 2){
        H_r_inv = solve(H_r)
      }
    }

    z_qa = sapply( x_plot, function(xx){
      sapply( y_plot, function(yy){
        beta_plot = numeric(length(beta_or))
        beta_plot[beta_x] = xx
        beta_plot[beta_y] = yy

        # Calculate betas not in plot
        if (ncol(X) > 2){
          if(!exists_max){
            beta_plot[-c( beta_x, beta_y)] = beta_or[-c( beta_x, beta_y)]
          }
          else {
            omega = - H_r_inv %*% (grad_r + H_rc %*% c( xx - beta_or[beta_x], yy - beta_or[beta_y]))
            omega = as.vector(omega)
            beta_plot[-c( beta_x, beta_y)] = beta_or[-c( beta_x, beta_y)] + omega
          }
        }
        logl_qapprox( y, X, beta_plot, beta_or)
      })
    })
    colnames(z_qa) = as.character(round( x_plot, 2))
    rownames(z_qa) = as.character(round( y_plot, 2))
  }

  xaxis = list(title = colnames(X)[beta_x])
  yaxis = list(title = colnames(X)[beta_y])
  zaxis = list(title = "logL")

  if (requireNamespace( "plotly", quietly = TRUE)){
    p = plotly::plot_ly(x = x_plot, y = y_plot, z = z)
    p = plotly::add_surface(p)
    if (add_qa){
      p = plotly::add_surface( p, z = z_qa)
    }
    plotly::layout( p, scene = list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis))
  }
  else {
    warning("Package 'plotly' must be loaded to produce the plot")
  }
}
