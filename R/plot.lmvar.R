#' @title Plot diagnostics for an 'lmvar' object
#'
#' @description This function produces 5 plots which should help to judge the goodness of an 'lmvar' fit.
#'
#' @param x Object of class 'lmvar'
#' @param which Integer vector slecting which of the 5 plots is produced
#' @param id.n Integer, the number of 'extreme' observations that are labelled in the plots
#' @param cex.id Numeric, scale-factor for the size of the observation labels in the plots
#' @param show Boolean, if TRUE the number of the plot is shown in the plot-title and the name
#' of \code{x} is shown in the label of the x-axis.
#' @param ... for compatibility with \code{\link[graphics]{plot}} generic
#'
#' @details The plots are intended to be a quick and easy way to get an impression of the
#' goodness-of-fit. The function is intended for an interactive R-session and users must hit <enter>
#' before each plot is deplayed. The following plots can be produced.
#' \enumerate{
#' \item A plot of the residuals \eqn{y - \mu} versus the fitted values \eqn{\mu}.
#' \item A QQ-plot, showing the z-score \eqn{(y - \mu) / \sigma} resulting from the fit versus
#' the z-score calculated from the sample quantiles. The sample quantiles are calculated as
#' \code{\link[stats]{ppoints}}(\eqn{n}) with \eqn{n} the number of obervations in \code{x}.
#' \item A histogram of the distribution of the quantiles of the response values. The quantiles
#' are calculated under the assumption that the response values are normally distributed with expected values
#' \eqn{\mu} and standard deviations \eqn{\sigma}.
#' \item A plot of the z-scores versus the fitted values.
#' \item A scale-location plot showing the square root of the absolute z-scores versus the fitted values.
#' }
#' If relevant, plots show the average y-value as a red line. This line is created by the function
#' \code{\link[graphics]{panel.smooth}}.
#' If relevant, plots show the expected average y-value as a dotted gray line.
#'
#' To suppress labelling of observations in the plots, set \code{id.n} to zero or a negative
#' value. If \code{id.n}  is set to a value equal to or larger than the number of observations
#' in \code{x}, all points in the plots are labelled.
#'
#' @return There is no return value. The function only shows plots in the graphics output device.
#'
#' @export
#'
#' @example R/examples/plot.lmvar_examples.R
#'

plot.lmvar <- function( x, which = c(1:3, 5), id.n = 3, cex.id = 0.75, show = TRUE, ...){

  object = x
  object_name = deparse(substitute(x))

  mu = fitted.lmvar( object, sigma = FALSE)
  res = residuals.lmvar(object)
  sigma = fitted.lmvar( object, mu = FALSE)
  z = res / sigma
  n = nobs.lmvar_no_fit(object)

  # ask for new page between plots
  oask <- grDevices::devAskNewPage(TRUE)
  on.exit(grDevices::devAskNewPage(oask))

  plot_number = 1
  if (plot_number %in% which){

    # Plot residuals versus fitted values
    xlab = "Fitted values"
    main = "Residuals vs fitted values"
    if (show){
      xlab = paste0( xlab, " (", object_name, ")")
      main = paste0( main, " (", plot_number, ")")
    }
    graphics::plot( mu, res, xlab = xlab, ylab = "Residuals", main = main)
    graphics::abline( 0, 0, lty = "dotted", col = "gray")
    graphics::panel.smooth(mu, res)

    # Add observation labels
    y_ordered = sort( abs(res), decreasing = TRUE, index.return = TRUE)
    i = min( id.n, length(y_ordered$ix))
    if (i > 0){
      i = y_ordered$ix[1:i]
      graphics::text( mu[i], res[i], labels = i, pos = 2, cex = cex.id)
    }
  }
  plot_number = 2
  if (plot_number %in% which){

    # QQ-plot
    y = sort(z)
    x = stats::ppoints(n)
    x = stats::qnorm(x)

    xlab = "z-scores from sample quantile"
    main = "Normal Q-Q plot"
    if (show){
      xlab = paste0( xlab, " (", object_name, ")")
      main = paste0( main, " (", plot_number, ")")
    }
    graphics::plot( x, y, main = main, xlab = xlab, ylab = "z-scores from fit")
    graphics::abline( 0, 1, lty = "dotted", col = "gray")

    # Add observation labels
    y_ordered = sort( abs(z), decreasing = TRUE, index.return = TRUE)
    i = min( id.n, length(y_ordered$ix))
    if (i > 0){
      obs = y_ordered$ix[1:i]
      i = which(z[obs] < 0)
      if (length(i) > 0){
        j = 1:length(i)
        graphics::text( x[j], y[j], labels = obs[i], pos = 4, cex = cex.id)
      }
      i = which(z[obs] >= 0)
      if (length(i) > 0){
        j = n:(n - length(i) + 1)
        graphics::text( x[j], y[j], labels = obs[i], pos = 2, cex = cex.id)
      }
    }
  }
  plot_number = 3
  if (plot_number %in% which){

    # Plot distribution of quantiles
    y = stats::pnorm( object$y, mean = mu, sd = sigma)

    xlab = "Quantile of response value"
    main = "Distribution of quantiles"
    if (show){
      xlab = paste0( xlab, " (", object_name, ")")
      main = paste0( main, " (", plot_number, ")")
    }
    graphics::hist( y, main = main, xlab = xlab, ylab = "Probability density", freq = FALSE)
    graphics::abline( 1, 0, lty = "dotted", col = "gray")
  }
  plot_number = 4
  if (plot_number %in% which){

    # Plot z-scores versus fitted values
    xlab = "Fitted values"
    main = "z-scores versus fitted values"
    if (show){
      xlab = paste0( xlab, " (", object_name, ")")
      main = paste0( main, " (", plot_number, ")")
    }
    graphics::plot( mu, z, xlab = xlab, ylab = "z-scores", main = main)
    graphics::abline( 0, 0 , lty = "dotted", col = "gray")
    graphics::panel.smooth(mu, z)

    # Add observation labels
    y_ordered = sort( abs(z), decreasing = TRUE, index.return = TRUE)
    i = min( id.n, length(y_ordered$ix))
    if (i > 0){
      i = y_ordered$ix[1:i]
      graphics::text( mu[i], z[i], labels = i, pos = 2, cex = cex.id)
    }
  }
  plot_number = 5
  if (plot_number %in% which){

    # Plot square root of absolute z-scores versus fitted values
    y = sqrt(abs(z))
    xlab = "Fitted values"
    main = "Scale - location"
    if (show){
      xlab = paste0( xlab, " (", object_name, ")")
      main = paste0( main, " (", plot_number, ")")
    }
    graphics::plot( mu, y, xlab = xlab, ylab = "Square root of absolute z-score",
                    main = main)
    graphics::abline( 2^{1/4} * gamma(3/4) / sqrt(pi), 0 , lty = "dotted", col = "gray")
    graphics::panel.smooth(mu, y)
  }
}
