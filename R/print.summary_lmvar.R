#' @title Print method for the summary of an 'lmvar' object.
#'
#' @description Print method for an object of the class 'summary_lmvar'. This object is created by \code{\link{summary.lmvar}}.
#'
#' @param x Object of class 'summary_lmvar'
#' @param ... For compatibility with \code{\link[base]{print}} generic.
#'
#' @export
#'
#' @seealso \code{\link{summary.lmvar}} for a summary of the fit present in an object of class 'lmvar'.
#'
#' @method print summary_lmvar
#'

print.summary_lmvar <- function( x, ...){

  cat ("Call: \n ")
  print( x$call)
  cat("\n")

  cat ("Standardized residuals: \n")
  print( round( x$residuals, 4))
  cat("\n")


  # Add aliased coefficients to the table of coefficients
  names_mu = character()
  names_sigma = character()
  n_mu = 0
  n_sigma = 0
  if (x$options$mu){
    names_mu = names(x$aliased_mu)
    n_mu = sum(x$aliased_mu)
  }
  if (x$options$sigma){
    names_sigma = beta_sigma_names( names_mu, names(x$aliased_sigma))
    n_sigma = sum(x$aliased_sigma)
  }
  b_names = c( names_mu, names_sigma)
  df = as.data.frame(x$coefficients)[ b_names,]
  rownames(df) = b_names

  # Number of aliased coefficients
  n = n_mu + n_sigma

  if (x$options$mu & !x$options$sigma){
    cat("Coefficients beta for mu:")
  }
  else if (!x$options$mu & x$options$sigma){
    cat("Coefficients beta for sigma:")
  }
  else {
    cat("Coefficients:")
  }

  if (n==0){
    cat("\n")
  }
  else {
    cat(" (", n, " not defined because of singularities)","\n", sep = "")
  }

  stats::printCoefmat( df, cs.ind = 1:2, tst.ind = 3, has.Pvalue = TRUE, P.values = TRUE)
  cat("\n")

  cat ("Standard deviations: \n")
  print( round( x$sigma, 4))

  if (x$intercept_sigma){
    cat("\n")
    cat("Comparison to model with constant variance (i.e. classical linear model)\n")
    cat("Log likelihood-ratio:", x$logLik_ratio, "\n")
    cat("Additional degrees of freedom:", x$df, "\n")
    cat("p-value for difference in deviance:", signif( x$p_value, 3), "\n")
  }
}
