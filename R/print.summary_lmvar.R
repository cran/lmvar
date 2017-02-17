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
  b_names = c( names(x$aliased_mu), beta_sigma_names( names(x$aliased_mu), names(x$aliased_sigma)))
  df = x$coefficients[ b_names,]
  rownames(df) = b_names

  # Number of aliased coefficients
  n = sum(x$aliased_mu) + sum(x$aliased_sigma)
  if (n==0){
    cat("Coefficients:\n")
  }
  else {
    cat("Coefficients: (", n, " not defined because of singularities)","\n", sep = "")
  }

  stats::printCoefmat( df, cs.ind = 1:2, tst.ind = 3, has.Pvalue = TRUE, P.values = TRUE)
  cat("\n")

  cat ("Standard deviations: \n")
  print( round( x$sigma, 4))
  cat("\n")

  cat("Comparison to model with constant variance (i.e. classical linear model)\n")
  cat("Log likelihood-ratio:", x$logLik_ratio, "\n")
  cat("Additional degrees of freedom:", x$df, "\n")
  cat("p-value for difference in deviance:", signif( x$p_value, 3), "\n")
}
