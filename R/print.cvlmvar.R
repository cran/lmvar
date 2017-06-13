#' @title Print method for an object of class 'cvlmvar'
#'
#' @description Print method for an object of class 'cvlmvar'. This object is created by the functions
#' \code{\link{cv.lm}} and \code{\link{cv.lmvar}}.
#'
#' @param x Object of class 'cvlmvar'
#' @param digits Integer, number of significant digits of standard deviations in printed output. Must
#' be larger than zero, or \code{NULL}.
#' @param ... For compatibility with \code{\link[base]{print}} generic.
#'
#' @details If \code{digits = NULL}, printed values are not rounded. Otherwise, all standard deviations are rounded to
#' \code{digits} significant digits. All corresponding mean values are rounded to a precision equal to the maximum
#' precision of the rounded value of the standard deviation.
#'
#' Examples of \code{print.cvlmvar} are provided in the examples of the functions \code{\link{cv.lm}} and
#' \code{\link{cv.lmvar}}.
#'
#' @method print cvlmvar
#'
#' @export
#'

print.cvlmvar <- function( x, digits = NULL, ...){

  format_info <- function( mean, sd, digits){

    order = order_10(sd)

    if (sd < 0.001 | sd >= 100000 | order > digits - 1){
      format = "e"
      digits_sd = digits - 1
      digits_mean = order_10(mean) - order + digits_sd
    }
    else {
      format = "f"
      digits_sd = digits - order - 1
      digits_mean = digits_sd
    }

    return(list( format = format, digits_mean = digits_mean, digits_sd = digits_sd))
  }

  # Check input
  if (!is.null(digits)){
    if (digits <= 0){
      stop("Argument 'digits' must be an integer larger than zero, or NULL")
    }
  }

  ks_test = "KS_distance" %in% names(x)

  if(!is.null(digits)){

    # Print with rounding

    mean_values = c( round_to_sd_accuracy( x$MAE$mean, x$MAE$sd, digits),
                     round_to_sd_accuracy( x$MSE$mean, x$MSE$sd, digits),
                     round_to_sd_accuracy( x$MSE_sqrt$mean, x$MSE_sqrt$sd, digits))

    if (ks_test){
      mean_values = c( mean_values,
                       round_to_sd_accuracy( x$KS_distance$mean, x$KS_distance$sd, digits),
                       round_to_sd_accuracy( x$KS_p.value$mean, x$KS_p.value$sd, digits))
    }

    v1 = mean_values[1]
    v2 = x$MAE$sd
    f = format_info( v1, v2, digits)
    cat("Mean absolute error              : ", formatC( v1, format = f$format, digits = f$digits_mean), "\n")
    cat("Sample standard deviation        : ", formatC( v2, format = f$format, digits = f$digits_sd), "\n\n")

    v1 = mean_values[2]
    v2 = x$MSE$sd
    f = format_info( v1, v2, digits)
    cat("Mean squared error               : ", formatC( v1, format = f$format, digits = f$digits_mean), "\n")
    cat("Sample standard deviation        : ", formatC( v2, format = f$format, digits = f$digits_sd), "\n\n")

    v1 = mean_values[3]
    v2 = x$MSE_sqrt$sd
    f = format_info( v1, v2, digits)
    cat("Square root of mean squared error: ", formatC( v1, format = f$format, digits = f$digits_mean), "\n")
    cat("Sample standard deviation        : ", formatC( v2, format = f$format, digits = f$digits_sd), "\n\n")

    if(ks_test){

      v1 = mean_values[4]
      v2 = x$KS_distance$sd
      f = format_info( v1, v2, digits)
      cat("Kolmogorov-Smirnov distance      : ", formatC( v1, format = f$format, digits = f$digits_mean), "\n")
      cat("Sample standard deviation        : ", formatC( v2, format = f$format, digits = f$digits_sd), "\n\n")

      v1 = mean_values[5]
      v2 = x$KS_p.value$sd
      f = format_info( v1, v2, digits)
      cat("Kolmogorov-Smirnov p-value       : ", formatC( v1, format = f$format, digits = f$digits_mean), "\n")
      cat("Sample standard deviation        : ", formatC( v2, format = f$format, digits = f$digits_sd), "\n\n")
    }
  }
  else {

    # Print without rounding

    cat("Mean absolute error              : ", x$MAE$mean, "\n")
    cat("Sample standard deviation        : ", x$MAE$sd, "\n\n")

    cat("Mean squared error               : ", x$MSE$mean, "\n")
    cat("Sample standard deviation        : ", x$MSE$sd, "\n\n")

    cat("Square root of mean squared error: ", x$MSE_sqrt$mean, "\n")
    cat("Sample standard deviation        : ", x$MSE_sqrt$sd, "\n\n")

    if (ks_test){

      cat("Kolmogorov-Smirnov distance      : ", x$KS_distance$mean, "\n")
      cat("Sample standard deviation        : ", x$KS_distance$sd, "\n\n")

      cat("Kolmogorov-Smirnov p-value       : ", x$KS_p.value$mean, "\n")
      cat("Sample standard deviation        : ", x$KS_p.value$sd, "\n")
    }
  }
}
