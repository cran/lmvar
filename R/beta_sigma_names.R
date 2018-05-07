#'
#' @title Unique names for beta_sigma
#'
#' @description Returns adapted names for the coefficients \eqn{\beta_\sigma} to distinguish them from the names
#' of the coefficients
#' \eqn{\beta_\mu}. This is a helper function which is used in situations where it is necessary or convenient
#' for the coefficient names of \eqn{\beta_\sigma}
#' to be different from \eqn{\beta_\mu}.
#'
#' @param beta_mu_names Character vector with the names of the coefficients \eqn{\beta_\mu}
#' @param beta_sigma_names Character vector with the names of the coefficients \eqn{\beta_\sigma}
#' @param ... Additional arguments, not used in the current implementation
#'
#' @return Named character vector with the names of the coefficients \eqn{\beta_\sigma}. The name of a vector element
#' is the original
#' name of the coefficient. The value is the adapted name. The name and the value are equal if no adaptation was needed.
#'
#' @details When the name of at least one coefficient in \eqn{\beta_\sigma} is equal to one of the names of the
#' coefficients in \eqn{\beta_\mu}, the string '_s' is
#' appended to the names of all coefficients in \eqn{\beta_\sigma}.
#' Otherwise, the names of the coefficients in \eqn{\beta_\sigma} are left unchanged.
#'
#' @export
#'
#' @example  R/examples/beta_sigma_names_examples.R
#'

beta_sigma_names <- function( beta_mu_names, beta_sigma_names, ...){

  if (is.null(beta_sigma_names)){
    return(character())
  }
  else
  {
    # Check if at least one name of beta_sigma is identical to beta_mu
    bool = any(is.element( beta_sigma_names, beta_mu_names))

    # Adapt the names of beta_sigma if at least one is identical to beta_mu
    if (bool){
      beta_sigma_names_new = sapply( beta_sigma_names, function(x){
        if (x != "(Intercept_s)"){
          x = paste0( x, "_s")

          # Keep adding "_s" until name is not amongst beta_mu names anymore
          while(x %in% beta_mu_names){
            x = paste0( x, "_s")
          }
        }
        return(x)
      })
    }
    else {
      beta_sigma_names_new = beta_sigma_names
    }

    names(beta_sigma_names_new) = beta_sigma_names

    return(beta_sigma_names_new)
  }
}
