#' @title  Forward / backward-step model selection
#'
#' @description Model selection by a forward / backward-stepping algorithm. The algorithm reduces the degrees of freedom of an existing
#' object containing a model fit. It searches for the subset of degrees of freedom that results in an optimal goodness-of-fit. This is the
#' subset for which a user-specified function reaches its minimum. The search is carried out by
#' alternately attempting to remove and insert degrees of freedom.
#'
#' @param object Object containing a fit to a specific model
#' @param fun User-specified function which measures the goodness-of-fit.
#' @param ... Further arguments for specific methods
#'
#' @return A list with the following members.
#' \itemize{
#' \item \code{object} An object which contains the model for which \code{fun} is minimized.
#' \item \code{fun} the minimum value of the user-specified function \code{fun}.
#' }
#'
#' @seealso \code{\link{fwbw.lm}} and \code{\link{fwbw.lmvar_no_fit}}
#'
#' @export
#'

fwbw <- function( object, fun, ...){
  UseMethod( "fwbw", object)
}
