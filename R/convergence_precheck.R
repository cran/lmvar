#' @title Pre-check model matrices for convergence issues
#'
#' @description The model matrices \eqn{X_\mu} and \eqn{X_\sigma} are checked to see if problems
#' with the convergence of the fit can be anticipated. If so, it is determined which columns
#' must be removed from \eqn{X_\sigma} to attempt to avoid convergence issues.
#'
#' @param y Numeric, response vector y
#' @param X_mu Model matrix for the expected values
#' @param X_sigma Model matrix for the standard deviations. This must be a full-rank matrix.
#'
#' @details A matrix can be of class 'matrix',
#' 'Matrix' or 'numeric' (in case it is a matrix of one column only).
#'
#' An intercept term must be included in the model matrices if the model is such.
#'
#' @return A list with the following members:
#' \itemize{
#' \item \code{column_numbers} The numbers of the columns of \code{X_sigma} that can be kept
#' \item \code{column_names} The names of the columns of \code{X_sigma} that can be kept
#' }
#' Numbers and names refer to the same columns. They are supplied both for convenience.
#'

convergence_precheck <- function( y, X_mu, X_sigma){

  # Store column names and replace them by column numbers
  column_names = colnames(X_sigma)
  colnames(X_sigma) = as.character(seq.int( 1, ncol(X_sigma)))

  # Set tolerance
  tol = max(dim(X_sigma)) * .Machine$double.eps

  # set number of observations
  n_obs = length(y)

  ## Try to reduce collinearity by a QR-decomposition

  # Loop over successive reductions in degrees of freedom
  continue = TRUE
  while (continue) {

    # QR decomposition of X_sigma
    qr_sigma = Matrix::qr(X_sigma)
    Q = Matrix::qr.Q(qr_sigma)

    # Calculate number of zeros in column of Q
    zeros = apply(Q, 2, function(x) {
      return(sum(abs(x) < tol))
    })
    zeros = sort(zeros, decreasing = TRUE, index.return = TRUE)

    # Loop over all possibilities to reduce rank of X_sigma by deleting rows
    i = 1
    success = FALSE
    at_end = FALSE
    while (!success & !at_end & zeros$x[i] > 0.5 * n_obs) {    # Do not attempt to remove more than half of the observations
      # Calculate rows of Q in which zeros occur
      j = which(abs(Q[, zeros$ix[i]]) < tol)

      # Check that X_mu_1 * beta = y_1 has a solution for beta
      # Use Rouche-Capelli theorem to check condition
      X = as.matrix(X_mu[-j,])     # Convert to matrix to avoid rankMatrix throwing errors oin some cases
      rank = Matrix::rankMatrix( X, method = "qr", warn.t = FALSE)
      rank_c = Matrix::rankMatrix( cbind(X, y[-j]), method = "qr", warn.t = FALSE)
      if (rank == rank_c) {
        success = TRUE
      }
      else {
        i = i + 1
        if (i > ncol(Q)){
          i = i - 1
          at_end = TRUE
        }
      }
    }

    # Remove columns from X_sigma
    if (success) {
      # Calculate column in X_sigma which corresponds to column in Q
      if (inherits( qr_sigma, "sparseQR")){
        j = qr_sigma@q[zeros$ix[i]] + 1
      }
      else {
        j = zeros$ix[i]
      }
      X_sigma = X_sigma[, -j]
    }
    else {
      continue = FALSE
    }
  }

  # Scrutinize X_sigma futher by looking at its hat matrix
  continue = TRUE
  while(continue){

    hatvalues = stats::hat( X_sigma, intercept = FALSE)
    rows = which(hatvalues > 1 - tol)

    success = FALSE
    i = 1
    while(!success & i <= length(rows)){

      # Check that X_mu_1 * beta = y_1 has a solution for beta
      X = X_mu[rows[i],]
      if (!(all(X == 0) & y[rows[i]] != 0)){
        success = TRUE
      }
      else {
        i = i + 1
      }
    }

    if (success){
      X_sigma = make_matrix_full_rank.default(X_sigma[-rows[i],])
    }
    else {
      continue = FALSE
    }
  }

  # Restore column names
  column_numbers = as.integer(colnames(X_sigma))
  column_names = column_names[column_numbers]
  colnames(X_sigma) = column_names

  outlist = list( column_numbers = column_numbers,
                  column_names = column_names)

  return(outlist)
}
