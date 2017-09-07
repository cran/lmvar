make_matrix_full_rank.sparseQR <- function( qX, X){

  # Function removes as many columns from the matrix as needed to make the matrix X full-rank

  # Input: qX: object of class 'sparseQR'. It must contain the qr-decomposition of X.
  #        X: matrix of class 'matrix' or 'Matrix'

  # Output: the full rank-matrix


  tol = max(dim(X)) * .Machine$double.eps

  qR = suppressWarnings(Matrix::qr.R(qX))
  diagR = Matrix::diag(qR)

  tol = tol * max(diagR)
  rank = sum(abs(diagR) >= tol)

  columns = sort(diagR, decreasing = TRUE, index.return = TRUE)
  columns = columns$ix[1:rank]
  columns = sort.int(qX@q[columns]) + 1
  return(X[, columns])
}
