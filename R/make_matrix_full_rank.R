
make_matrix_full_rank <- function( X){

  # Function removes as many columns from the matrix as needed to make the matrix X full-rank

  # Input: X: matrix

  # Output: the full rank-matrix

  qX = qr(as.matrix(X))

  X = X[, qX$pivot[1:qX$rank]]

  return(X)
}
