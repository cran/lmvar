
make_matrix_full_rank.qr <- function( qX, X){

  # Function removes as many columns from the matrix as needed to make the matrix X full-rank

  # Input: qX: object of class 'qr'. It must contain the qr-decomposition of X.
  #        X: matrix of class 'matrix' or 'Matrix'

  # Output: the full rank-matrix

  if (qX$rank != 1){
    X = X[, qX$pivot[1:qX$rank]]
  }
  else {
    name = colnames(X)[qX$pivot[1]]
    X = as.matrix(X[, qX$pivot[1]])
    colnames(X) = name
  }
  return(X)
}
