
make_matrix_full_rank <- function( X){

  # Function removes as many columns from the matrix as needed to make the matrix X full-rank

  # Input: X: matrix

  # Output: the full rank-matrix

  qX = qr(as.matrix(X))

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
