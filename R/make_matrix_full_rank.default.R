
make_matrix_full_rank.default <- function(X){

  # Function removes as many columns from the matrix as needed to make the matrix X full-rank

  # Input: X: matrix of class 'matrix' or any of the classes of the package 'Matrix'

  # Output: the full rank-matrix

  qX = Matrix::qr(X)

  make_matrix_full_rank( qX, X)

}
