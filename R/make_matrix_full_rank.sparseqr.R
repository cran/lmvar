make_matrix_full_rank.sparseQR <- function( qX, X){

  # Function removes as many columns from the matrix as needed to make the matrix X full-rank

  # Input: qX: object of class 'sparseQR'. It must contain the qr-decomposition of X.
  #        X: matrix of class 'matrix' or 'Matrix'

  # Output: the full rank-matrix

  #  Check whether intercept term is present
  intercept = FALSE
  if (colnames(X)[1] %in% c( "(Intercept)", "(Intercept_s)")){
    intercept = TRUE
    col_name = colnames(X)[1]
  }

  tol = max(dim(X)) * .Machine$double.eps

  qR = suppressWarnings(Matrix::qr.R(qX))
  diagR = Matrix::diag(qR)

  tol = tol * max(diagR)
  rank = sum(abs(diagR) >= tol)

  columns = sort(diagR, decreasing = TRUE, index.return = TRUE)
  columns = columns$ix[1:rank]
  columns = sort.int(qX@q[columns]) + 1

  X = X[,columns]

  # Check whether intercept term is amongst selected columns
  if (!intercept | columns[1] == 1){
    return(X)
  }
  else {
    beta = Matrix::solve( Matrix::t(X) %*% X, Matrix::colSums(X))
    beta = as.numeric(beta)
    i = which.max(abs(beta))
    X = cbind( 1, X[,-i])
    colnames(X)[1] = col_name
    return(X)
  }
}

