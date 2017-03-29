matrix_column_names <- function( X, X_sigma, intercept_mu = TRUE, intercept_sigma = TRUE){

  # Helper function to set column names of matrices

  # Input: X: model matrix for expected values
  #        X_sigma: model matrix for log of standard deviations
  #        X_intc_only: Boolean, if TRUE the matix X consists of an intercept term only
  #        X_sigma_intc_only: Boolean, if TRUE the matix X_sigma consists of an intercept term only
  #        intercept_mu: Boolean, if TRUE the matrix X contains an intercept term. This must be the first column
  #        intercept_sigma: Boolean, if TRUE the matrix X_sigma contains an intercept term. This must be the first column

  # set the names of the matrix X
  i_first = 1
  if (intercept_mu){
    colnames(X)[1] = "(Intercept)"
    i_first = 2
  }
  l = ncol(X)
  if (l >= i_first){
    if (all(colnames(X)[i_first:l] == "")){
      colnames(X)[i_first:l] = sapply(1:(l-i_first+1), function(x){paste('v', as.character(x), sep="")})
    }
  }

  # Set the names of the matrix X_sigma
  i_first = 1
  if (intercept_sigma){
    colnames(X_sigma)[1] = "(Intercept_s)"
    i_first = 2
  }
  l = ncol(X_sigma)
  if (l >= i_first){
    if (all(colnames(X_sigma)[i_first:l] == "")){
      colnames(X_sigma)[i_first:l] = sapply(1:(l-i_first+1), function(x){paste('v', as.character(x), "_s", sep="")})
    }
  }

  return( list( colnames_X = colnames(X),
                colnames_X_sigma = colnames(X_sigma)))
}
