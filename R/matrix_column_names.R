matrix_column_names <- function( X, X_sigma){

  # Helper function to set column names of matrices

  colnames(X)[1] = "(Intercept)"
  l = ncol(X)
  if (l > 1){
    if (all(colnames(X)[2:l] == "")){
      colnames(X)[2:l] = sapply(2:l, function(x){paste('v', as.character(x-1), sep="")})
    }
  }

  colnames(X_sigma)[1] = "(Intercept_s)"
  l = ncol(X_sigma)
  if (l > 1){
    if (all(colnames(X_sigma)[2:l] == "")){
      colnames(X_sigma)[2:l] = sapply(2:l, function(x){paste('v', as.character(x-1), "_s", sep="")})
    }
  }

  return( list( colnames_X = colnames(X),
                colnames_X_sigma = colnames(X_sigma)))
}
