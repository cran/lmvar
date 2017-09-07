
round_to_sd_accuracy <- function( values, sd, digits){

  # Function rounds 'values' to a number of digits given by:
  # 1) if sd is smaller than one: the first significant digit of 'sd' plus 'digits', or
  # 2) if sd is larger than one: the first power of 10 of 'sd' plus 'digits'.
  #
  # Input: values: Numeric, can be a single value or a vector of values
  #        sd:    Single numeric value, must be stricly positive
  #        digits: Integer, number of digits. Must be larger than zero
  #
  # Output: Numeric vector with the rounded values of 'values'.

  if (is.na(sd)){
    return(rep( NA, length(values)))
  }
  else {
    n = - order_10(sd)
    return(round( values, n + digits - 1))
  }
}
