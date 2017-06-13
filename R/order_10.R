
order_10 <- function(x){

  # Function calculates the integer n such that x / 10^n is in the interval (-10,10)
  #
  # Input: x: a single numeric value
  #
  # Output: n

  n = 0
  x = abs(x)
  if (x <= 1){
    while (x < 1){
      n = n - 1
      x = x * 10
    }
  }
  else {
    while(x >= 10){
      n = n + 1
      x = x / 10
    }
  }
  return(n)
}
