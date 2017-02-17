gaussian_var <- function(){

  # Function returns (part of) a family object for an 'lmvar' type of fit

  rlist = list( family = "gaussian with non-constant variances",
                link = "identity")
  return (structure( rlist, class = "family"))
}
