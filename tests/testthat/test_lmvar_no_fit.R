context("lmvar_no_fit contructor")

test_that("lmvar_no_fit is identical to the subclass of lmvar", {

  no_fit = lmvar_no_fit( fit$y, fit$X_mu[,-1], fit$X_sigma[,-1], intercept_mu = fit$intercept_mu,
                         intercept_sigma = fit$intercept_sigma, slvr_options = fit$slvr_options, control = fit$control)

  check_names = names(no_fit)
  i = which(check_names == "call")
  check_names = check_names[-i]

  for (el in check_names){
    expect_identical( no_fit[[el]], fit[[el]])
  }
})
