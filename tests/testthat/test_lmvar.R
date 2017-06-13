context("lmvar constructor")

test_that("missing arguments are handled correctly", {

  y = fit$y
  fittest = lmvar(y)
  fitlm = lm( y~1, as.data.frame(fit$X_mu))
  expect_equal( coef(fittest)[1], coef(fitlm))
  expect_lt( abs(as.numeric(exp(coef(fittest)[2])) - summary(fitlm)$sigma), 0.002)

})

test_that("no errors occur when working with class Matrix",{

  M_mu = Matrix::Matrix(fit$X_mu)
  M_mu = M_mu[, -1]
  M_sigma = Matrix::Matrix(fit$X_sigma)
  M_sigma = M_sigma[, -1]
  fittest = lmvar( fit$y, M_mu, M_sigma)
  expect_equal( coef(fittest), coef(fit))

  expect_error( summary(fittest), NA)

})

test_that("results compare with standard linear regression",{

  y = fit$y
  M_mu = fit$X_mu
  M_mu = M_mu[,-1]

  fitlm = lm( y~., as.data.frame(as.matrix(M_mu)))
  fitlmvar = lmvar( y, M_mu)

  # compare betas and standard errors
  coeflm = coef(summary(fitlm))
  coeflmvar = as.matrix(coef(summary(fitlmvar)))
  coeflmvar = coeflmvar[ rownames(coeflm),]

  expect_equal( coeflm[,1], coeflmvar[,1])
  expect_equal( signif( coeflm[,2], 2), signif( coeflmvar[,2], 2))

  # Compare sigmas
  sigmalmvar = predict( fitlmvar, mu = FALSE)[1]
  variance = vcov( fitlmvar, mu = FALSE)[1,1]
  sdev = sigmalmvar * sqrt(exp(variance) - 1)
  sigmalm = summary(fitlm)$sigma
  expect_lt( sigmalmvar - sdev, sigmalm)
  expect_gt( sigmalmvar + sdev, sigmalm)
})

test_that("betas are as expected", {

  coeff = coef(summary(fit))

  beta = c( -2, 1, 3, -1.5, -1.1, 0.8, -0.5)

  v = sapply( 1:nrow(coeff), function(i){

    beta_f = coeff[i,1]
    sterr_f = coeff[i,2]
    expect_lt( beta_f - 1.5 * sterr_f, beta[i])
    expect_gt( beta_f + 1.5 * sterr_f, beta[i])
  })
})

test_that("no errors occur when matrix becomes vector", {

  M_mu = fit$X_mu[,2]
  M_sigma = fit$X_sigma[,2]

  expect_error( suppressWarnings(lmvar( fit$y, M_mu, M_sigma)), NA)

})

test_that("coefficient names are as expected", {

  X_mu = fit$X_mu[,-1]
  X_sigma = fit$X_sigma[,-1]

  names_mu = c("A", "B", "C")
  names_sigma = c("D", "E")
  colnames(X_mu) = names_mu
  colnames(X_sigma) = names_sigma

  fittest = lmvar( fit$y, X_mu, X_sigma)
  expect_equal(names(coef( fittest, sigma = FALSE)), c( "(Intercept)", names_mu))
  expect_equal(names(coef( fittest, mu = FALSE)), c( "(Intercept_s)", names_sigma))

  fittest = lmvar( fit$y, X_mu, X_sigma, intercept_mu = FALSE, intercept_sigma = FALSE)
  expect_equal(names(coef( fittest, sigma = FALSE)), names_mu)
  expect_equal(names(coef( fittest, mu = FALSE)), names_sigma)

  fittest = lmvar(fit$y)
  expect_equal(names(coef( fittest, sigma = FALSE)), "(Intercept)")
  expect_equal(names(coef( fittest, mu = FALSE)), "(Intercept_s)")
})

test_that("control options have effect", {

  # Request solver log
  y = fit$y
  fittest = lmvar(y, control = list(slvr_log = TRUE))
  expect_true("slvr_log" %in% names(fittest))

  # Fit model without solution
  c1 = rep(1, 7)
  c2 = c( 0, 1, 1, 0, 0, 0, 0)
  c3 = c( 0, 0, 0, 1, 1, 0, 0)
  c4 = c( 0, 0, 0, 0, 0, 1, 1)
  X = matrix( c( c1, c2, c3, c4), nrow = 7)

  set.seed(5678)
  y = rnorm( 7, mean = 1)
  beta_start = c(1, -1, -1, -1)

  # Test that warning about Hessian appears
  fittest = capture_warnings(lmvar( y, X_sigma = X, intercept_sigma = FALSE, slvr_options = list(start = -10 * beta_start),
                                    control = list(running_diagnostics = FALSE)))
  expect_true(any(grepl( "Log-likelihood appears not to be at a maximum!", fittest)))

  # Test that warning about Hessian does not appear
  fittest = capture_warnings(lmvar( y, X_sigma = X, intercept_sigma = FALSE, slvr_options = list(start = -10 * beta_start),
                            control = list( check_hessian = FALSE, running_diagnostics = FALSE)))
  expect_false(any(grepl( "Log-likelihood appears not to be at a maximum!", fittest)))
})
