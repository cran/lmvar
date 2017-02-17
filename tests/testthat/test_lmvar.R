context("lmvar constructor")

test_that("missing arguments are handled correctly", {

  y = fit$y
  fittest = suppressWarnings(lmvar(y))
  fitlm = lm( y~1, as.data.frame(fit$X_mu))
  expect_equal( coef(fittest)[1], coef(fitlm))
  expect_lt( abs(as.numeric(exp(coef(fittest)[2])) - summary(fitlm)$sigma), 0.002)

})

test_that("no errors occur when working with class Matrix",{

  M_mu = Matrix::Matrix(fit$X_mu)
  M_mu = M_mu[, -1]
  M_sigma = Matrix::Matrix(fit$X_sigma)
  M_sigma = M_sigma[, -1]
  fittest = suppressWarnings(lmvar( fit$y, M_mu, M_sigma))
  expect_equal( coef(fittest), coef(fit))

  expect_error( summary(fittest), NA)

})

test_that("results compare with standard linear regression",{

  y = fit$y
  M_mu = fit$X_mu
  M_mu = M_mu[, 2:ncol(M_mu)]

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
