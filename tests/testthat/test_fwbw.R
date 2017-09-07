context("Forward-backward model selection")

test_that("Selection retrieves original model (lmvar)", {

  set.seed(8901)

  # create 'lmvar' object
  X = fit$X_mu
  X = cbind( X, matrix( sample( -9:9, 6 * nrow(X), replace = TRUE), ncol = 6))
  colnames(X) = c( colnames(X)[1:4], paste( "v", 4:9, sep = ""))

  X_sigma = fit$X_sigma
  X_sigma = cbind( X_sigma, matrix( sample( -1:1, 3 * nrow(X), replace = TRUE), ncol = 3))
  colnames(X_sigma) = c( colnames(X_sigma)[1:3], paste( "v", 3:5, sep = ""))

  beta = c( coef( fit, sigma = FALSE), c( 2, -3, 0.5, 1.5, 2.5, -1))
  mu = X %*% beta

  beta_sigma = 0.5 * c( coef( fit, mu = FALSE), c( 2, -3, 0.5))
  sigma = exp( X_sigma %*% beta_sigma)

  # add spurious terms to model matrix
  X_spur = matrix( sample( -9:9, 5 * nrow(X), replace = TRUE), ncol = 5)
  colnames(X_spur) = paste( "vspur", 1:5, sep = "")
  X_spur = cbind( X, X_spur)

  X_sigma_spur = matrix( sample( -1:1, 3 * nrow(X), replace = TRUE), ncol = 3)
  colnames(X_sigma_spur) = paste( "vspur", 1:3, sep = "")
  X_sigma_spur = cbind( X_sigma, X_sigma_spur)

  y = rnorm( nrow(X), mean = mu, sd = sigma)

  fit_lmvar = lmvar( y, X[,-1], X_sigma[,-1])
  fit_test = lmvar(y, X_spur[,-1], X_sigma_spur[,-1])

  # Forward-stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$X_mu)
  expect_equal( cols, colnames(X))

  cols = colnames(fwbw$object$X_sigma)
  expect_equal( cols, colnames(X_sigma))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lmvar))

  # Backward-stepping
  fwbw = fwbw( fit_test, BIC)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$X_mu)
  expect_equal( cols, colnames(X))

  cols = colnames(fwbw$object$X_sigma)
  expect_equal( cols, colnames(X_sigma))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lmvar))
})

test_that("Selection retrieves original model (lm and lmvas as lm)", {

  set.seed(1234)

  # create 'lm' object
  X = fit$X_mu
  X = cbind( X, matrix( sample( -9:9, 6 * nrow(X), replace = TRUE), ncol = 6))
  colnames(X) = c( colnames(X)[1:4], paste( "v", 4:9, sep = ""))

  beta = c( coef( fit, sigma = FALSE), c( 2, -3, 0.5, 1.5, 2.5, -1))
  mu = X %*% beta

  y = rnorm( nrow(X), mean = mu, sd = 0.5)

  fit_lm = lm( y ~., as.data.frame(X[,-1]))

  # add spurious terms to model matrix
  X_spur = matrix( sample( -9:9, 5 * nrow(X), replace = TRUE), ncol = 5)
  colnames(X_spur) = paste( "vspur", 1:5, sep = "")
  X_spur = cbind( X, X_spur)
  fit_test = lm(y ~ ., as.data.frame(X_spur[,-1]), x = TRUE, y = TRUE)

  # Backward-stepping
  fwbw = fwbw( fit_test, BIC)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$x)
  expect_equal( cols, colnames(X))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lm))

  # Forward-stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$x)
  expect_equal( cols, colnames(X))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lm))

  # check lmvar as well
  fit_test = lmvar( y, X_spur[,-1])

  # Backward-stepping
  fwbw = fwbw( fit_test, BIC)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$X_mu)
  expect_equal( cols, colnames(X))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lm))

  # Forward-stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$X_mu)
  expect_equal( cols, colnames(X))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lm))
})

test_that("Case in which all degrees of freedom are removed (lmvar)", {

  v1_unif = runif( nobs(fit), -10, 10)
  v2_unif = runif( nobs(fit), -10, 10)
  v3_unif = runif( nobs(fit), -10, 10)
  v4_unif = runif( nobs(fit), -10, 10)
  X_mu = cbind( v1_unif, v2_unif, v3_unif, v4_unif)

  X_sigma = model.matrix( ~ . - 1 + .*., data = as.data.frame(fit$X_sigma[,-1]))

  fit_test = lmvar( fit$y, X_mu, X_sigma)
  fit_expected = lmvar(fit$y)

  fwbw = fwbw( fit_test, BIC)

  expect_equal( fwbw$object$coefficients_mu, fit_expected$coefficients_mu)
  expect_equal( fwbw$object$coefficients_sigma, fit_expected$coefficients_sigma)
  expect_equal( fwbw$object$logLik, fit_expected$logLik)
  expect_equivalent( fwbw$object$X_mu, fit_expected$X_mu)
  expect_equivalent( fwbw$object$X_sigma, fit_expected$X_sigma)
  expect_equal( fwbw$fun, BIC(fit_expected))
})

test_that("Case in which all degrees of freedom are removed (lm)", {

  v1_unif = runif( nobs(fit), -10, 10)
  v2_unif = runif( nobs(fit), -10, 10)
  v3_unif = runif( nobs(fit), -10, 10)
  v4_unif = runif( nobs(fit), -10, 10)
  X_mu = cbind( v1_unif, v2_unif, v3_unif, v4_unif)

  fit_test = lm( fit$y ~ ., as.data.frame(X_mu), x = TRUE, y = TRUE)
  fit_expected = lm(fit$y ~ 1, x = TRUE)

  fwbw = fwbw( fit_test, BIC)

  expect_equal( fwbw$object$coefficients, fit_expected$coefficients)
  expect_equal( logLik(fwbw$object), logLik(fit_expected))
  expect_equivalent( fwbw$object$x, fit_expected$x)
  expect_equal( fwbw$fun, BIC(fit_expected))
})

test_that("No intercept term (lmvar)", {

  # Generate response vector for model without intercept terms
  X_mu = fit$X_mu[,-1]
  X_sigma = fit$X_sigma[,-1]
  beta_mu = coef( fit, sigma = FALSE)[-1]
  beta_sigma = coef( fit, mu = FALSE)[-1]

  mu = as.numeric(X_mu %*% beta_mu)
  sigma = exp(as.numeric(X_sigma %*% beta_sigma))
  y = rnorm( nobs(fit), mean = mu, sd = sigma)

  fit_lmvar = lmvar( y, X_mu, X_sigma, intercept_mu = FALSE, intercept_sigma = FALSE)

  # add spurious interaction terms to model matrices
  X_mu_spur = model.matrix( ~ . - 1 + .*., data = as.data.frame(X_mu))
  X_sigma_spur = model.matrix( ~ . - 1 + .*., data = as.data.frame(X_sigma))
  fit_test = lmvar( y, X_mu_spur, X_sigma_spur, intercept_mu = FALSE, intercept_sigma = FALSE)

  # backward stepping
  fwbw = fwbw( fit_test, BIC)

  # check whether original model has been retrieved
  cols_mu = colnames(fwbw$object$X_mu)
  cols_sigma = colnames(fwbw$object$X_sigma)
  expect_equal( cols_mu, colnames(X_mu))
  expect_equal( cols_sigma, colnames(X_sigma))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lmvar))

  # forward stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)

  # check whether original model has been retrieved
  cols_mu = colnames(fwbw$object$X_mu)
  cols_sigma = colnames(fwbw$object$X_sigma)
  expect_equal( cols_mu, colnames(X_mu))
  expect_equal( cols_sigma, colnames(X_sigma))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lmvar))
})

test_that("No intercept term (lm)", {

  # Generate response vector for model without intercept terms
  X_mu = fit$X_mu[,-1]
  beta_mu = coef(fit, sigma = FALSE)[-1]

  mu = as.numeric(X_mu %*% beta_mu)
  y = rnorm( nobs(fit), mean = mu, sd = 0.5)

  fit_lm = lm( y ~ . - 1, as.data.frame(X_mu))

  # add spurious interaction terms to model matrices
  X_mu_spur = model.matrix( ~ . - 1 + .*., data = as.data.frame(X_mu))
  fit_test = lm( y ~ . - 1, as.data.frame(X_mu_spur), x = TRUE, y = TRUE)

  # backward stepping
  fwbw = fwbw( fit_test, BIC)

  # check whether original model has been retrieved
  cols_mu = colnames(fwbw$object$x)
  expect_equal( cols_mu, colnames(X_mu))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lm))

  # forward stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)

  # check whether original model has been retrieved (which is not the case)
  cols_mu = colnames(fwbw$object$x)
  expect_equal( cols_mu, c( "v1", "v3"))
})

test_that("Unusual variable names (lm)", {

  set.seed(1234)

  # create 'lm' object
  X = fit$X_mu
  X = cbind( X, matrix( sample( -9:9, 6 * nrow(X), replace = TRUE), ncol = 6))
  colnames(X) = c( colnames(X)[1:4], "(v4)", "_v5_", "=v6=", " v7 ", "^v8^", "0v90")

  beta = c( coef( fit, sigma = FALSE), c( 2, -3, 0.5, 1.5, 2.5, -1))
  mu = X %*% beta

  y = rnorm( nrow(X), mean = mu, sd = 2.5)
  fit_lm = lm (y ~ ., as.data.frame(as.matrix(X[,-1])))

  # add spurious terms to model matrix
  X_spur = model.matrix( ~ . + 0 + v1 * ., as.data.frame(as.matrix(X[,-1])))
  fit_test = lm(y ~ ., as.data.frame(X_spur), x = TRUE, y = TRUE)

  # forward stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)
  expect_equal( BIC(fwbw$object), BIC(fit_lm))

  # backward stepping
  fwbw = fwbw( fit_test, BIC)
  expect_equal( BIC(fwbw$object), BIC(fit_lm))

  # remove intercept
  X = X[,-1]
  beta = beta[-1]
  mu = X %*% beta

  y = rnorm( nrow(X), mean = mu, sd = 2.5)
  fit_lm = lm (y ~ . + 0, as.data.frame(as.matrix(X)))

  # add spurious terms to model matrix
  X_spur = model.matrix( ~ . + 0 + v1 * ., as.data.frame(as.matrix(X)))
  fit_test = lm(y ~ . + 0, as.data.frame(X_spur), x = TRUE, y = TRUE)

  # forward stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)
  expect_equal( BIC(fwbw$object), BIC(fit_lm))

  # backward stepping
  fwbw = fwbw( fit_test, BIC)
  expect_equal( BIC(fwbw$object), BIC(fit_lm))

})

test_that("Unusual variable names (lmvar)", {

  set.seed(8901)

  # create 'lmvar' object
  X = fit$X_mu
  X = cbind( X, matrix( sample( -9:9, 6 * nrow(X), replace = TRUE), ncol = 6))
  colnames(X) = c( colnames(X)[1:4], "(v4)", "_v5_", "=v6=", " v7 ", "^v8^", "0v90")

  X_sigma = fit$X_sigma
  X_sigma = cbind( X_sigma, matrix( sample( -1:1, 3 * nrow(X), replace = TRUE), ncol = 3))
  colnames(X_sigma) = c( colnames(X_sigma)[1:3], "(v3)", "_v4_", "=v5=")

  beta = c( coef( fit, sigma = FALSE), c( 2, -3, 0.5, 1.5, 2.5, -1))
  mu = X %*% beta

  beta_sigma = 0.5 * c( coef( fit, mu = FALSE), c( 2, -3, 0.5))
  sigma = exp( X_sigma %*% beta_sigma)

  # add spurious terms to model matrix
  X_spur = matrix( sample( -9:9, 5 * nrow(X), replace = TRUE), ncol = 5)
  colnames(X_spur) = paste( "vspur", 1:5, sep = "")
  X_spur = cbind( X, X_spur)

  X_sigma_spur = matrix( sample( -1:1, 3 * nrow(X), replace = TRUE), ncol = 3)
  colnames(X_sigma_spur) = paste( "vspur", 1:3, sep = "")
  X_sigma_spur = cbind( X_sigma, X_sigma_spur)

  y = rnorm( nrow(X), mean = mu, sd = sigma)

  fit_lmvar = lmvar( y, X[,-1], X_sigma[,-1])
  fit_test = lmvar(y, X_spur[,-1], X_sigma_spur[,-1])

  # Forward-stepping
  fwbw = fwbw( fit_test, BIC, fw = TRUE)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$X_mu)
  expect_equal( cols, colnames(X))

  cols = colnames(fwbw$object$X_sigma)
  expect_equal( cols, colnames(X_sigma))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lmvar))

  # Backward-stepping
  fwbw = fwbw( fit_test, BIC)

  # check whether original model has been retrieved
  cols = colnames(fwbw$object$X_mu)
  expect_equal( cols, colnames(X))

  cols = colnames(fwbw$object$X_sigma)
  expect_equal( cols, colnames(X_sigma))

  # check minimum BIC value
  expect_equal( fwbw$fun, BIC(fit_lmvar))
})
