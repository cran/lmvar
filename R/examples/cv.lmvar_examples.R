# Create an object of class 'lmvar'. We use a model matrix obtained from the 'cats' dataframe,
# arbitrary parameter vectors beta and a generated response vector y for the purpose of the
# example.
library(MASS)

X = model.matrix(~ Sex + Bwt, cats)
beta_mu = c(-0.1, 0.3, 4)
beta_sigma = c(-0.5, -0.1, 0.3)

mu = X %*% beta_mu
log_sigma = X %*% beta_sigma

y = rnorm( nrow(X), mean = mu, sd = exp(log_sigma))

fit = lmvar(y, X_mu = X[,-1], X_sigma = X[,-1])

# Carry out a cross-validation
cv = cv.lmvar(fit)

# Carry out a cross-validation including a Kolmogorov-Smirnov test
cv = cv.lmvar(fit, ks_test = TRUE)

# Carry out a cross-validation with 5 folds and control the random numbers used
cv = cv.lmvar(fit, k = 5, seed = 5483)

# Carry out a cross-validation and exclude observations 5, 11 and 20 from the calculation of
# the error statistics
cv = cv.lmvar(fit, exclude = c(5, 11, 20))

# Carry out a cross-validation and specify the maximization routine and maximum number of iterations
cv = cv.lmvar(fit, slvr_options = list( method = "NR", control = list(iterlim = 500)))

# Use option 'log = TRUE' if you fit the log of the response vector and require error estimates for
# the response vector itself

fit = lmvar(log(y), X_mu = X[,-1], X_sigma = X[,-1])
cv = cv.lmvar(fit, log = TRUE)

# Print 'cv' using the print-method print.cvlmvar
cv

# Print 'cv' with a specified number of digits
print(cv, digits = 2)

