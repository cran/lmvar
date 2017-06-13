# Create an object of class 'lm'. We use an arbitrary model matrix, an arbitrary parameter
# vector beta and a generated response vector y for the purpose of the example.
library(MASS)

X = model.matrix(~ Sex + Bwt, cats)
beta_mu = c(-0.1, 0.3, 4)

mu = X %*% beta_mu

y = rnorm( nrow(X), mean = mu, sd = 0.5)

fit = lm(y ~ ., as.data.frame(X[,-1]), x = TRUE, y = TRUE)

# Carry out a cross-validation
cv = cv.lm(fit)

# Carry out a cross-validation including a Kolmogorov-Smirnov test
cv = cv.lm(fit, ks_test = TRUE)

# Carry out a cross-validation with 5 folds and control the random numbers used
cv = cv.lm(fit, k = 5, seed = 5483)

# Use option 'log = TRUE' if you fit the log of the response vector and require error estimates for
# the response vector itself

fit = lm(log(y) ~ ., as.data.frame(X[,-1]), x = TRUE, y = TRUE)
cv = cv.lm(fit, log = TRUE)

# Print 'cv' using the print-method print.cvlmvar
cv

# Print 'cv' with a specified number of digits
print(cv, digits = 2)
