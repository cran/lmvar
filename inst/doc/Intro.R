## ------------------------------------------------------------------------
library(lmvar)

# As example we use the dataset 'attenu' from the library 'datasets'. The dataset contains
# the response variable 'accel' and two explanatory variables 'mag'  and 'dist'.
library(datasets)

# Create the model matrix for the expected values
X = cbind(attenu$mag, attenu$dist)
colnames(X) = c("mag", "dist")

# Create the model matrix for the standard deviations.
X_s = cbind(attenu$mag, 1 / attenu$dist)
colnames(X_s) = c("mag", "dist_inv")

# Carry out the fit
fit = lmvar(attenu$accel, X, X_s)

## ------------------------------------------------------------------------
summary(fit)

## ------------------------------------------------------------------------
sigma = fitted(fit, mu = FALSE)
hist(sigma)

## ------------------------------------------------------------------------
hist(residuals(fit))

## ------------------------------------------------------------------------
dfree(fit, sigma = FALSE)

## ------------------------------------------------------------------------
nobs(fit)
logLik(fit)
AIC(fit)

## ------------------------------------------------------------------------
coef(fit)

## ------------------------------------------------------------------------
coef(fit, mu = FALSE)

## ------------------------------------------------------------------------
vcov(fit, sigma = FALSE)

