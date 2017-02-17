# As example we use the dataset 'attenu' from the library 'datasets'. The dataset contains
# the response variable 'accel' and two explanatory variables 'mag'  and 'dist'.
library(datasets)

# For more info on the data, study the dataset
help("attenu")

# Create the model matrix for the expected values
X = cbind(attenu$mag, attenu$dist)
colnames(X) = c("mag", "dist")

# Create the model matrix for the standard deviations. The standard deviation
# is large for small distances and small for large distances. The use of 'dist'
# as explanatory variable makes the beta for the intercept term blow up.
# Therefore we use '1 / dist' as explanatory variable
X_s = cbind(attenu$mag, 1 / attenu$dist)
colnames(X_s) = c("mag", "dist_inv")

# Carry out the fit
fit_lmvar  = lmvar(attenu$accel, X, X_s)

# Inspect the results. Note from the p-value for the difference in
# deviance that this fit appears to be significantly better than
# a classical linear fit
summary(fit_lmvar)

# Carry out a classical linear fit for comparison
fit_lm = lm(attenu$accel ~ mag + dist, attenu)

# A comparison of the AIC values also favours the fit with 'lmvar'
AIC(fit_lm)
AIC(fit_lmvar)
