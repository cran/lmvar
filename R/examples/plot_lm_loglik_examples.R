\dontrun{

# Carry out a linear regression with the 'iris' data set
fit = lm( Petal.Length ~ Species, data = iris, x = TRUE, y = TRUE)
X = fit$x
y = fit$y

# We center the plot at the maximum-likelihood
beta_or = coef(fit)

# Plot the maximum log-likelihood
lmvar:::plot_lm_loglik( y, X, beta_or = beta_or, beta_x = "(Intercept)",
                        beta_y = "Speciesversicolor")

# Plot against the two species
lmvar:::plot_lm_loglik( y, X, beta_or = beta_or, beta_x = "Speciesversicolor",
                        beta_y = "Speciesvirginica")

# Increase the resolution
lmvar:::plot_lm_loglik( y, X, beta_or = beta_or, beta_x = "Speciesversicolor",
                        beta_y = "Speciesvirginica", plot_points = 40)

# Remove the intercept term from the model matrix and fit again
XX = X[,-1]
fit = lm( y ~ . - 1, data = as.data.frame(XX))

# Estimate the effect of adding an intercept term in a quadratic approximation and compare
# with exact result
beta_or = c( 0, coef(fit))
lmvar:::plot_lm_loglik( y, X, beta_or = beta_or, beta_x = 1, beta_y = "Speciesversicolor",
                        add_qa = TRUE, plot_points = 40, plot_width = 5)

# Note that in the last case the quadratic approximation has no maximum. Hence the beta for
# "Speciesvirginica" is kept at beta_or[3] in the calculation of the surface of the
# quadratic approximation.
}
