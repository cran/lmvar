
<!-- README.md is generated from README.Rmd. Please edit that file -->
When to use the 'lmvar' package?
--------------------------------

The 'lmvar' package fits a linear model in which the variance can be different for different observations. The assumption of a constant variance inherent in a classical linear model, is dropped in 'lmvar'. I.e., 'lmvar' fits a heteroscedastic model whereas the classical linear model is homoscedastic. In those cases where the restriction of a constant variance is too restrictive in the linear model, 'lmvar' provides an alternative.

Hence 'lmvar' is an alternative for the function `lm` in the package 'stats'.

How to use the 'lmvar' package?
-------------------------------

The main function in the package is `lmvar`. Its basic usage requires three arguments:

-   a vector of response values
-   a model matrix for the expected values of the responses
-   a model matrix for the logarithms of the standard deviations of the responses.

As an example we use the data frame 'cats' in the package 'MASS'. We regress the cats heart weight 'Hwt' onto their body weight 'Bwt'.

``` r
require(MASS)
#> Loading required package: MASS
require(lmvar)
#> Loading required package: lmvar

# Create the model matrix. Do not include an intercept term: it will be added by 'lmvar'
X = model.matrix( ~ Bwt - 1, cats)

# Perform the fit. Use the same matrix for the expected values and the standard deviations
fit = lmvar( cats$Hwt, X, X)

# Print a summary of the fit
summary(fit)
#> Call: 
#>  lmvar(y = cats$Hwt, X_mu = X, X_sigma = X)
#> 
#> Standardized residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -2.4159 -0.7261 -0.0655  0.6703  2.5998 
#> 
#> Coefficients:
#>                Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)   -0.012179   0.679820 -0.0179  0.985707    
#> Bwt            3.904310   0.258964 15.0766 < 2.2e-16 ***
#> (Intercept_s) -0.522274   0.337044 -1.5496  0.121244    
#> Bwt_s          0.315837   0.121843  2.5922  0.009537 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Standard deviations: 
#>    Min     1Q Median     3Q    Max 
#> 1.1156 1.2265 1.3916 1.5422 2.0330 
#> 
#> Comparison to model with constant variance (i.e. classical linear model)
#> Log likelihood-ratio: 4.069774 
#> Additional degrees of freedom: 1 
#> p-value for difference in deviance: 0.00433
```

The low p-value in the summary is a strong indication that 'lmvar' does a better fit than 'lm'. The spread of the heart weight values grows with increasing body weight, resulting in a larger standard deviation of the heart weight with increasing body weight.

What is in the 'lmvar' package?
-------------------------------

The function `lmvar` produces an object of class 'lmvar'. In the example above, this is the object `fit`. The package provides a number of utility functions to extract information from an 'lmvar' object. The example above demonstrates the utility function `summary`. There are also utility functions to obtain the fitted coefficients, the expected values, the standard deviations, the log-likelihood, the degrees of freedom, the covariance matrix, the AIC value, confidence intervals, and more. To view the whole list of functions, use `help(package = "lmvar")`.

Further reading
---------------

The package comes with two vignettes: 'Intro' and 'Math'. 'Intro' gives an introduction to the package and is helpful reading for new users. 'Math' gives a detailed mathematical background to the package. The vignette 'Intro' can be viewed with `vignette( "Intro", package = "lmvar")`, and likewise for 'Math'.
