Version 1.2.1
-------------

* Fix show-stopper that prevents 'cv.lmvar' from running.

Version 1.2.0
-------------

* The function 'cv.lmvar' is added to carry out cross-validations. For comparison with a classical linear model, 
the analogous function 'cv.lm' is provided. 

* Arguments of the function 'lmvar' not related to the definition of the model, are bundled in a single argument 'control' of class 'list'.

* The function 'maxLik' from the package 'maxLik' is used to find the maximum likelihood. This allows one to use all methods (exept for 'BHHH') and options supported by 'maxLik'. 

* Occurrence of errors in the calculation of the log-likelihood during the iterative search for the maximum log-liklihood, no longer stop the execution of the function 'lmvar'. Corrective measures and warnings are left to 'maxLik'.

* All arguments of the function 'lmvar' are included in the output object.

* The summary-overview prints the number of observations and the degrees of freedom of the fit.


Version 1.1.0
-------------

* To find the maximum likelihood, a set of non-linear equations must be solved. This is done by the function 'maxNR' from the package 'maxLik'. It results in faster and more robust solves, compared to the previous version of 'lmvar'. 

* It is possible to pass on options to 'maxNR' from the call to 'lmvar'. For this, 'lmvar' has a new argument 'slvr_options'.

* It is possible to request the result log of 'maxNR'. For this, 'lmvar' has a new argument 'slvr_log'.

* It is possible to suppress the intercept terms in the model for mu and the model for log sigma. For this, 'lmvar' has the new arguments 'intercept_mu' and 'intercept_sigma'.

* The functions 'fitted.lmvar' and 'predict.lmvar' support the interval type 'prediction'.

* Input matrices X_mu and X_sigma to 'lmvar' that have only one column, can be of type 'numeric'.

* The matrix of coefficients shown by 'summary', can be restricted to only the coefficients 'beta_mu' or only the coefficients 'beta_sigma'. For this, 'summary.lmvar' has the new arguments 'mu' and 'sigma'. This option is useful when the length of the vector 'beta_mu' and/or 'beta_sigma' is large.

* The function 'df.residuals.lmvar' has been removed. The 'residual degrees of freedom' are a useful concept in a classical linear model where it specifies the degrees of freedom of the Student t-distribution for e.g. confidence intervals of beta. This distribution plays no role in the model of 'lmvar' which uses aymptotically normal distributions.

* Improved documentation, in particular the README and the vignettes 'Intro' and 'Math'

* Bug fixes

Version 1.0.0
-------------

First release-version.
