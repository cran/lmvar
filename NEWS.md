Version 1.4.0
-------------

* Introduce the class 'lmvar_no_fit'. This class is like the class 'lmvar', but without members that
are the result of a model fit. The constructor of an object of class 'lmvar_no_fit' is the function 'lmvar_no_fit()'.  
  The class 'lmvar' is an extension of the class 'lmvar_no_fit'. This means that wherever an object of class 'lmvar_no_fit' is required, an object of class 'lmvar' can be used as well.   
  The class 'lmvar_no_fit' was motivated by situations in which 'lmvar()' does not converge for a model with many degrees of freedom
and one wants to resort to 'fwbw()' to obtain a subset of degrees of freedom which does converge.  

* Modify appropriate functions (such as 'dfree()' and 'fwbw.lmvar()') such that they take a 'lmvar_no_fit' object as input.

* Add control option 'remove_df_sigma' to function 'lmvar()'. In cases where the 'lmvar' fit does not converge, switching on this option
may restore convergence. Update both vingettes to give information about this option. 

* Change the algorithm of the functions 'fwbw.lm()' and 'fwbw.lmvar_no_fit()'. The change makes the insert and remove step of the algorithm symmetric. Rather than inserting degrees of freedom one-by-one, it will attempt to insert a percentage of the left-out degrees of freedom, just like the algorithm attempts to remove a percentage of the degrees of freedom that are included.  
  This makes it feasible to run the functions with the option `fw = TRUE`, even for models with potentially many degrees of freedom. However, the outcome of the functions can be different compared to previous versions of the 'lmvar' package.

* Rename control option 'running_diagnostics' to 'monitor' for function 'lmvar()'. This is consistent with the function 'fwbw()'. 

* Make calculation of $\beta_\mu$ in 'lmvar()' more robust.
 
* Minor change in print format of object of class 'cvlmvar' to be more consistent with accepted terminology.

* Fix bug which can cause intercept term to be removed when matrix of class 'Matrix' is made full-rank. 

* Fix bug in 'lmvar()' which can cause wrong beta to be reported during run when control option 'monitor' is set to TRUE.

* Fix bug in 'fwbw.lmvar_no_fit()' which causes an error when an initial estimate for $\beta_\sigma$ was specified in `object`.

Version 1.3.0
-------------

* New functions 'fwbw.lm' and 'fwbw.lmvar' for model-selection by means of a forward- / backward-step algorithm. 

* A user-specified function can be cross-validated in 'cv.lm' and 'cv.lmvar'.

* The function 'lmvar' allows one to solve the model under the contraint of minimum standard deviations for all
observations.

* The k fits in a k-fold cross-validation are executed in parallel to boost performance in case of large models (in 'cv.lm' and 'cv.lmvar').

* Performance improvements for large sparse model matrices of class 'dgCMatrix'

* Fix bug in 'cv.lm' and 'cv.lmvar' in case model matrices have one column only.

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
