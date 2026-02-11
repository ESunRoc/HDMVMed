`HDMVMed` is an R package for fitting and performing inference on the high-dimensional multivariate 
mediation model of Sun et al. (2026), developed by the authors.
`HDMVMed` provides estimates for the causal direct effects (DEs), total indirect effects (TIDEs), and
mediator-wise partial indirect effects (PIDEs). Bootstrap p-values as well as BCa and 
percentile-based bootstrap confidence intervals are additionally provided following the framework of
Sun et al. (2026). For the simulation code and OMEI application results discussed in the paper,
please refer to the repo [https://github.com/ESunRoc/HDMVMediation](https://github.com/ESunRoc/HDMVMediation) 
and accompanying GitHub website, [esunroc.github.io/HDMVMediation/](esunroc.github.io/HDMVMediation/).


# Install
`HDMVMed` is installable via GitHub:
```
devtools::install_github("ESunRoc/HDMVMed")
```

Note that `HDMVMed` makes heavy usage of the `MSGLasso` package of Li et al. (2016). This package
was removed from CRAN and therefore needs to be installed from the source files, which are available
via the CRAN archive. This is true for `coxed` package, as well, from which we use the `bca` function. 
Specifically, download the tarballs [`MSGLasso_2.1.tar.gz`](https://cran.r-project.org/src/contrib/Archive/MSGLasso/) and [`coxed_0.3.3.tar.gz`](https://cran.r-project.org/src/contrib/Archive/coxed/)
and then run
```
# install.packages("rtools")
install.packages(path_to_MSGLasso_2.1.tar.gz, repos = NULL, type = "source")
install.packages(path_to_coxed_0.3.3.tar.gz, repos = NULL, type = "source")
```
to install the necessary packages.



# References

Li, Y., Nan, B., and Zhu, J. (2016). MSGLasso: Multivariate Sparse Group Lasso for the Multivariate Multiple Linear Regression with an Arbitrary
Group Structure. R package version 2.1, <https://CRAN.R-project.org/package=MSGLasso>.

Sun, E., Xiao, J., and Wu, T.T. (2026). Causal Mediation Analysis for Multiple Outcomes and High-dimensional Mediators: 
Identification, Inference, and Application. Biometrics. Submitted.
