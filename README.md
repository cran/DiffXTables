The 'DiffXTables' R package
===============================

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/DiffXTables)](https://cran.r-project.org/package=DiffXTables)
[![CRAN_latest_release_date](https://www.r-pkg.org/badges/last-release/DiffXTables)](https://cran.r-project.org/package=DiffXTables)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/DiffXTables)](https://cran.r-project.org/package=DiffXTables)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/DiffXTables)](https://cran.r-project.org/package=DiffXTables)



### Overview

Pattern heterogeneity between two variables across conditions is often fundamental to a scientific inquiry. For example, a biologist could ask whether co-expression between two genes in a cancer cell has been modified from a normal cell. The 'DiffXTables' R package answers such questions via evaluating statistical evidence for distributional changes in the involved variables based on observed data without using a parametric mathematical model.

The package provides statistical methods for hypothesis testing of differences in the underlying distributions across two or more contingency tables. They include five statistical tests: 

1. The comparative chi-squared test: Its statistical foundation is first established in [(Song et al. 2014) <doi:10.1093/nar/gku086>](https://doi.org/10.1093/nar/gku086). It is later extended to identify differential patterns in networks [(Zhang et al. 2015) <doi:10.1093/nar/gkv358>](https://doi.org/10.1093/nar/gkv358).
2. The Sharma-Song test: The test detects differential departure from independence via second-order difference in the joint distributions underlying two or more contingency tables. The test is fully described in [(Sharma et al. 2021) <doi:10.1093/bioinformatics/btab240>](https://doi.org/10.1093/bioinformatics/btab240).
3. The heterogeneity test: It is described in (Zar, 2010) and widely appears in textbooks. In contrast to the above comparative chi-squared test, it is not always powerful as demonstrated by examples in the package vignette.
4. The marginal-change test: It is determines the first-order (marginal) differences across conditions [(Sharma et al. 2020) <doi:10.1145/3388440.3412485>](https://doi.org/10.1145/3388440.3412485).
5. The strength test: It determines the strength of association of all the conditions [(Sharma et al. 2020) <doi:10.1145/3388440.3412485>](https://doi.org/10.1145/3388440.3412485).

The package also provides a comparative type analysis of difference in association across contingency tables to reveal the highest order of their differences. 

Their null test statistics all follow an asymptotically chi-squared null distribution. These options test for heterogeneous patterns that differ in either the first order (marginal) or the second order (joint distribution deviation from product of marginals). Second-order differences may reveal more fundamental changes than first-order differences across heterogeneous patterns.

### When to use the package

This package takes a model-free approach without assuming an underlying parametric model for the relationship between variables, in contrast to differential correlation based on differences between linear models. Its input is contingency tables that store the counts or frequencies of discrete variables. Thus, continuous variables need to be discretized before using the tests. One option to do discretization is via optimal univariate clustering provided by the ['Ckmeans.1d.dp'](https://cran.r-project.org/package=Ckmeans.1d.dp) R package.

### To download and install the package

```{r}
install.packages("DiffXTables")
```
