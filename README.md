The 'DiffXTables' R package
===============================

The quesiton that whether the relationship between two variables has changed across conditions is often fundamental to a scientific inquiry. For example, a biologist could ask whether the relationship between two genes in a cancer cell has been modified from a normal cell. The 'DiffXTables' R package answers such questions via evaluating statistical evidence for distributional changes in the involved variables from data.

The package provides statistical methods for hypothesis testing of differences in the underlying distributions across two or more contingency tables. They include three statistical tests: 

1. The comparative chi-squared test: Its statistical foundation is first established in [(Song et al, 2014) <doi:10.1093/nar/gku086>](https://doi.org/10.1093/nar/gku086). It is later extended to identify differential patterns in networks [(Zhang et al, 2015) <doi:10.1093/nar/gkv358>](https://doi.org/10.1093/nar/gkv358).
2. The Sharma-Song test: A manuscript describing its theoretical foundation is being submitted for peer review.
3. The heterogeneity test: It is described in (Zar, 2010) and widely appears in textbooks. In contrast to the above comparative chi-squared test, it is not always powerful as demonstrated by examples in the package vignette.

They all have an asymptotically chi-squared null distribution. These options allow one to test for patterns that differ in first and second orders in their distributions, related to the marginal and joint distributions of the given contingency tables.

This package takes a model-free approach without assuming an underlying parametric model for the relationship between variables, in contrast to differential correlation based on linear models. Its input is contingency tables that store the counts or frequencies of discrete variables. Thus, continuous variables need to be discretized before using the tests. One option to do discretization is via optimal univariate clustering provided by the ['Ckmeans.1d.dp'](https://cran.r-project.org/package=Ckmeans.1d.dp) R package.


### To download and install the package

```{r}
install.packages("DiffXTables")
```
