# COMBSS-VIGNETTE

Here we provide a first vignette to run COMBSS, a novel algorithm for best subset selection for a linear regression model, using R. For a vignette that runs COMBSS in Python, we refer to https://github.com/saratmoka/COMBSS-Python-VIGNETTE.

This vignette reproduces some replications from the simulation study presented in our article:

> Moka S, Liquet B, Zhu H, and Muller S (2022). COMBSS: Best Subset Selection via Continuous Optimization. *arXiv*, https://doi.org/10.48550/arxiv.2205.02617.


## Getting started

In this short vignette we use the following R packages

```
library(cPCG)
library(fields)
library(mvtnorm)
```

##  COMBSS in a low-dimensional setup

In this example, we consider a dataset with n = 100 samples and p = 20 predictors, of which 10 are active predictors.

This analysis is presented [here](/Low_dimensional_example.md).
 

## COMBSS in a high-dimensional setup

In this example, we consider a dataset with n = 100 samples and p = 1000 predictors, of which 10 are active predictors.

This analysis is presented [here](/High_dimensional_example.md).

## COMBSS in an ultra high-dimensional setup

In this example, we consider a dataset with n = 100 samples and p = 10,000 predictors, of which 3 are active predictors.

This analysis is presented [here](/Ultra_High_dimensional_example.md).
