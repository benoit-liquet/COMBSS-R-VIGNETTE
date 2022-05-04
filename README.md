# COMBSS-VIGNETTE

Here we provide a first vignette to run COMBSS, a novel algorithm for best subset selection for a linear regression model, using R. For a vignette that runs COMBSS in Python, we refer to https://github.com/saratmoka/COMBSS-Python-VIGNETTE.

We mainly use R code here. In addition, we show how to use the python version of COMBSS from R.

This vignette reproduces some replications from the simulation study presented in our article:

> Moka S, Liquet B, Zhu H, and Muller S (2022). COMBSS: Best Subset Selection via Continuous Optimization. *Submitted to arXiv*, 36 pages.


## Getting started

In this short vignette we use the following R packages

```
library(cPCG)
library(fields)
library(mvtnorm)
```

##  COMBSS in a low-dimensional setup

- We consider in this simulated example, a data set of n=100 samples and p=20 predictors, where 10 are active predictors.

- This analysis is presented [here](/Low_dimensional_example.md)
 

## COMBSS in a high-dimensional setup

- We consider in this simulated example, a data set of n=100 samples and p=1000 predictors, where 10 are active predictors.

- This analysis is presented [here](/High_dimensional_example.md)

## COMBSS in an ultra high-dimensional setup

- We consider in this simulated example, a data set of n=100 samples and p=10,000 predictors, where only 3 are active predictors.

- This analysis is presented [here](/Ultra_High_dimensional_example.md)
