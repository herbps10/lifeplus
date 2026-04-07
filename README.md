# lifeplus
<!-- badges: start -->
[![R-CMD-check](https://github.com/herbps10/lifeplus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/herbps10/lifeplus/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview
This package implements a modular Bayesian modeling framework for estimating and projecting life expectancy in multiple populations.

## Installation

You can install the development version of `lifeplus` from [Github](https://github.com/herbps10/lifeplus) with:

```r
# install.packages("pak")
pak::pak("herbps10/lifeplus")
```

## Model Framework
Let $c = 1, \dots, C$ index the areas (typically countries) and $t = 1, \dots, T$ index time points (typically years). Let $y_{c,t}$ be the observed life expectancy in area $c$ at time $t$. The model takes the form

$$
\mathrm{logit}\left(y_{c,t}\right) = g\left(y_{c,t-1}, \beta_c\right) + \delta_{c,t} + \epsilon_{c,t}
$$

where

- $g$ is the _transition model_ with area-specific parameters $\beta_c$.
- $\delta_{c,t}$ is the _shock model_.
- $\epsilon_{c,t}$ is a residual modeled by a _data model_.

Multiple options are available for each of the models.

### Transition functions

- **Double logistic**: `transition_model_double_logistic`. The transition is modeled by a parametric double logistic function (Raftery et al. 2013). 

### Shock terms

- **None**: `shock_model_none`. All shock terms are fixed to zero: $\delta_{c,t} = 0$.
- **Regularized horseshoe**: `shock_model_regularized_horseshoe`. Shock terms are modeled with a regularized horseshoe prior (Piironen and Vehtari 2017). 

### Data models
- *Normal*: `data_model_normal`. The residuals are modeled as $\epsilon_{c,t} \sim N\left(0, \sigma^2\right)$.
- *Outlier*:  `data_model_outlier` . Any absolute life expectancy differences larger than a user-specified threshold are ignored. The remaining observations are modeled as $\epsilon_{c,t} \sim N(\left(0, \sigma^2\right)$. 


## References

- Raftery, A.E., Chunn, J.L., Gerland, P. et al. [Bayesian probabilistic projections of life expectancy for all countries](https://doi.org/10.1007/s13524-012-0193-x). 
  Demography 50, 777–801 (2013).
- Juho Piironen, Aki Vehtari. [Sparsity information and regularization in the horseshoe and other shrinkage priors](https://doi.org/10.1214/17-EJS1337SI).
  Electronic Journal of Statistics, 11(2), 5018-5051, (2017)  
