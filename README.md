# lifeplus
<!-- badges: start -->
[![R-CMD-check](https://github.com/herbps10/lifeplus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/herbps10/lifeplus/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview
This package implements a modular Bayesian modeling framework for estimating and projecting life expectancy in multiple populations.

## Model Framework
Let $c = 1, \dots, C$ index the area and $t = 1, \dots, T$ index time points. Let $y_{c,t}$ be the observed life expectancy in area $c$ at time $t$. The model takes the form

$$
\mathrm{logit}(y_{c,t}) = g(y_{c,t-1}, \beta_c) + \delta_{c,t} + \epsilon_{c,t}
$$

where

- $g$ is the \textit{transition function} with area-specific parameters $\beta_c$
- $\delta_{c,t}$ is the \textit{shock term}
- $\epsilon_{c,t}$ is a residual modeled by a \textit{data model}.

### Transition functions

- *Double logistic* (`transition_modeldouble_logistic`).

### Shock terms

- *None* (`shock_model_none`). All shock terms are fixed to zero: $\delta_{c,t} = 0$.
- *Regularized horseshoe* (`shock_model_regularized_horseshoe`). 

### Data models
- *Normal* (`data_model_normal`). $\epsilon_{c,t} \sim N(0, \sigma^2)$.
- *Outlier* (`data_model_outlier`). 
