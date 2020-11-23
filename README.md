
# propensityml <a href='https://github.com/ygeunkim/propensityml'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/ygeunkim/propensityml.svg?branch=master)](https://travis-ci.com/ygeunkim/propensityml)
[![Codecov test
coverage](https://codecov.io/gh/ygeunkim/propensityml/branch/master/graph/badge.svg)](https://codecov.io/gh/ygeunkim/propensityml?branch=master)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

## Overview

This is an R package to help the [SKKU modern statistical methods
project](https://github.com/ygeunkim/psweighting-ml). It is basically
based on the paper

[Lee, B. K., Lessler, J., & Stuart, E. A. (2010). Improving propensity
score weighting using machine learning. Statistics in
Medicine, 29(3), 337–346.
doi:10.1002/sim.3782](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.3782)

## Installation

``` r
# install.packages("remotes")
remotes::install_github("ygeunkim/propensityml")
```

## Usage

`propensityml` package aims at estimating propensity score with machine
learning methods as in the paper mentioned above.

``` r
library(propensityml)
```

The package provides simulation function that generates the dataset in
the paper:

[Setoguchi, S., Schneeweiss, S., Brookhart, M. A., Glynn, R. J., & Cook,
E. F. (2008). Evaluating uses of data mining techniques in propensity
score estimation: a simulation study. Pharmacoepidemiology and Drug
Safety, 17(6), 546–555
https://doi.org/10.1002/pds.1555](https://onlinelibrary.wiley.com/doi/abs/10.1002/pds.1555)

and additional toy datasets. Consider simulation.

The most simplest scenario, i.e. additivity and linearity model:

``` r
(x <- sim_outcome(1000, covmat = build_covariate()))
#>       w1      w2 w3      w4 w5 w6      w7 w8 w9     w10 exposure outcome_prob
#>    1:  1  0.0298  0  1.6870  0  1  0.5546  0  1  0.0615        0     1.42e-55
#>    2:  0  2.4582  1  0.9966  0  1 -1.0959  1  1  1.7027        1     1.17e-32
#>    3:  0  0.3381  0 -0.0133  0  1 -0.6201  0  0  0.5083        0     1.55e-11
#>    4:  0 -0.1344  1 -0.4582  0  1 -0.6519  1  0  0.7542        0     1.39e-01
#>    5:  0 -0.3925  0  1.5564  0  0 -0.2308  0  1 -0.1060        0     2.28e-02
#>   ---                                                                        
#>  996:  0  1.5070  1 -1.0451  1  1 -0.7617  0  0 -0.0549        1     2.99e-32
#>  997:  1  0.5003  0  1.5290  1  1  1.1412  0  1  0.0915        0     2.29e-06
#>  998:  1 -0.4310  1 -0.1737  1  0  1.7049  1  0  0.3598        0     1.98e-02
#>  999:  0 -0.1910  0  0.3357  1  0 -1.5107  1  0 -0.2657        1     2.51e-32
#> 1000:  1 -0.4443  1  0.8992  1  0  0.0272  0  1 -0.5417        1     1.18e-30
```

``` r
(fit_rf <- 
  x %>% 
  ps_rf(exposure ~ . -outcome_prob, data = .))
#> 
#> Call:
#>  randomForest(formula = formula, data = .) 
#>                Type of random forest: classification
#>                      Number of trees: 500
#> No. of variables tried at each split: 3
#> 
#>         OOB estimate of  error rate: 48.9%
#> Confusion matrix:
#>     0   1 class.error
#> 0 268 235       0.467
#> 1 254 243       0.511
```

We have defined the class named `propmod` for some usage.

``` r
class(fit_rf)
#> [1] "propmod"
```

Estimating propensity score:

``` r
estimate_ps(fit_rf) %>% head()
#>     1     2     3     4     5     6 
#> 0.391 0.661 0.469 0.359 0.544 0.472
```
