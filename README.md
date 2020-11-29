
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
#>       w1      w2 w3     w4 w5 w6      w7 w8 w9    w10 exposure       y
#>    1:  0 -0.0234  1  1.397  0  0 -0.2294  0  1 -0.865        0 -101.50
#>    2:  0 -0.8632  0 -0.144  1  0  1.1669  1  0  2.907        1    1.41
#>    3:  1 -2.4124  1 -0.224  0  0 -1.7911  1  0  1.695        0  174.52
#>    4:  1 -0.7639  0 -0.838  1  0 -1.1463  0  0 -1.809        0  -65.72
#>    5:  1  1.1810  0  0.352  0  1  0.8446  0  1 -0.344        0   -3.22
#>   ---                                                                 
#>  996:  0  0.9541  0  1.736  1  1  2.3953  0  1 -0.510        1  -73.39
#>  997:  1 -1.0332  0  0.703  0  0 -0.1423  0  1 -0.143        1    6.69
#>  998:  1 -2.3697  1  0.068  1  0  0.4174  0  0  0.291        0   -2.68
#>  999:  1  0.8074  1 -0.834  0  1 -0.7798  0  0 -0.430        0   -3.21
#> 1000:  1 -0.3249  0  1.176  0  1 -0.0487  0  1 -0.394        1  -88.73
#>       exposure_prob
#>    1:         0.500
#>    2:         0.849
#>    3:         0.588
#>    4:         0.538
#>    5:         0.104
#>   ---              
#>  996:         0.992
#>  997:         0.352
#>  998:         0.155
#>  999:         0.168
#> 1000:         0.935
```

``` r
(fit_rf <- 
  x %>% 
  ps_rf(exposure ~ . - y - exposure_prob, data = .))
#> 
#> Call:
#>  randomForest(formula = formula, data = .) 
#>                Type of random forest: classification
#>                      Number of trees: 500
#> No. of variables tried at each split: 3
#> 
#>         OOB estimate of  error rate: 47.2%
#> Confusion matrix:
#>     0   1 class.error
#> 0 266 232       0.466
#> 1 240 262       0.478
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
#> 0.433 0.584 0.500 0.548 0.446 0.364
```
