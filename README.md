
# propensityml <a href='https://github.com/ygeunkim/propensityml'><img src='man/figures/logo.png' align="right" height="139" /></a>

[![Travis build
status](https://travis-ci.com/ygeunkim/propensityml.svg?branch=master)](https://travis-ci.com/ygeunkim/propensityml)
[![Codecov test
coverage](https://codecov.io/gh/ygeunkim/propensityml/branch/master/graph/badge.svg)](https://codecov.io/gh/ygeunkim/propensityml?branch=master)

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
#>       w1      w2 w3     w4 w5 w6      w7 w8 w9     w10 exposure outcome_prob
#>    1:  1 -0.9993  0  0.206  1  0 -0.4554  0  1  0.7953        1     1.07e-08
#>    2:  0 -0.2123  0 -0.937  1  0 -1.4709  0  0  0.1445        0     5.75e-01
#>    3:  1 -1.2956  1 -1.335  1  0 -0.1212  1  0  1.1384        1     1.00e+00
#>    4:  0 -0.7127  1 -0.785  1  0 -1.2592  0  0 -2.1676        0     1.00e+00
#>    5:  1  0.0527  1  0.463  1  0 -2.0587  1  0  0.4541        0     3.78e-01
#>   ---                                                                       
#>  996:  1  0.9268  1 -1.253  1  1  1.6927  1  0 -1.0818        1     2.27e-33
#>  997:  1 -0.3305  0  0.324  1  0 -1.1944  0  1 -0.3728        0     1.00e+00
#>  998:  0  0.0556  1  0.448  1  0 -1.0255  0  1  0.0524        0     4.66e-34
#>  999:  0 -0.4900  0  0.387  1  0 -0.3385  1  1 -1.2081        0     3.94e-32
#> 1000:  1  0.6380  0  1.559  0  1 -0.0819  1  1  1.0133        1     1.73e-51
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
#>         OOB estimate of  error rate: 53.9%
#> Confusion matrix:
#>     0   1 class.error
#> 0 292 233       0.444
#> 1 306 169       0.644
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
#> 0.377 0.464 0.283 0.473 0.543 0.567
```
