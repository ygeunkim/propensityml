#' Covariate Balance
#'
#' @description
#' gives covariate balance summary
#' @param data data
#' @param col_name column name of standardized difference means between treatment and control groups. By default, "\code{balance}"
#' @param treatment column name of treatment
#' @param trt_indicator value that indicates the unit is treated, e.g. 1 or TRUE
#' @param outcome outcome variable included in the data. It should be specified because it is not covariate.
#' @import data.table
#' @importFrom stats setNames var
#' @export
compute_balance <- function(data, col_name = "balance", treatment, trt_indicator = 1, outcome) {
  data %>%
    tidy_moment(treatment = treatment, col_exclude = outcome) %>%
    .[,
      .(balance = diff(value[moment == "mean"]) / sqrt(value[moment == "var" & get(treatment) == trt_indicator])),
      by = variable] %>%
    setNames(c("variable", col_name))
}

compute_moment <- function(x) {
  if (is.factor(x)) x <- as.numeric(levels(x))[x]
  list(mean = mean(x), var = var(x))
}

tidy_moment <- function(data, treatment, col_exclude) {
  data[,
       unlist(lapply(.SD, compute_moment)) %>% as.list(),
       by = treatment,
       .SDcols = -col_exclude] %>%
    melt(id.vars = treatment) %>%
    .[,
      c("variable", "moment") := tstrsplit(variable, ".", fixed = TRUE)]
}


#' Average Standardized Absolute Mean Distance
#'
#' @description
#' computes average standardized absolute mean distance (ASAM)
#' @param data data
#' @param treatment column name of treatment
#' @param trt_indicator value that indicates the unit is treated, e.g. 1 or TRUE
#' @param outcome outcome variable included in the data. It should be specified because it is not covariate.
#' @references Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
#' @details
#' For each covariate,
#' compute absolute (standardized difference of means between treatment and control groups), and take average.
#' Denote that standardization is done by sd of treatment group covariates.
#'
#' Lower ASAM means that treatment and control groups are more similar w.r.t. the given covariates.
#' @import data.table
#' @export
compute_asam <- function(data, treatment, trt_indicator = 1, outcome) {
  data %>%
    compute_balance(col_name = "balance", treatment = treatment, trt_indicator = trt_indicator, outcome = outcome) %>%
    .[, balance] %>%
    abs() %>%
    mean()
}
