#' Covariate Balance
#'
#' @description
#' gives covariate balance summary
#' @param data data
#' @param col_name column name of standardized difference means between treatment and control groups. By default, "\code{balance}"
#' @param treatment column name of treatment
#' @param trt_indicator value that indicates the unit is treated, e.g. 1 or TRUE
#' @param outcome outcome variable included in the data. It should be specified because it is not covariate.
#' @param exclude Additional columns to exlude
#' @import data.table
#' @importFrom stringr str_remove_all
#' @importFrom stats setNames var
#' @export
compute_balance <- function(data, col_name = "balance", treatment, trt_indicator = 1, outcome, exclude = NULL) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  old <- names(data)
  setnames(data, old, str_remove_all(old, pattern = "\\."))
  data %>%
    tidy_moment(treatment = treatment, col_exclude = c(outcome, exclude)) %>%
    .[,
      .(balance = diff(value[moment == "mean"]) / ( sqrt(sum(value[moment == "var"]) / 2) ) ),
      by = variable] %>%
    # .[,
    #   .(balance = diff(value[moment == "mean"]) / sqrt(value[moment == "var" & get(treatment) == trt_indicator])),
    #   by = variable] %>%
    setNames(c("variable", col_name))
}

compute_moment <- function(x) {
  if (is.factor(x)) x <- as.numeric(levels(x))[x]
  list(mean = mean(x), var = var(x))
}

tidy_moment <- function(data, treatment, with_melt = NULL, col_exclude) {
  data[,
       unlist(lapply(.SD, compute_moment)) %>% as.list(),
       by = treatment,
       .SDcols = -col_exclude] %>%
    melt(id.vars = c(treatment, with_melt)) %>%
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
#' @param exclude Additional columns to exlude
#' @param object A \code{propmod} object if already fitted.
#' @param formula If not, write a \link[stats]{formula} to be fitted. Remember that you don't have to worry about group variable. \link[data.table]{.SD} do exclude `by`.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param weighting Weighting methods, IPW or SIPW
#' @param mc Indicator column name for group if exists, e.g. c("mcname", "scenario")
#' \itemize{
#'  \item MC simulation
#'  \item Scenario
#' }
#' @param parallel parallelize some operation
#' @param ... Additional arguments of fitting functions
#' @references Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
#' @details
#' For each covariate,
#' compute absolute (standardized difference of means between treatment and control groups), and take average.
#' Denote that standardization is done by sd of treatment group covariates.
#'
#' Lower ASAM means that treatment and control groups are more similar w.r.t. the given covariates.
#' @import data.table foreach
#' @export
compute_asam <- function(data, treatment, trt_indicator = 1, outcome, exclude = NULL,
                         object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"),
                         weighting = c("IPW", "SIPW"), mc = NULL, parallel = FALSE, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  # IPW or SIPW------------------------
  weighting <- match.arg(weighting)
  if (weighting == "IPW") { # ipw
    wt <-
      data %>%
      compute_ipw(
        treatment = treatment, trt_indicator = trt_indicator, outcome = outcome,
        object = object, formula = formula, method = method, mc = mc, ...
      )
  } else { # sipw
    wt <-
      data %>%
      compute_sipw(
        treatment = treatment, trt_indicator = trt_indicator, outcome = outcome,
        object = object, formula = formula, method = method, mc = mc, ...
      )
  }
  if (!is.null(mc)) {
    data <- merge(data, wt, by = mc)
  } else {
    data <- data[,
                 (weighting) := wt]
  }
  # Balancing--------------------------------------------------------------
  wt_vars <- paste0(weighting, c(".mean", ".var"))
  left_formula <- paste0(c(treatment, "variable", wt_vars), collapse = "+")
  formul <- paste0(c(left_formula, "moment"), collapse = "~")
  formul <- as.formula(formul)
  # covariate column names--------------------------
  covariate_name <- names(data)
  covariate_name <- setdiff(covariate_name, c(outcome, treatment, exclude, mc, weighting))
  # covariate columns that are factor---------------
  cols_fct <- sapply(data, class)[covariate_name]
  cols_fct <- names(cols_fct[cols_fct == "factor"])
  # change to numeric-------------------------------
  for (col in cols_fct) set(data, j = col, value = as.numeric(levels(data[[col]]))[data[[col]]])
  # weighing----------------------------------------
  for (col in covariate_name) set(data, j = col, value = data[[col]] * data[[weighting]])
  # balance dt--------------------------------------
  if (!is.null(mc)) {
    if (parallel) {
      balance_dt <- foreach(mc_id = data[,get(mc)] %>% unique(), .combine = rbind) %dopar% {
        data[get(mc) == mc_id] %>%
          .[, .SD, .SDcols = -mc] %>%
          compute_balance(
            treatment = treatment,
            trt_indicator = trt_indicator,
            outcome = outcome,
            exclude = c(exclude, weighting)
          ) %>%
          .[,
            .(asam = mean(balance))] %>%
          .[,
            MC := mc_id]
      }
    } else {
      balance_dt <- foreach(mc_id = data[,get(mc)] %>% unique(), .combine = rbind) %do% {
        data[get(mc) == mc_id] %>%
          .[, .SD, .SDcols = -mc] %>%
          compute_balance(
            treatment = treatment,
            trt_indicator = trt_indicator,
            outcome = outcome,
            exclude = c(exclude, weighting)
          ) %>%
          .[,
            .(asam = mean(balance))] %>%
          .[,
            MC := mc_id]
      }
    }
    return(balance_dt)
  }
  balance_dt <-
    data %>%
    compute_balance(
      treatment = treatment,
      trt_indicator = trt_indicator,
      outcome = outcome,
      exclude = c(exclude, weighting)
    )
  balance_dt[, .(asam = mean(balance))]
}
