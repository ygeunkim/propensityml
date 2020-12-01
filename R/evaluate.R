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
                         object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), weighting = c("IPW", "SIPW"), ...) {
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
        object = object, formula = formula, method = method, ...
      )
  } else { # sipw
    wt <-
      data %>%
      compute_sipw(
        treatment = treatment, trt_indicator = trt_indicator, outcome = outcome,
        object = object, formula = formula, method = method, ...
      )
  }
  data[,
       (weighting) := wt]
  # Balancing--------------------------------------------------------------
  wt_vars <- paste0(weighting, c(".mean", ".var"))
  left_formula <- paste0(c(treatment, "variable", wt_vars), collapse = "+")
  formul <- paste0(c(left_formula, "moment"), collapse = "~")
  formul <- as.formula(formul)
  # covariate column names--------------------------
  covariate_name <- names(data)
  covariate_name <- setdiff(covariate_name, c(outcome, treatment, exclude, weighting))
  # covariate columns that are factor---------------
  cols_fct <- sapply(data, class)[covariate_name]
  cols_fct <- names(cols_fct[cols_fct == "factor"])
  # change to numeric-------------------------------
  for (col in cols_fct) set(data, j = col, value = as.numeric(levels(data[[col]]))[data[[col]]])
  # weighing----------------------------------------
  for (col in covariate_name) set(data, j = col, value = data[[col]] * data[[weighting]])
  # balance dt--------------------------------------
  data %>%
    compute_balance(
      treatment = treatment,
      trt_indicator = trt_indicator,
      outcome = outcome,
      exclude = c(exclude, weighting)
    ) %>%
    .[,
      .(asam = mean(abs(balance)))]
}

#' ASAM by group
#'
#' @description
#' computes average standardized absolute mean distance (ASAM) in MC setting
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
#' @param mc_col Indicator column name for MC simulation if exists
#' @param sc_col Indicator column name for various scenarios if exists
#' @param parallel parallelize some operation
#' @param mc_core The number of cores to use for MC simulation
#' @param ... Additional arguments of fitting functions
#' @references Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
#' @details
#' For each covariate,
#' compute absolute (standardized difference of means between treatment and control groups), and take average.
#' Denote that standardization is done by sd of treatment group covariates.
#'
#' Lower ASAM means that treatment and control groups are more similar w.r.t. the given covariates.
#' @import data.table foreach
#' @importFrom parallel mclapply
#' @export
mc_asam <- function(data, treatment, trt_indicator = 1, outcome, exclude = NULL,
                    object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"),
                    weighting = c("IPW", "SIPW"), mc_col, sc_col = NULL, parallel = FALSE, mc_core = 1, ...) {
  if (missing(mc_col)) stop("Use this function only in MC simulation setting")
  mc_list <- data[, get(mc_col)] %>% unique()
  if (!is.null(sc_col)) {
    # Multiple scenarios------------------------------------------------------------------------
    sc_list <- data[, get(sc_col)] %>% unique()
    if (!parallel) {
      balance_dt <- foreach(sc_id = sc_list, .combine = rbind) %do% {
        sc_dt <- copy(data[get(sc_col) == sc_id, .SD, .SDcols = -sc_col])
        # Average ASAM for each MC-------------------
        mclapply(
          mc_list,
          function(id) {
            compute_asam(
              sc_dt[get(mc_col) == id],
              treatment = treatment, trt_indicator = trt_indicator, outcome = outcome, exclude = c(exclude, mc_col),
              object = object, formula = formula, method = method, weighting = weighting, ...
            )
          },
          mc.cores = mc_core
        ) %>%
          rbindlist() %>%
          .[,
            `:=`(
              MC = sc_dt[,get(mc_col)] %>% unique(),
              SC = sc_id
            )] %>%
          .[]
      }
    } else {
      # parallel with foreach-----------------------------------------------------
      balance_dt <- foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %dopar% {
        sc_dt <- copy(data[get(sc_col) == sc_id, .SD, .SDcols = -sc_col])
        mclapply(
          mc_list,
          function(id) {
            compute_asam(
              sc_dt[get(mc_col) == id],
              treatment = treatment, trt_indicator = trt_indicator, outcome = outcome, exclude = c(exclude, mc_col),
              object = object, formula = formula, method = method, weighting = weighting, ...
            )
          },
          mc.cores = mc_core
        ) %>%
          rbindlist() %>%
          .[,
            `:=`(
              MC = sc_dt[,get(mc_col)] %>% unique(),
              SC = sc_id
            )] %>%
          .[]
      }
    }
    # return the result of various scenario--------------
    return(
      balance_dt[,
                 .(asam = mean(asam)),
                 by = SC]
    )
  }
  # One scenario-------------------------------------------------------------------------------
  mclapply(
    data[,get(mc_col)] %>% unique(),
    function(id) {
      compute_asam(
        data[get(mc_col) == id],
        treatment = treatment, trt_indicator = trt_indicator, outcome = outcome, exclude = c(exclude, mc_col),
        object = object, formula = formula, method = method, weighting = weighting, ...
      )
    },
    mc.cores = mc_core
  ) %>%
    rbindlist() %>%
    .[,
      MC := data[,get(mc_col)] %>% unique()] %>%
    .[,
      .(asam = mean(asam))]
}

