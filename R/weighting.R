#' Add estimated propensity score to a data frame
#'
#' @description
#' adds propensity score to a data frame
#' @param data A data frame to be used.
#' @param object A \code{propmod} object if already fitted.
#' @param formula If not, write a \link[stats]{formula} to be fitted. Remember that you don't have to worry about group variable. \link[data.table]{.SD} do exclude `by`.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param var The name of the propensity score column.
#' @param mc Indicator column name for MC simulation if exists.
#' @param ... Additional arguments of fitting functions
#' @import data.table
#' @export
add_propensity <- function(data, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), var = "propensity", mc = NULL, ...) {
  method <- match.arg(method)
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  if (!is.null(object)) {
    data[[var]] <- estimate_ps(object, ...)
    data
    # data[,
    #      (var) := estimate_ps(object, ...)]
    # data[]
  } else {
    if (method == "logit") {
      # data[[var]] <-
      #   data %>%
      #   ps_glm(formula, data = ., ...) %>%
      #   estimate_ps()
      # data
      data[,
           (var) := ps_glm(formula, data = .SD, ...) %>%
             estimate_ps(),
           by = mc]
      data[]
    } else if (method == "rf") {
      # data[[var]] <-
      #   data %>%
      #   ps_rf(formula, data = ., ...) %>%
      #   estimate_ps()
      # data
      data[,
           (var) := ps_rf(formula, data = .SD, ...) %>%
             estimate_ps(),
           by = mc]
      data[]
    } else if (method == "cart") {
      # data[[var]] <-
      #   data %>%
      #   ps_cart(formula, data = ., ...) %>%
      #   estimate_ps()
      # data
      data[,
           (var) := ps_cart(formula, data = .SD, ...) %>%
             estimate_ps(),
           by = mc]
      data[]
    } else if (method == "SVM") {
      # data[[var]] <-
      #   data %>%
      #   ps_svm(formula, data = ., ...) %>%
      #   estimate_ps()
      # data
      data[,
           (var) := ps_svm(formula, data = .SD, ...) %>%
             estimate_ps(),
           by = mc]
      data[]
    }
  }
}

# weighting------------------------------------

add_ipw_wt <- function(data, treatment, trt_indicator = 1, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc = NULL, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  data %>%
    add_propensity(object = object, formula = formula, method = method, mc = mc, ...) %>%
    .[,
      treatment := ifelse(get(treatment) == trt_indicator, 1, 0)] %>%
    .[,
      ipw_wt := treatment / propensity - (1 - treatment) / (1 - propensity)] %>%
    .[]
  # data %>%
  #   add_propensity(object = object, formula = formula, method = method, ...) %>%
  #   mutate(treatment = ifelse(!!sym(treatment) == trt_indicator, 1, 0)) %>%
  #   mutate(ipw_wt = treatment / propensity - (1 - treatment) / (1 - propensity))
}

#' Estimation of Inverse Probability Weighting
#'
#' @description
#' estimates inverse probability weighting (IPW) based on propensity score estimates
#' @param data A data frame to be used
#' @param treatment Treatment variable name
#' @param trt_indicator Value that indicates the unit is treated
#' @param outcome Outcome variable name
#' @param object A \code{propmod} object if already fitted.
#' @param formula If not, write a \link[stats]{formula} to be fitted. Remember that you don't have to worry about group variable. \link[data.table]{.SD} do exclude `by`.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param mc Indicator column name for group if exists, e.g. c("mcname", "scenario")
#' \itemize{
#'  \item MC simulation
#'  \item Scenario
#' }
#' @param ... Additional arguments of fitting functions
#' @import data.table
#' @export
compute_ipw <- function(data, treatment, trt_indicator = 1, outcome, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc = NULL, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  data %>%
    add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, mc = mc, ...) %>%
    .[,
      .(IPW = mean(ipw_wt * get(outcome))),
      by = mc]
  # data %>%
  #   add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, ...) %>%
  #   summarise(
  #     IPW = mean(ipw_wt * !!sym(outcome))
  #   )
}

#' Estimation of Stabilized Inverse Probability Weighting
#'
#' @description
#' estimates stabilized inverse probability weighting (SIPW) based on propensity score estimates
#' @param data A data frame to be used
#' @param treatment Treatment variable name
#' @param trt_indicator Value that indicates the unit is treated
#' @param outcome Outcome variable name
#' @param object A \code{propmod} object if already fitted.
#' @param formula If not, write a \link[stats]{formula} to be fitted. Remember that you don't have to worry about group variable. \link[data.table]{.SD} do exclude `by`.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param mc Indicator column name for group if exists, e.g. c("mcname", "scenario")
#' \itemize{
#'  \item MC simulation
#'  \item Scenario
#' }
#' @param ... Additional arguments of fitting functions
#' @import data.table
#' @export
compute_sipw <- function(data, treatment, trt_indicator = 1, outcome, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc = NULL, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  data %>%
    add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, mc = mc, ...) %>%
    .[,
      sipw_wt := ipw_wt / sum(ipw_wt),
      by = treatment] %>%
    .[,
      .(
        SIPW = sum(treatment * sipw_wt * get(outcome) - (1 - treatment) * sipw_wt * get(outcome))
      ),
      by = mc]
  # data %>%
  #   add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, ...) %>%
  #   group_by(treatment) %>%
  #   mutate(sipw_wt = ipw_wt / sum(ipw_wt)) %>%
  #   ungroup() %>%
  #   summarise(
  #     SIPW = sum(
  #       treatment * sipw_wt * !!sym(outcome) - (1 - treatment) * sipw_wt * !!sym(outcome)
  #     )
  #   )
}

#' Inverse Probability Treatment Weighting
#'
#' @description
#' fits weighted regression of outcome on treatment
#' @param data A data frame to be used
#' @param treatment Treatment variable name
#' @param trt_indicator Value that indicates the unit is treated
#' @param object A \code{propmod} object if already fitted.
#' @param formula If not, write a \link[stats]{formula} to be fitted. Remember that you don't have to worry about group variable. \link[data.table]{.SD} do exclude `by`.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param mc Indicator column name for MC simulation if exists.
#' @param ... Additional arguments of fitting functions
#' @details
#' This functions add a column by
#' \deqn{\frac{trt_i}{\hat{e}_i} + \frac{1- trt_i}{1 - \hat{e}_i}}
#' @references Pirracchio, R., Petersen, M. L., & Laan, M. van der. (2015). \emph{Improving Propensity Score Estimators’ Robustness to Model Misspecification Using Super Learner}. American Journal of Epidemiology, 181(2), 108–119. \url{https://doi.org/10.1093/aje/kwu253}
#' @seealso
#' \code{\link{add_propensity}}
#' @import data.table
#' @export
add_iptw <- function(data, treatment, trt_indicator = 1, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc = NULL, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  data %>%
    add_propensity(object = object, formula = formula, method = method, mc = mc, ...) %>%
    .[,
      treatment := ifelse(get(treatment) == trt_indicator, 1, 0)] %>%
    .[,
      iptw := treatment / propensity + (1 - treatment) / (1 - propensity)] %>%
    .[,
      `:=`(treatment = NULL, propensity = NULL)] %>%
    .[]
}


