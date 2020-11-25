#' Add estimated propensity score to a data frame
#'
#' @description
#' adds propensity score to a data frame
#' @param data A data frame to be used.
#' @param object A \code{propmod} object if already fitted.
#' @param formula If not, write a \link[stats]{formula} to be fitted.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param var The name of the propensity score column.
#' @param ... Additional arguments of fitting functions
#' @import data.table
#' @export
add_propensity <- function(data, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), var = "propensity", ...) {
  method <- match.arg(method)
  if (!is.null(object)) {
    # data[[var]] <- estimate_ps(object, ...)
    data[,
         (var) := estimate_ps(object, ...)]
    data[]
  } else {
    if (method == "logit") {
      # data[[var]] <-
      #   data %>%
      #   ps_glm(formula, data = ., ...) %>%
      #   estimate_ps()
      data[,
           (var) := ps_glm(formula, data = .SD, ...) %>%
             estimate_ps()]
      data[]
    } else if (method == "rf") {
      # data[[var]] <-
      #   data %>%
      #   ps_rf(formula, data = ., ...) %>%
      #   estimate_ps()
      data[,
           (var) := ps_glm(formula, data = .SD, ...) %>%
             estimate_ps()]
      data[]
    } else if (method == "cart") {
      # data[[var]] <-
      #   data %>%
      #   ps_cart(formula, data = ., ...) %>%
      #   estimate_ps()
      data[,
           (var) := ps_cart(formula, data = .SD, ...) %>%
             estimate_ps()]
      data[]
    } else if (method == "SVM") {
      # data[[var]] <-
      #   data %>%
      #   ps_svm(formula, data = ., ...) %>%
      #   estimate_ps()
      data[,
           (var) := ps_svm(formula, data = .SD, ...) %>%
             estimate_ps()]
      data[]
    }
  }
}

# weighting------------------------------------

#' Estimation of Inverse Probability Weighting
#'
#' @description
#' estimates inverse probability weighting (IPW) based on propensity score estimates
#' @param data A data frame to be used
#' @param treatment Treatment variable name
#' @param trt_indicator Value that indicates the unit is treated
#' @param outcome Outcome variable name
#' @param object A \code{propmod} object if already fitted.
#' @param formula If not, write a \link[stats]{formula} to be fitted.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param ... Additional arguments of fitting functions
#' @importFrom dplyr mutate summarise
#' @importFrom rlang sym
#' @export
compute_ipw <- function(data, treatment, trt_indicator = 1, outcome, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), ...) {
  data %>%
    add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, ...) %>%
    summarise(
      IPW = mean(ipw_wt * !!sym(outcome))
    )
}

add_ipw_wt <- function(data, treatment, trt_indicator = 1, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), ...) {
  data %>%
    add_propensity(object = object, formula = formula, method = method, ...) %>%
    mutate(treatment = ifelse(!!sym(treatment) == trt_indicator, 1, 0)) %>%
    mutate(ipw_wt = treatment / propensity - (1 - treatment) / (1 - propensity))
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
#' @param formula If not, write a \link[stats]{formula} to be fitted.
#' @param method Estimating methods
#' \itemize{
#'  \item "logit" - \code{\link{ps_glm}}
#'  \item "rf" - \code{\link{ps_rf}}
#'  \item "cart" - \code{\link{ps_cart}}
#'  \item "SVM" - \code{\link{ps_svm}}
#' }
#' @param ... Additional arguments of fitting functions
#' @importFrom dplyr mutate group_by ungroup summarise
#' @importFrom rlang sym
#' @export
compute_sipw <- function(data, treatment, trt_indicator = 1, outcome, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), ...) {
  data %>%
    add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, ...) %>%
    group_by(treatment) %>%
    mutate(sipw_wt = ipw_wt / sum(ipw_wt)) %>%
    ungroup() %>%
    summarise(
      SIPW = sum(
        treatment * sipw_wt * !!sym(outcome) - (1 - treatment) * sipw_wt * !!sym(outcome)
      )
    )
}






