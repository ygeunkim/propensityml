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
#' @param mc_col Indicator column name for MC simulation if exists
#' @param sc_col Indicator column name for various scenarios if exists
#' @param parallel parallelize some operation
#' @param ... Additional arguments of fitting functions
#' @import data.table foreach
#' @export
add_propensity <- function(data, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), var = "propensity", mc_col = NULL, sc_col = NULL, parallel = FALSE, ...) {
  method <- match.arg(method)
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data)
    setDT(data)
  }
  if (!is.null(object)) {
    data[[var]] <- estimate_ps(object, ...)
    data
  }
  # for each method-------------------------------------
  switch(
    method,
    "logit" = {
      # glm------------------------------------
      if (!is.null(mc_col) & !is.null(sc_col)) {
        if (parallel) {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, get(mc_col)] %>% unique(), .combine = rbind) %dopar% {
              # sc_dt <- copy(data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)])
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)] %>%
                .[,
                  (var) := ps_glm(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        } else {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, .SD, .SDcols = -sc_col] %>% .[,get(mc_col)] %>% unique(), .combine = rbind) %do% {
              # sc_dt <- copy(data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)])
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)] %>%
                .[,
                  (var) := ps_glm(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        }
      } else {
        data[,
             (var) := ps_glm(formula, data = .SD, ...) %>%
               estimate_ps(),
             by = mc_col]
        data[]
      }
      #---------------------------------------
    },
    "rf" = {
      # randomForest---------------------------
      if (!is.null(mc_col) & !is.null(sc_col)) {
        if (parallel) {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, get(mc_col)] %>% unique(), .combine = rbind) %dopar% {
              # sc_dt <- copy(data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)])
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)] %>%
                .[,
                  (var) := ps_rf(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        } else {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, .SD, .SDcols = -sc_col] %>% .[,get(mc_col)] %>% unique(), .combine = rbind) %do% {
              # sc_dt <- copy(data[get(sc_col) == sc_id, .SD, .SDcols = -(sc_col, mc_col)])
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -sc_col] %>%
                .[,
                  (var) := ps_rf(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        }
      } else {
        data[,
             (var) := ps_rf(formula, data = .SD, ...) %>%
               estimate_ps(),
             by = mc_col]
        data[]
      }
      #---------------------------------------
    },
    "cart" = {
      # rpart--------------------------------------
      if (!is.null(mc_col) & !is.null(sc_col)) {
        if (parallel) {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, get(mc_col)] %>% unique(), .combine = rbind) %dopar% {
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)] %>%
                .[,
                  (var) := ps_cart(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        } else {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, .SD, .SDcols = -sc_col] %>% .[,get(mc_col)] %>% unique(), .combine = rbind) %do% {
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)] %>%
                .[,
                  (var) := ps_cart(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        }
      } else {
        data[,
             (var) := ps_cart(formula, data = .SD, ...) %>%
               estimate_ps(),
             by = mc_col]
        data[]
      }
      #---------------------------------------
    },
    "SVM" = {
      # svm---------------------------------------
      if (!is.null(mc_col) & !is.null(sc_col)) {
        if (parallel) {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, get(mc_col)] %>% unique(), .combine = rbind) %dopar% {
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)] %>%
                .[,
                  (var) := ps_svm(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        } else {
          foreach(sc_id = data[,get(sc_col)] %>% unique(), .combine = rbind) %:%
            foreach(mc_id = data[get(sc_col) == sc_id, .SD, .SDcols = -sc_col] %>% .[,get(mc_col)] %>% unique(), .combine = rbind) %do% {
              data[get(sc_col) == sc_id & get(mc_col) == mc_id, .SD, .SDcols = -c(sc_col, mc_col)] %>%
                .[,
                  (var) := ps_svm(formula, data = .SD, ...) %>%
                    estimate_ps()] %>%
                .[,
                  (mc_col) := mc_id] %>%
                .[,
                  (sc_col) := sc_id]
            }
        }
      } else {
        data[,
             (var) := ps_svm(formula, data = .SD, ...) %>%
               estimate_ps(),
             by = mc_col]
        data[]
      }
      #---------------------------------------
    }
  )
}

# weighting------------------------------------

add_ipw_wt <- function(data, treatment, trt_indicator = 1, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc_col = NULL, sc_col = NULL, parallel = FALSE, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  data %>%
    add_propensity(object = object, formula = formula, method = method, mc_col = mc_col, sc_col = sc_col, parallel = parallel, ...) %>%
    .[,
      treatment := ifelse(get(treatment) == trt_indicator, 1, 0)] %>%
    .[,
      ipw_wt := treatment / propensity - (1 - treatment) / (1 - propensity)] %>%
    .[]
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
#' @param mc_col Indicator column name for MC simulation if exists
#' @param sc_col Indicator column name for various scenarios if exists
#' @param parallel parallelize some operation
#' @param ... Additional arguments of fitting functions
#' @import data.table
#' @export
compute_ipw <- function(data, treatment, trt_indicator = 1, outcome, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc_col = NULL, sc_col = NULL, parallel = FALSE, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  # for single group--------------
  mc <- c(mc_col, sc_col)
  if (any(is.null(mc))) mc <- NULL
  #-------------------------------
  data %>%
    add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, mc_col = mc_col, sc_col = sc_col, parallel = parallel, ...) %>%
    .[,
      .(IPW = mean(ipw_wt * get(outcome))),
      by = mc]
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
#' @param mc_col Indicator column name for MC simulation if exists
#' @param sc_col Indicator column name for various scenarios if exists
#' @param parallel parallelize some operation
#' @param ... Additional arguments of fitting functions
#' @import data.table
#' @export
compute_sipw <- function(data, treatment, trt_indicator = 1, outcome, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc_col = NULL, sc_col = NULL, parallel = FALSE, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  # for single group--------------
  mc <- c(mc_col, sc_col)
  if (any(is.null(mc))) mc <- NULL
  #-------------------------------
  data %>%
    add_ipw_wt(treatment = treatment, trt_indicator = trt_indicator, object = object, formula = formula, method = method, mc_col = mc_col, sc_col = sc_col, parallel = parallel, ...) %>%
    .[,
      sipw_wt := ipw_wt / sum(ipw_wt),
      by = treatment] %>%
    .[,
      .(
        SIPW = sum(treatment * sipw_wt * get(outcome) - (1 - treatment) * sipw_wt * get(outcome))
      ),
      by = mc]
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
#' @param mc_col Indicator column name for MC simulation if exists
#' @param sc_col Indicator column name for various scenarios if exists
#' @param parallel parallelize some operation
#' @param ... Additional arguments of fitting functions
#' @details
#' This functions add a column by
#' \deqn{\frac{trt_i}{\hat{e}_i} + \frac{1- trt_i}{1 - \hat{e}_i}}
#' @references Pirracchio, R., Petersen, M. L., & Laan, M. van der. (2015). \emph{Improving Propensity Score Estimators’ Robustness to Model Misspecification Using Super Learner}. American Journal of Epidemiology, 181(2), 108–119. \url{https://doi.org/10.1093/aje/kwu253}
#' @seealso
#' \code{\link{add_propensity}}
#' @import data.table
#' @export
add_iptw <- function(data, treatment, trt_indicator = 1, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), mc_col = NULL, sc_col = NULL, parallel = FALSE, ...) {
  if (is.data.table(data)) {
    data <- copy(data)
  } else {
    data <- copy(data %>% data.table())
  }
  data %>%
    add_propensity(object = object, formula = formula, method = method, mc_col = mc_col, sc_col = sc_col, parallel = parallel, ...) %>%
    .[,
      treatment := ifelse(get(treatment) == trt_indicator, 1, 0)] %>%
    .[,
      iptw := treatment / propensity + (1 - treatment) / (1 - propensity)] %>%
    # .[,
    #   `:=`(treatment = NULL, propensity = NULL)] %>%
    .[,
      treatment := NULL] %>%
    .[]
}


