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
#' @export
add_propensity <- function(data, object = NULL, formula = NULL, method = c("logit", "rf", "cart", "SVM"), var = "propensity", ...) {
  method <- match.arg(method)
  if (!is.null(object)) {
    data[[var]] <- estimate_ps(object, ...)
    data
  } else {
    if (method == "logit") {
      data[[var]] <-
        data %>%
        ps_glm(formula, data = ., ...) %>%
        estimate_ps()
      data
    } else if (method == "rf") {
      data[[var]] <-
        data %>%
        ps_rf(formula, data = ., ...) %>%
        estimate_ps()
      data
    } else if (method == "cart") {
      data[[var]] <-
        data %>%
        ps_cart(formula, data = ., ...) %>%
        estimate_ps()
    } else if (method == "SVM") {
      data[[var]] <-
        data %>%
        ps_svm(formula, data = ., ...) %>%
        estimate_ps()
    }
  }
}

# weighting------------------------------------

