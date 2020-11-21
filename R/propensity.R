#' Fitting Logistic Regression for Propensity Score
#'
#' @description
#' fits logistic regression model before estimating propensity scores
#' @param formula an object \link[stats]{formula} to be fitted. Response should be treatment.
#' @param data optional data frame.
#' @param ... For additional options of \link[stats]{glm}.
#' @return \link[propensityml]{propmod} class, a list with model and its name
#' \itemize{
#'  \item model - \link[stats]{glm} model
#'  \item name - "glm"
#'  \item data - dataset
#' }
#' @references Rosenbaum, P. R., & Rubin, D. B. (1983). \emph{The central role of the propensity score in observational studies for causal effects}. Biometrika, 70(1), 41-55. \url{https://doi.org/10.1093/biomet/70.1.41}
#' @details
#' Propensity score is
#' \deqn{e(X) = P(Z_i = 1 \mid X_i = x)},
#' which is the conditional probability of receiving treatment.
#' Naturally, logit model is the easiest way to estimate the score.
#' @export
ps_glm <- function(formula, data, ...) {
  result <-
    data %>%
    glm(formula, data = ., family = binomial, ...)
  res <- structure(list(
    model = result,
    name = "glm",
    data = data
  ))
  class(res) <- "propmod"
  result
  return(invisible(res))
}

#' Fitting Random Forests for Propensity Score
#'
#' @description
#' fits random forests before estimating propensity scores
#' @param formula an object \link[stats]{formula} to be fitted. Response should be treatment.
#' @param data optional data frame. REMEMBER that treatment should be `TRUE` or `1`.
#' @param ... For additional options of \link[randomForest]{randomForest}.
#' @return \link[propensityml]{propmod} class, a list with model and its name
#' \itemize{
#'  \item model - \link[randomForest]{randomForest}
#'  \item name - "rf"
#'  \item data - dataset
#' }
#' @references Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
#' @importFrom randomForest randomForest
#' @export
ps_rf <- function(formula, data, ...) {
  result <-
    data %>%
    randomForest(formula, data = ., ...)
  res <- structure(list(
    model = result,
    name = "rf",
    data = data
  ))
  class(res) <- "propmod"
  result
  return(invisible(res))
}

#' Fitting CART for Propensity Score
#'
#' @description
#' fits classification tree before estimating propensity scores
#' @param formula an object \link[stats]{formula} to be fitted. Response should be treatment.
#' @param data optional data frame.
#' @param ... For additional options of \link[rpart]{rpart}.
#' @return \link[rpart]{rpart}
#' @references Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
#' @importFrom rpart rpart
#' @export
ps_cart <- function(formula, data, ...) {
  result <-
    data %>%
    rpart(formula, data = ., method = "class", ...)
  res <- structure(list(
    model = result,
    name = "cart",
    data = data
  ))
  class(res) <- "propmod"
  result
  return(invisible(res))
}

# Estimate-------------------------------------

#' Estimation of Propensity Score
#'
#' @description
#' estimates propensity score based on the given model
#' @param object fitted \link[propensityml]{propmod} object
#' @param ... additional arguments for \link[stats]{predict}
#' @export
estimate_ps <- function(object, ...) {
  if (object$name == "glm") {
    return(predict(object$model, type = "response"))
  } else if (object$name == "rf") {
    trt_lev <-
      ifelse(
        which(c("1", "TRUE") %in% object$model$classes) == 1,
        1,
        TRUE
      )
    pred <- predict(object$model, type = "prob")
    return(pred[, colnames(pred) == trt_lev])
  }
}

# Propensity score model class-----------------

#' `propmod` class
#'
#' @description
#' The \code{propmod} class ready for computing propensity score
#' @name propmod-class
#' @rdname propmod-class
#' @aliases propmod cvar-class
#' @importFrom methods setOldClass
#' @exportClass propmod
setOldClass("propmod")

