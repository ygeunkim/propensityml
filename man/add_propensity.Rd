% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighting.R
\name{add_propensity}
\alias{add_propensity}
\title{Add estimated propensity score to a data frame}
\usage{
add_propensity(
  data,
  object = NULL,
  formula = NULL,
  method = c("logit", "rf", "cart", "SVM"),
  var = "propensity",
  ...
)
}
\arguments{
\item{data}{A data frame to be used.}

\item{object}{A \code{propmod} object if already fitted.}

\item{formula}{If not, write a \link[stats]{formula} to be fitted.}

\item{method}{Estimating methods
\itemize{
\item "logit" - \code{\link{ps_glm}}
\item "rf" - \code{\link{ps_rf}}
\item "cart" - \code{\link{ps_cart}}
\item "SVM" - \code{\link{ps_svm}}
}}

\item{var}{The name of the propensity score column.}

\item{...}{Additional arguments of fitting functions}
}
\description{
adds propensity score to a data frame
}