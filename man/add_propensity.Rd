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
  mc_col = NULL,
  sc_col = NULL,
  parallel = FALSE,
  ...
)
}
\arguments{
\item{data}{A data frame to be used.}

\item{object}{A \code{propmod} object if already fitted.}

\item{formula}{If not, write a \link[stats]{formula} to be fitted. Remember that you don't have to worry about group variable. \link[data.table]{.SD} do exclude \code{by}.}

\item{method}{Estimating methods
\itemize{
\item "logit" - \code{\link{ps_glm}}
\item "rf" - \code{\link{ps_rf}}
\item "cart" - \code{\link{ps_cart}}
\item "SVM" - \code{\link{ps_svm}}
}}

\item{var}{The name of the propensity score column.}

\item{mc_col}{Indicator column name for MC simulation if exists}

\item{sc_col}{Indicator column name for various scenarios if exists}

\item{parallel}{parallelize some operation}

\item{...}{Additional arguments of fitting functions}
}
\description{
adds propensity score to a data frame
}
