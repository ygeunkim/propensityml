% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propensity.R
\name{ps_rf}
\alias{ps_rf}
\title{Fitting Random Forests for Propensity Score}
\usage{
ps_rf(formula, data, ...)
}
\arguments{
\item{formula}{an object \link[stats]{formula} to be fitted. Response should be treatment.}

\item{data}{optional data frame. REMEMBER that treatment should be \code{TRUE} or \code{1}.}

\item{...}{For additional options of \link[randomForest]{randomForest}.}
}
\value{
\code{\link{propmod}} class, a list with model and its name
\itemize{
\item model - \link[randomForest]{randomForest}
\item name - "rf"
\item data - dataset
}
}
\description{
fits random forests before estimating propensity scores
}
\examples{
fit <- chemical \%>\% ps_rf(poisox ~ age + sex, data = .)
}
\references{
Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
}
