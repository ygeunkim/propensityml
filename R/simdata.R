#' Potential Outcome Framework Toy Dataset
#'
#' @description
#' Toy example for potential outcome framework. The table construction is given by:
#'
#' @docType data
#'
#' @format A data frame with 8 rows and 5 variables:
#' \describe{
#'   \item{y0}{outcome when the unit is untreated}
#'   \item{y1}{outcome when the unit is treated}
#'   \item{z}{binary treatment}
#'   \item{y}{outcome}
#'   \item{n}{observed count}
#' }
"outcome_frame"

#' Cohort Study Data about POISOX
#'
#' @description
#' Cohort study interesting in the effect of exposure to an industrial chemical, POISOX, on subsequent blood pressure and mortality.
#'
#' @docType data
#'
#' @format A data frame with 5000 rows and 5 variables:
#' \describe{
#'   \item{age}{age of each subject}
#'   \item{sex}{sex of each subject, 0 if male, 1 if female}
#'   \item{poisox}{0 if unexposed, 1 if exposed}
#'   \item{mortal}{0 if alive, 1 if dead}
#'   \item{blood}{subsequent blood pressure, after - prior}
#' }
"chemical"

#' Simulation dataset
#'
#' @description
#' generates dataset with various scenarios
#' @param n sample size
#' @param covmat Covariance matrix of the covariates
#' @references Setoguchi, S., Schneeweiss, S., Brookhart, M. A., Glynn, R. J., & Cook, E. F. (2008). \emph{Evaluating uses of data mining techniques in propensity score estimation: a simulation study}. Pharmacoepidemiology and Drug Safety, 17(6), 546–555 \url{https://doi.org/10.1002/pds.1555}
#' @references Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
#' @details
#' This function reproduces the setting in the paper in Setoguchi et al. and Lee et al.
#' First generate binary covariates (w1, w3, w5, w6, w8, w9), and continuous covariates (w2, w4, w7, w10).
#' \enumerate{
#' \item Generate 10-dim multivariate normal v1, v3, v5, v6, v8, v9 (corresponding to binary), and w7, w10
#' }
#' @seealso \code{\link{build_covariate}}
#' @importFrom mvtnorm rmvnorm
#' @import data.table
#' @export
sim_outcome <- function(n, covmat = build_covariate()) {
  x <-
    rmvnorm(n, sigma = covmat %>% as.matrix()) %>%
    data.table()
  setnames(
    x,
    old = paste0("V", 1:10),
    new = paste0("w", 1:10)
  )
  x
  #------------------------
}

#' Build covariance matrix of the covariates
#'
#' @description
#' builds covariance matrix of the onfounders, exposure predictors, and outcome predictors
#' @param sig variance of each covariate
#' @param confounder_exposure covariance between confounder and exposure predictors, vector in order of column
#' @param confounder_outcome covariance between confounder and outcome predictors, vector in order of column
#' @param exposure_outcome covariance between exposure predictors and outcome predictors, vector in order of column
#' @param value_cor are given values (except `sig`) are correlation? (default = TRUE)
#' @references Setoguchi, S., Schneeweiss, S., Brookhart, M. A., Glynn, R. J., & Cook, E. F. (2008). \emph{Evaluating uses of data mining techniques in propensity score estimation: a simulation study}. Pharmacoepidemiology and Drug Safety, 17(6), 546–555 \url{https://doi.org/10.1002/pds.1555}
#' @details
#' This function builds covariance matrix of the covariates, which are
#' confounders (w1, w2, w3, w4), exposure predictors (w5, w6, w7), and outcome predictors(w8, w9, w10).
#' In each group, there is no correlation.
#' On the other hand, you can specify picewise correlation between e.g. confounder-exposure predictor.
#' @seealso \code{\link{sim_outcome}}
#' @importFrom Matrix bdiag symmpart
build_covariate <- function(sig = rep(1L, 10),
                            confounder_exposure = c(.2, 0, 0, 0, .9, 0, 0, 0, 0, 0, 0, 0),
                            confounder_outcome = c(0, 0, 0, 0, 0, 0, .2, 0, 0, 0, .9, 0),
                            exposure_outcome = rep(0L, 9),
                            value_cor = TRUE) {
  confounder <- diag(sig[1:4])
  exposure <- diag(sig[5:7])
  outcome <- diag(sig[8:10])
  S <- bdiag(confounder, exposure, outcome)
  x <- S
  x[5:7, 1:4] <- confounder_exposure * 2
  x[8:10, 1:4] <- confounder_outcome * 2
  x[8:10, 5:7] <- exposure_outcome * 2
  x <- symmpart(x)
  if (value_cor) {
    sqrt(S) %*% x %*% sqrt(S)
  } else {
    x
  }
}
