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


# Setoguchi paper-------------------------------

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
#' @export
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
    sqrt(S) %*% x %*% sqrt(S) %>% as.matrix()
  } else {
    x %>% as.matrix()
  }
}

#' Simulating covariates
#'
#' @description
#' generates covariates in the references
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
sim_covariate <- function(n, covmat = build_covariate()) {
  x <-
    rmvnorm(n, sigma = covmat) %>%
    data.table()
  new_col <- paste0("w", 1:10)
  setnames(
    x,
    old = paste0("V", 1:10),
    new = new_col
  )
  # dichotomize w1, w3, w5, w6, w8, w9------------------
  binary <- paste0("w", c(1, 3, 5, 6, 8, 9))
  x[,
    (binary) := lapply(.SD, function(x) {
      ifelse(x > mean(x), 1, 0)
    }),
    .SDcols = binary] %>%
    .[,
      w0 := 1]
  setcolorder(x, c("w0", new_col))
  x[]
}

#' Simulating Dataset for Various Scenarios
#'
#' @description
#' generates a dataset for various scenarios
#' @param n sample size
#' @param covmat Covariance matrix of the covariates
#' @param scenario scenarios
#' @param b coefficients for confounder and exposure predictors
#' @param a coefficients in outcome model
#' @param gam coefficient of exposure
#' @references Setoguchi, S., Schneeweiss, S., Brookhart, M. A., Glynn, R. J., & Cook, E. F. (2008). \emph{Evaluating uses of data mining techniques in propensity score estimation: a simulation study}. Pharmacoepidemiology and Drug Safety, 17(6), 546–555 \url{https://doi.org/10.1002/pds.1555}
#' @references Lee, B. K., Lessler, J., & Stuart, E. A. (2010). \emph{Improving propensity score weighting using machine learning. Statistics in Medicine}. Statistics in Medicine, 29(3), 337-346. \url{https://doi.org/10.1002/sim.3782}
#' @details
#' About scenarios:
#' \itemize{
#'  \item A: additivity and linearity
#'  \item B: mild non-linearity
#'  \item C: moderate non-linearity
#'  \item D: mild non-additivity
#'  \item E: mild non-additivity and non-linearity
#'  \item F: moderate non-linearity
#'  \item F: moderate non-additivity and non-linearity
#' }
#' See Appendix of Setoguchi et al.
#' @import data.table
#' @export
sim_outcome <- function(n, covmat = build_covariate(), scenario = LETTERS[1:7],
                        b = c(0, .8, -.25, .6, -.4, -.8, -.5, .7),
                        a = c(-3.85, .3, -.36, -73, -.2, .71, -.19, .26),
                        gam = -.4) {
  scenario <- match.arg(scenario)
  x <- sim_covariate(n, covmat)
  covariate <- paste0("w", 0:7)
  confounder <- paste0("w", 1:4)
  out_cov <- paste0("w", 8:10)
  if (scenario == "A") {
    x[,
      exposure_prob := Reduce("+", b * .SD[, .SD, .SDcols = covariate])]
  } else if (scenario == "B") {
    b <- c(b, b[2])
    x[,
      w2w2 := w2 * w2] %>%
      .[,
        exposure_prob := Reduce("+", b * .SD[, .SD, .SDcols = c(covariate, "w2w2")])]
  } else if (scenario == "C") {
    b <- c(b, b[2], b[4], b[7])
    x[,
      `:=` (
        w2w2 = w2 * w2,
        w4w4 = w4 * w4,
        w7w7 = w7 * w7
      )] %>%
      .[,
        exposure_prob := Reduce("+", b * .SD[, .SD, .SDcols = c(covariate, "w2w2", "w4w4", "w7w7")])]
  } else if (scenario == "D") {
    b <- c(b, .5 * b[1], .7 * b[2], .5 * b[4], .5 * b[5])
    x[,
      `:=` (
        w1w3 = w1 * w3,
        w2w4 = w2 * w4,
        w4w5 = w4 * w5,
        w5w6 = w5 * w6
      )] %>%
      .[,
        exposure_prob := Reduce("+", b * .SD[, .SD, .SDcols = c(covariate, "w1w3", "w2w4", "w4w5", "w5w6")])]
  } else if (scenario == "E") {
    b <- c(b, b[2], .5 * b[1], .7 * b[2], .5 * b[4], .5 * b[5])
    x[,
      `:=` (
        w2w2 = w2 * w2,
        w1w3 = w1 * w3,
        w2w4 = w2 * w4,
        w4w5 = w4 * w5,
        w5w6 = w5 * w6
      )] %>%
      .[,
        exposure_prob := Reduce("+", b * .SD[, .SD, .SDcols = c(covariate, "w2w2", "w1w3", "w2w4", "w4w5", "w5w6")])]
  } else if (scenario == "F") {
    b <- c(b, .5 * b[1], .7 * b[2], .5 * b[3], .7 * b[4], .5 * b[5], .5 * b[1], .7 * b[2], .5 * b[3], .5 * b[4], .5 * b[5])
    x[,
      `:=` (
        w1w3 = w1 * w3,
        w2w4 = w2 * w4,
        w3w5 = w3 * w5,
        w4w6 = w4 * w6,
        w5w7 = w5 * w7,
        w1w6 = w1 * w6,
        w2w3 = w2 * w3,
        w3w4 = w3 * w4,
        w4w5 = w4 * w5,
        w5w6 = w5 * w6
      )] %>%
      .[,
        exposure_prob := Reduce("+", b * .SD[, .SD, .SDcols = c(covariate, "w1w3", "w2w4", "w3w5", "w4w6", "w5w7", "w1w6", "w2w3", "w3w4", "w4w5", "w5w6")])]
  } else if (scenario == "G") {
    b <- c(
      b,
      b[2], b[4], b[7], # C
      .5 * b[1], .7 * b[2], .5 * b[4], .5 * b[5], # D
      .5 * b[3], .7 * b[4], .5 * b[5], .5 * b[1], .7 * b[2], .5 * b[3] # F
    )
    x[,
      `:=` (
        # C-----------
        w2w2 = w2 * w2,
        w4w4 = w4 * w4,
        w7w7 = w7 * w7,
        # D-----------
        w1w3 = w1 * w3,
        w2w4 = w2 * w4,
        w4w5 = w4 * w5,
        w5w6 = w5 * w6,
        # F-----------
        w3w5 = w3 * w5,
        w4w6 = w4 * w6,
        w5w7 = w5 * w7,
        w1w6 = w1 * w6,
        w2w3 = w2 * w3,
        w3w4 = w3 * w4
      )] %>%
      .[,
        exposure_prob := Reduce("+", b * .SD[, .SD, .SDcols = c(covariate, "w2w2", "w4w4", "w7w7", "w1w3", "w2w4", "w4w5", "w5w6", "w3w5", "w4w6", "w5w7", "w1w6", "w2w3", "w3w4")])]
  }
  x[,
    exposure_prob := (1 + exp(-exposure_prob))^(-1)] %>%
    .[,
      `:=`(
        exposure = ifelse(exposure_prob > runif(n), 1, 0),
        exposure_prob = NULL
      )] %>%
    .[,
      outcome_prob := Reduce("+", c(a, gam) * .SD[, .SD, .SDcols = c(confounder, out_cov, "exposure")])] %>%
    .[,
      outcome_prob := (1 + exp(-outcome_prob))^(-1)]
  bin_cols <- c("exposure", paste0("w", c(1, 3, 5, 6, 8, 9)))
  x[,
    (bin_cols) := lapply(.SD, factor),
    .SDcols = bin_cols]
  x[, .SD, .SDcols = c(paste0("w", 1:10), "exposure", "outcome_prob")]
}
