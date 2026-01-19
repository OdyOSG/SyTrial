#' Build Large-Scale Propensity Score Model
#'
#' @description
#' Constructs propensity score model using Cyclops for high-dimensional
#' covariate adjustment. Supports regularization for variable selection.
#'
#' @param covariateData CovariateData object from FeatureExtraction
#' @param treatmentVariable Binary treatment indicator (1 = treatment, 0 = control)
#' @param prior Prior specification for regularization:
#'   \itemize{
#'     \item "laplace": L1 regularization (LASSO)
#'     \item "normal": L2 regularization (Ridge)
#'   }
#' @param priorVariance Variance of the prior (controls regularization strength)
#' @param control Cyclops control object
#'
#' @return A propensity score model object
#' @export
#'
#' @examples
#' \dontrun{
#' psModel <- buildPropensityModel(
#'   covariateData = covariateData,
#'   treatmentVariable = treatment,
#'   prior = "laplace",
#'   priorVariance = 0.1
#' )
#' }
buildPropensityModel <- function(covariateData,
                                 treatmentVariable,
                                 prior = "laplace",
                                 priorVariance = 0.1,
                                 control = NULL) {

  checkmate::assertClass(covariateData, "CovariateData")
  checkmate::assertNumeric(treatmentVariable, lower = 0, upper = 1)
  checkmate::assertChoice(prior, c("laplace", "normal"))

  message("Building large-scale propensity score model...")

  if (is.null(control)) {
    control <- Cyclops::createControl(
      noiseLevel = "quiet",
      cvType = "auto",
      startingVariance = priorVariance,
      tolerance = 2e-07,
      cvRepetitions = 10,
      threads = max(1, parallel::detectCores() - 1)
    )
  }

  # Convert covariate data to Cyclops format
  cyclopsData <- Cyclops::convertToCyclopsData(
    outcomes = treatmentVariable,
    covariates = covariateData$covariates,
    modelType = "lr",
    addIntercept = TRUE
  )

  # Fit model with regularization
  psModel <- Cyclops::fitCyclopsModel(
    cyclopsData = cyclopsData,
    prior = Cyclops::createPrior(
      priorType = prior,
      variance = priorVariance,
      exclude = 0  # Don't regularize intercept
    ),
    control = control
  )

  message(sprintf("Model fitted with %d non-zero coefficients",
                                  sum(coef(psModel) != 0)))

  # Calculate propensity scores
  ps <- stats::predict(psModel, cyclopsData)

  # Add diagnostics
  psModel$propensityScores <- ps
  psModel$treatmentVariable <- treatmentVariable
  psModel$covariateData <- covariateData

  # Calculate AUC
  auc <- calculateAUC(ps, treatmentVariable)
  message(sprintf("Propensity score model AUC: %.3f", auc))

  psModel$auc <- auc

  class(psModel) <- c("SyTrialPSModel", class(psModel))

  return(psModel)
}

#' Calculate Overlap Weights
#'
#' @description
#' Calculates overlap weights (Li et al. 2018) which optimize balance
#' and precision by focusing on the population with equipoise.
#'
#' @param propensityScores Vector of propensity scores
#' @param treatment Binary treatment indicator
#' @param epsilon Small value used to truncate propensity scores away from 0 and 1
#'
#' @return Vector of overlap weights
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' treatment <- rbinom(n, 1, 0.4)
#' propensityScores <- pmin(pmax(rbeta(n, 2, 5), 1e-6), 1 - 1e-6)
#' w <- calculateOverlapWeights(propensityScores, treatment)
#' summary(w)
calculateOverlapWeights <- function(propensityScores, treatment, epsilon = 1e-6) {

  checkmate::assertNumeric(propensityScores, lower = 0, upper = 1)
  checkmate::assertNumeric(treatment, lower = 0, upper = 1)
  checkmate::assertNumber(epsilon, lower = 0, finite = TRUE)

  message("Calculating overlap weights...")

  # Truncate propensity scores away from 0 and 1 to avoid zero weights
  ps <- pmin(pmax(propensityScores, epsilon), 1 - epsilon)
  if (any(ps != propensityScores)) {
    message(sprintf("Truncated propensity scores to [%.2e, %.2e]", epsilon, 1 - epsilon))
  }

  # Overlap weights: w = (1-ps) for treated, ps for controls
  weights <- ifelse(treatment == 1,
                    1 - ps,
                    ps)

  # Normalize weights within treatment groups
  weightsNormalized <- weights
  weightsNormalized[treatment == 1] <- weights[treatment == 1] / mean(weights[treatment == 1])
  weightsNormalized[treatment == 0] <- weights[treatment == 0] / mean(weights[treatment == 0])

  message(sprintf("Effective sample size: %.1f (treated: %.1f, control: %.1f)",
                                  sum(weightsNormalized),
                                  sum(weightsNormalized[treatment == 1]),
                                  sum(weightsNormalized[treatment == 0])))

  return(weightsNormalized)
}

#' Calculate Inverse Probability of Treatment Weights (IPTW)
#'
#' @param propensityScores Vector of propensity scores
#' @param treatment Binary treatment indicator
#' @param stabilized Logical; use stabilized weights (default: TRUE)
#' @param epsilon Small value used to truncate propensity scores away from 0 and 1
#'
#' @return Vector of IPTW weights
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' treatment <- rbinom(n, 1, 0.4)
#' propensityScores <- pmin(pmax(rbeta(n, 2, 5), 1e-6), 1 - 1e-6)
#' w <- calculateIPTW(propensityScores, treatment, stabilized = TRUE)
#' summary(w)
calculateIPTW <- function(propensityScores, treatment, stabilized = TRUE, epsilon = 1e-6) {

  checkmate::assertNumeric(propensityScores, lower = 0, upper = 1)
  checkmate::assertNumeric(treatment, lower = 0, upper = 1)
  checkmate::assertFlag(stabilized)
  checkmate::assertNumber(epsilon, lower = 0, finite = TRUE)

  message("Calculating IPTW weights...")

  # Truncate propensity scores away from 0 and 1 to avoid infinite weights
  ps <- pmin(pmax(propensityScores, epsilon), 1 - epsilon)
  if (any(ps != propensityScores)) {
    message(sprintf("Truncated propensity scores to [%.2e, %.2e]", epsilon, 1 - epsilon))
  }

  # Standard IPTW
  weights <- ifelse(treatment == 1,
                    1 / ps,
                    1 / (1 - ps))

  # Stabilized weights
  if (stabilized) {
    pTreatment <- mean(treatment)
    weights <- ifelse(treatment == 1,
                      pTreatment / ps,
                      (1 - pTreatment) / (1 - ps))
  }

  # Trim extreme weights (optional - at 99th percentile)
  weightCap <- quantile(weights, 0.99)
  weightsTrimmed <- pmin(weights, weightCap)

  nTrimmed <- sum(weights > weightCap)
  if (nTrimmed > 0) {
    message(sprintf("Trimmed %d extreme weights at %.2f", nTrimmed, weightCap))
  }

  message(sprintf("Effective sample size: %.1f", sum(weightsTrimmed)^2 / sum(weightsTrimmed^2)))

  return(weightsTrimmed)
}

# Internal function: Calculate AUC
calculateAUC <- function(predictions, labels) {
  order_idx <- order(predictions, decreasing = TRUE)
  predictions <- predictions[order_idx]
  labels <- labels[order_idx]

  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)

  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }

  tp <- cumsum(labels == 1)
  fp <- cumsum(labels == 0)

  tpr <- tp / n_pos
  fpr <- fp / n_neg

  auc <- sum(diff(c(0, fpr)) * (c(0, tpr[-length(tpr)]) + tpr) / 2)

  return(auc)
}
