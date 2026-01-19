#' Perform Targeted Maximum Likelihood Estimation (TMLE)
#'
#' @description
#' Implements TMLE for causal inference with double robustness property.
#' Integrates with Super Learner for flexible machine learning.
#'
#' @param data Data frame containing treatment, outcome, and covariates
#' @param treatment Character; name of treatment variable
#' @param outcome Character; name of outcome variable
#' @param covariates Character vector; names of covariate columns
#' @param family Outcome family ("binomial", "gaussian", "survival")
#' @param gbound Bounds for propensity score truncation (default: c(0.01, 0.99))
#' @param Q.SL.library Super Learner library for outcome model
#' @param g.SL.library Super Learner library for treatment model
#' @param cvQinit Logical; cross-validate outcome model (default: TRUE)
#'
#' @return A TMLE result object with ATE, confidence intervals, and diagnostics
#' @export
#'
#' @examples
#' \dontrun{
#' tmleResult <- performTMLE(
#'   data = analysisData,
#'   treatment = "treatment",
#'   outcome = "outcome",
#'   covariates = covariateNames,
#'   family = "binomial",
#'   Q.SL.library = c("SL.glm", "SL.glmnet", "SL.ranger"),
#'   g.SL.library = c("SL.glm", "SL.glmnet")
#' )
#' }
performTMLE <- function(data,
                        treatment,
                        outcome,
                        covariates,
                        family = "binomial",
                        gbound = c(0.01, 0.99),
                        Q.SL.library = c("SL.glm", "SL.glmnet", "SL.ranger"),
                        g.SL.library = c("SL.glm", "SL.glmnet"),
                        cvQinit = TRUE) {

  checkmate::assertDataFrame(data)
  checkmate::assertChoice(treatment, colnames(data))
  checkmate::assertChoice(outcome, colnames(data))
  checkmate::assertSubset(covariates, colnames(data))
  checkmate::assertChoice(family, c("binomial", "gaussian", "survival"))

  message("Performing Targeted Maximum Likelihood Estimation (TMLE)...")

  # Check if required packages are available
  if (!requireNamespace("tmle", quietly = TRUE)) {
    stop("Package 'tmle' is required for TMLE estimation. Please install it.")
  }

  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package 'SuperLearner' is required for Super Learner. Please install it.")
  }

  # Prepare data
  A <- data[[treatment]]
  Y <- data[[outcome]]
  W <- as.matrix(data[, covariates, drop = FALSE])

  # Handle missing data
  completeIdx <- complete.cases(cbind(A, Y, W))
  if (sum(!completeIdx) > 0) {
    ParallelLogger::logWarn(sprintf("Removing %d rows with missing data", sum(!completeIdx)))
    A <- A[completeIdx]
    Y <- Y[completeIdx]
    W <- W[completeIdx, , drop = FALSE]
  }

  # Fit TMLE
  message("Fitting TMLE with Super Learner...")

  tmleResult <- tryCatch({
    tmle::tmle(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      family = family,
      gbound = gbound,
      cvQinit = cvQinit,
      verbose = FALSE
    )
  }, error = function(e) {
    ParallelLogger::logError(sprintf("TMLE estimation failed: %s", e$message))
    return(NULL)
  })

  if (is.null(tmleResult)) {
    stop("TMLE estimation failed. Check data and model specifications.")
  }

  # Extract results
  ate <- tmleResult$estimates$ATE$psi
  ateSE <- sqrt(tmleResult$estimates$ATE$var.psi)
  ci <- tmleResult$estimates$ATE$CI
  pValue <- tmleResult$estimates$ATE$pvalue

  message(sprintf("TMLE ATE: %.4f (95%% CI: %.4f to %.4f, p = %.4f)",
                                  ate, ci[1], ci[2], pValue))

  # Package results
  result <- list(
    ate = ate,
    se = ateSE,
    ci = ci,
    pValue = pValue,
    tmleObject = tmleResult,
    treatment = treatment,
    outcome = outcome,
    covariates = covariates,
    family = family,
    n = length(A),
    nTreated = sum(A == 1),
    nControl = sum(A == 0)
  )

  class(result) <- "SyTrialTMLE"

  return(result)
}

#' Perform G-Computation for Causal Inference
#'
#' @description
#' Implements parametric G-computation (standardization) for causal effect estimation.
#'
#' @param data Data frame containing treatment, outcome, and covariates
#' @param treatment Character; name of treatment variable
#' @param outcome Character; name of outcome variable
#' @param covariates Character vector; names of covariate columns
#' @param family Outcome family ("binomial", "gaussian", "survival")
#' @param outcomeModel Formula or model specification (optional)
#' @param nBootstrap Number of bootstrap iterations for CI (default: 1000)
#'
#' @return A G-computation result object
#' @export
performGComputation <- function(data,
                                treatment,
                                outcome,
                                covariates,
                                family = "binomial",
                                outcomeModel = NULL,
                                nBootstrap = 1000) {

  checkmate::assertDataFrame(data)
  checkmate::assertChoice(treatment, colnames(data))
  checkmate::assertChoice(outcome, colnames(data))
  checkmate::assertSubset(covariates, colnames(data))

  message("Performing G-Computation...")

  # Build outcome model
  if (is.null(outcomeModel)) {
    formulaStr <- sprintf("%s ~ %s + %s", outcome, treatment, paste(covariates, collapse = " + "))
    outcomeModel <- as.formula(formulaStr)
  }

  # Fit outcome model
  if (family == "binomial") {
    fit <- glm(outcomeModel, data = data, family = binomial())
  } else if (family == "gaussian") {
    fit <- glm(outcomeModel, data = data, family = gaussian())
  } else {
    stop("Unsupported family for G-computation")
  }

  # Predict under treatment
  dataTreated <- data
  dataTreated[[treatment]] <- 1
  predTreated <- predict(fit, newdata = dataTreated, type = "response")

  # Predict under control
  dataControl <- data
  dataControl[[treatment]] <- 0
  predControl <- predict(fit, newdata = dataControl, type = "response")

  # Calculate ATE
  ate <- mean(predTreated - predControl)

  message(sprintf("G-Computation ATE: %.4f", ate))

  # Bootstrap confidence intervals
  message(sprintf("Computing bootstrap CI with %d iterations...", nBootstrap))

  bootstrapATEs <- numeric(nBootstrap)

  for (i in 1:nBootstrap) {
    bootIdx <- sample(nrow(data), replace = TRUE)
    bootData <- data[bootIdx, ]

    bootFit <- tryCatch({
      if (family == "binomial") {
        glm(outcomeModel, data = bootData, family = binomial())
      } else {
        glm(outcomeModel, data = bootData, family = gaussian())
      }
    }, error = function(e) NULL)

    if (!is.null(bootFit)) {
      bootDataTreated <- bootData
      bootDataTreated[[treatment]] <- 1
      bootPredTreated <- predict(bootFit, newdata = bootDataTreated, type = "response")

      bootDataControl <- bootData
      bootDataControl[[treatment]] <- 0
      bootPredControl <- predict(bootFit, newdata = bootDataControl, type = "response")

      bootstrapATEs[i] <- mean(bootPredTreated - bootPredControl)
    } else {
      bootstrapATEs[i] <- NA
    }
  }

  bootstrapATEs <- bootstrapATEs[!is.na(bootstrapATEs)]
  ci <- quantile(bootstrapATEs, c(0.025, 0.975))
  se <- sd(bootstrapATEs)

  message(sprintf("G-Computation ATE: %.4f (95%% CI: %.4f to %.4f)",
                                  ate, ci[1], ci[2]))

  result <- list(
    ate = ate,
    se = se,
    ci = ci,
    outcomeModel = fit,
    treatment = treatment,
    outcome = outcome,
    covariates = covariates,
    family = family,
    n = nrow(data),
    nTreated = sum(data[[treatment]] == 1),
    nControl = sum(data[[treatment]] == 0),
    bootstrapATEs = bootstrapATEs
  )

  class(result) <- "SyTrialGComp"

  return(result)
}
