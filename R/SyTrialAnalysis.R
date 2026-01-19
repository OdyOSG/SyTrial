#' Run Complete SyTrial Analysis Pipeline
#'
#' @description
#' Executes the complete synthetic control arm analysis workflow, from
#' cohort definition through causal inference and diagnostics.
#'
#' @param syTrialConnection SyTrialConnection object
#' @param treatmentCohortId Treatment cohort definition ID
#' @param controlCohortId Control cohort definition ID (may lack index dates)
#' @param imputationMethod Index date imputation method
#' @param covariateSettings Covariate extraction settings
#' @param causalMethod Causal inference method ("tmle", "gcomp", "iptw")
#' @param weightingMethod Weighting method ("overlap", "iptw", "matching")
#' @param outcomeType Outcome type ("binary", "continuous", "survival")
#' @param outcomeVariable Name of outcome variable
#' @param timeVariable Name of time variable (for survival outcomes)
#' @param eventVariable Name of event indicator (for survival outcomes)
#' @param outputDir Directory for diagnostic outputs
#'
#' @return A SyTrialResult object
#' @export
#'
#' @examples
#' \dontrun{
#' # Connect to database
#' syTrialConn <- createSyTrialConnection(
#'   connectionDetails = connectionDetails,
#'   cdmDatabaseSchema = "cdm_synpuf",
#'   cohortDatabaseSchema = "results",
#'   cohortTable = "cohort"
#' )
#'
#' # Run analysis
#' result <- runSyTrialAnalysis(
#'   syTrialConnection = syTrialConn,
#'   treatmentCohortId = 1,
#'   controlCohortId = 2,
#'   imputationMethod = "risk_set_sampling",
#'   causalMethod = "tmle",
#'   weightingMethod = "overlap",
#'   outcomeType = "survival",
#'   outcomeVariable = "death",
#'   timeVariable = "time_to_event",
#'   eventVariable = "event_indicator",
#'   outputDir = "./sytrial_results"
#' )
#'
#' # View results
#' print(result)
#' summary(result)
#'
#' # Disconnect
#' disconnectSyTrial(syTrialConn)
#' }
runSyTrialAnalysis <- function(syTrialConnection,
                               treatmentCohortId,
                               controlCohortId,
                               imputationMethod = "risk_set_sampling",
                               covariateSettings = NULL,
                               causalMethod = "tmle",
                               weightingMethod = "overlap",
                               outcomeType = "binary",
                               outcomeVariable = NULL,
                               timeVariable = NULL,
                               eventVariable = NULL,
                               outcomeData = NULL,
                               outputDir = "./sytrial_results") {

  # Input validation
  checkmate::assertClass(syTrialConnection, "SyTrialConnection")
  checkmate::assertInt(treatmentCohortId)
  checkmate::assertInt(controlCohortId)
  checkmate::assertChoice(causalMethod, c("tmle", "gcomp", "iptw"))
  checkmate::assertChoice(weightingMethod, c("overlap", "iptw", "matching"))
  checkmate::assertChoice(outcomeType, c("binary", "continuous", "survival"))

  # Create output directory
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }
  message(paste0("=", strrep("=", 60)))
  message("SyTrial Analysis Pipeline Started")
  message(paste0("=", strrep("=", 60)))

  startTime <- Sys.time()

  # Step 1: Impute control index dates
  message("\n[Step 1/6] Imputing control cohort index dates...")
  controlWithIndex <- imputeControlIndexDates(
    syTrialConnection = syTrialConnection,
    treatmentCohortId = treatmentCohortId,
    controlCohortId = controlCohortId,
    imputationMethod = imputationMethod
  )

  # Step 2: Extract features
  message("\n[Step 2/6] Extracting high-dimensional covariates...")
  if (is.null(covariateSettings)) {
    covariateSettings <- createDefaultCovariateSettings()
  }

  covariateData <- extractSCAFeatures(
    syTrialConnection = syTrialConnection,
    treatmentCohortId = treatmentCohortId,
    controlWithIndex = controlWithIndex,
    covariateSettings = covariateSettings
  )

  # Step 3: Build propensity score model
  message("\n[Step 3/6] Building propensity score model...")

  # Extract treatment indicator (aligned to covariateData rows)
  subjectIds <- covariateData$covariates$rowId
  treatment <- extractTreatmentIndicator(
    syTrialConnection = syTrialConnection,
    treatmentCohortId = treatmentCohortId,
    subjectIds = subjectIds
  )

  psModel <- buildPropensityModel(
    covariateData = covariateData,
    treatmentVariable = treatment,
    prior = "laplace",
    priorVariance = 0.1
  )

  propensityScores <- psModel$propensityScores

  # Step 4: Calculate weights
  message("\n[Step 4/6] Calculating weights...")
  weights <- switch(
    weightingMethod,
    "overlap" = calculateOverlapWeights(propensityScores, treatment),
    "iptw" = calculateIPTW(propensityScores, treatment, stabilized = TRUE),
    "matching" = performMatching(propensityScores, treatment)
  )

  # Step 5: Causal inference
  message("\n[Step 5/6] Performing causal inference...")

  # Prepare analysis data
  analysisData <- prepareAnalysisData(
    syTrialConnection = syTrialConnection,
    treatmentCohortId = treatmentCohortId,
    controlWithIndex = controlWithIndex,
    covariateData = covariateData,
    outcomeVariable = outcomeVariable,
    timeVariable = timeVariable,
    eventVariable = eventVariable,
    outcomeData = outcomeData
  )

  causalResult <- switch(
    causalMethod,
    "tmle" = performTMLE(
      data = analysisData,
      treatment = "treatment",
      outcome = outcomeVariable,
      covariates = names(analysisData)[!names(analysisData) %in% c("treatment", outcomeVariable, timeVariable, eventVariable)],
      family = ifelse(outcomeType == "binary", "binomial", "gaussian")
    ),
    "gcomp" = performGComputation(
      data = analysisData,
      treatment = "treatment",
      outcome = outcomeVariable,
      covariates = names(analysisData)[!names(analysisData) %in% c("treatment", outcomeVariable, timeVariable, eventVariable)],
      family = ifelse(outcomeType == "binary", "binomial", "gaussian")
    ),
    "iptw" = performIPTWAnalysis(
      data = analysisData,
      treatment = "treatment",
      outcome = outcomeVariable,
      weights = weights,
      family = if (outcomeType == "continuous") stats::gaussian() else stats::binomial()
    )
  )

  # Step 6: Generate diagnostics
  message("\n[Step 6/6] Generating diagnostic reports...")

  # Package results
  result <- list(
    treatmentCohortId = treatmentCohortId,
    controlCohortId = controlCohortId,
    imputationMethod = imputationMethod,
    causalMethod = causalMethod,
    weightingMethod = weightingMethod,
    outcomeType = outcomeType,
    n = nrow(analysisData),
    nTreated = sum(treatment == 1),
    nControl = sum(treatment == 0),
    covariateData = covariateData,
    treatment = treatment,
    propensityScores = propensityScores,
    psAUC = psModel$auc,
    weights = weights,
    covariates = names(analysisData)[!names(analysisData) %in% c("treatment", outcomeVariable, timeVariable, eventVariable)],
    ate = causalResult$ate,
    se = causalResult$se,
    ci = causalResult$ci,
    pValue = causalResult$pValue,
    causalResult = causalResult,
    analysisData = analysisData,
    survivalData = if (outcomeType == "survival") {
      data.frame(
        time = analysisData[[timeVariable]],
        event = analysisData[[eventVariable]],
        treatment = treatment
      )
    } else {
      NULL
    },
    outputDir = outputDir,
    executionTime = difftime(Sys.time(), startTime, units = "mins")
  )

  class(result) <- "SyTrialResult"

  # Generate diagnostics
  diagnostics <- generateDiagnosticReport(result, outputDir)
  result$diagnostics <- diagnostics

  # Save result object
  saveRDS(result, file.path(outputDir, "sytrial_result.rds"))

  endTime <- Sys.time()
  executionTime <- difftime(endTime, startTime, units = "mins")

  message(paste0("\n", strrep("=", 60)))
  message(sprintf("SyTrial Analysis Completed in %.2f minutes", executionTime))
  message(paste0("=", strrep("=", 60)))

  try(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER"), silent = TRUE)

  return(result)
}

# Internal helper functions

extractTreatmentIndicator <- function(syTrialConnection, treatmentCohortId, subjectIds) {
  checkmate::assertClass(syTrialConnection, "SyTrialConnection")
  checkmate::assertInt(treatmentCohortId)
  checkmate::assertVector(subjectIds, any.missing = FALSE)

  # Extract treatment subjects
  sqlTreatment <- "
    SELECT subject_id
    FROM @cohort_database_schema.@cohort_table
    WHERE cohort_definition_id = @treatment_cohort_id
  "

  sqlTreatment <- SqlRender::render(
    sqlTreatment,
    cohort_database_schema = syTrialConnection$cohortDatabaseSchema,
    cohort_table = syTrialConnection$cohortTable,
    treatment_cohort_id = treatmentCohortId
  )

  sqlTreatment <- SqlRender::translate(sqlTreatment, targetDialect = syTrialConnection$dbms)
  treatmentSubjects <- DatabaseConnector::querySql(syTrialConnection$connection, sqlTreatment)

  # Ensure consistent column naming across DBMS / DatabaseConnector settings
  treatmentSubjects <- SqlRender::snakeCaseToCamelCase(treatmentSubjects)
  if (!("subjectId" %in% names(treatmentSubjects))) {
    stop("Expected subjectId column in treatment cohort query result.")
  }

  treatedSubjectIds <- treatmentSubjects$subjectId
  treatedSubjectIds <- treatedSubjectIds[!is.na(treatedSubjectIds)]

  # Return indicator aligned to subjectIds order
  return(as.integer(subjectIds %in% treatedSubjectIds))
}

prepareAnalysisData <- function(syTrialConnection, treatmentCohortId, controlWithIndex,
                                covariateData, outcomeVariable, timeVariable, eventVariable,
                                outcomeData = NULL) {

  # Extract covariates as data frame
  # Note: this converts to a (potentially large) dense matrix; for large studies consider
  # using sparse workflows end-to-end.
  covSparse <- FeatureExtraction::covariateDataToSparseMatrix(covariateData)
  covDF <- as.data.frame(as.matrix(covSparse))

  # Add subject IDs aligned to covariate rows (rowId are the subject IDs)
  covDF$subjectId <- covariateData$covariates$rowId

  # Add treatment indicator aligned to covariate rows
  covDF$treatment <- extractTreatmentIndicator(
    syTrialConnection = syTrialConnection,
    treatmentCohortId = treatmentCohortId,
    subjectIds = covDF$subjectId
  )

  # Outcome data must be provided explicitly
  if (is.null(outcomeData)) {
    stop("No outcomeData provided. Please supply outcomeData with at least subjectId and the requested outcome/time/event columns.")
  }
  checkmate::assertDataFrame(outcomeData)
  checkmate::assertNames(names(outcomeData), must.include = "subjectId")

  requiredCols <- c(outcomeVariable, timeVariable, eventVariable)
  requiredCols <- requiredCols[!is.null(requiredCols)]
  missingCols <- setdiff(requiredCols, names(outcomeData))
  if (length(missingCols) > 0) {
    stop(sprintf(
      "outcomeData is missing required columns: %s",
      paste(missingCols, collapse = ", ")
    ))
  }

  covDF <- merge(covDF, outcomeData, by = "subjectId", all.x = TRUE, sort = FALSE)

  # Ensure outcome columns are present after merge
  if (!is.null(outcomeVariable) && any(is.na(covDF[[outcomeVariable]]))) {
    stop("Missing outcome values after merging outcomeData. Ensure outcomeData covers all subjects in the covariate data.")
  }
  if (!is.null(timeVariable) && any(is.na(covDF[[timeVariable]]))) {
    stop("Missing time values after merging outcomeData. Ensure outcomeData covers all subjects in the covariate data.")
  }
  if (!is.null(eventVariable) && any(is.na(covDF[[eventVariable]]))) {
    stop("Missing event values after merging outcomeData. Ensure outcomeData covers all subjects in the covariate data.")
  }

  return(covDF)
}

performMatching <- function(propensityScores, treatment) {
  # Simplified 1:1 matching
  # In practice, use MatchIt or similar

  treatedIdx <- which(treatment == 1)
  controlIdx <- which(treatment == 0)

  weights <- rep(0, length(treatment))

  for (i in treatedIdx) {
    # Find nearest control
    distances <- abs(propensityScores[controlIdx] - propensityScores[i])
    matchIdx <- controlIdx[which.min(distances)]

    weights[i] <- 1
    weights[matchIdx] <- 1
  }

  return(weights)
}

performIPTWAnalysis <- function(data, treatment, outcome, weights, family = stats::binomial()) {
  # Weighted outcome regression
  formula <- as.formula(paste(outcome, "~", treatment))

  fit <- glm(formula, data = data, weights = weights, family = family)

  # Extract treatment effect
  if (!(treatment %in% names(coef(fit)))) {
    stop(sprintf("Treatment coefficient '%s' not found in model coefficients.", treatment))
  }
  ate <- unname(coef(fit)[treatment])
  se <- sqrt(vcov(fit)[treatment, treatment])
  ci <- ate + c(-1.96, 1.96) * se
  pValue <- 2 * pnorm(-abs(ate / se))

  result <- list(
    ate = ate,
    se = se,
    ci = ci,
    pValue = pValue,
    model = fit
  )

  class(result) <- "SyTrialIPTW"

  return(result)
}

#' Print SyTrial Result
#'
#' @param x SyTrialResult object
#' @param ... Additional arguments
#' @export
print.SyTrialResult <- function(x, ...) {
  cat("\n")
  cat("=" %s% strrep("=", 60), "\n")
  cat("SyTrial Analysis Results\n")
  cat("=" %s% strrep("=", 60), "\n\n")

  cat(sprintf("Treatment Cohort ID: %d\n", x$treatmentCohortId))
  cat(sprintf("Control Cohort ID: %d\n", x$controlCohortId))
  cat(sprintf("Imputation Method: %s\n", x$imputationMethod))
  cat(sprintf("Causal Method: %s\n", x$causalMethod))
  cat(sprintf("Weighting Method: %s\n\n", x$weightingMethod))

  cat("Sample Size:\n")
  cat(sprintf("  Total: %d\n", x$n))
  cat(sprintf("  Treatment: %d\n", x$nTreated))
  cat(sprintf("  Control: %d\n\n", x$nControl))

  cat("Propensity Score Model:\n")
  cat(sprintf("  AUC: %.3f\n", x$psAUC))
  cat(sprintf("  Effective Sample Size: %.1f\n\n", sum(x$weights)^2 / sum(x$weights^2)))

  cat("Causal Effect Estimate:\n")
  cat(sprintf("  ATE: %.4f\n", x$ate))
  cat(sprintf("  SE: %.4f\n", x$se))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", x$ci[1], x$ci[2]))
  cat(sprintf("  P-value: %.4f\n\n", x$pValue))

  cat(sprintf("Execution Time: %.2f minutes\n", x$executionTime))
  cat(sprintf("Output Directory: %s\n", x$outputDir))

  cat("\n")
  cat("=" %s% strrep("=", 60), "\n")
}

#' Summary of SyTrial Result
#'
#' @param object SyTrialResult object
#' @param ... Additional arguments
#' @export
summary.SyTrialResult <- function(object, ...) {
  print(object)

  cat("\nCovariate Balance Summary:\n")
  # Add balance summary here

  cat("\nDiagnostic Files:\n")
  cat(sprintf("  - %s\n", file.path(object$outputDir, "love_plot.png")))
  cat(sprintf("  - %s\n", file.path(object$outputDir, "ps_distribution.png")))
  if (!is.null(object$survivalData)) {
    cat(sprintf("  - %s\n", file.path(object$outputDir, "kaplan_meier.png")))
  }
  cat(sprintf("  - %s\n", file.path(object$outputDir, "summary_statistics.csv")))
  cat(sprintf("  - %s\n", file.path(object$outputDir, "sytrial_result.rds")))
}

