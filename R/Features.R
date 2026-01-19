#' Extract High-Dimensional Covariates for SCA
#'
#' @description
#' Extracts comprehensive covariate set using HADES FeatureExtraction package.
#' Includes demographics, conditions, drugs, procedures, measurements, and
#' derived risk scores.
#'
#' @param syTrialConnection SyTrialConnection object
#' @param treatmentCohortId Treatment cohort ID
#' @param controlWithIndex Control cohort data frame with imputed index dates
#' @param covariateSettings FeatureExtraction covariate settings (optional)
#' @param aggregated Logical; use aggregated covariates (default: FALSE)
#'
#' @return A CovariateData object
#' @export
#'
#' @examples
#' \dontrun{
#' covariateData <- extractSCAFeatures(
#'   syTrialConnection = syTrialConn,
#'   treatmentCohortId = 1,
#'   controlWithIndex = controlData,
#'   covariateSettings = createDefaultCovariateSettings()
#' )
#' }
extractSCAFeatures <- function(syTrialConnection,
                               treatmentCohortId,
                               controlWithIndex,
                               covariateSettings = NULL,
                               aggregated = FALSE) {

  checkmate::assertClass(syTrialConnection, "SyTrialConnection")

  if (is.null(covariateSettings)) {
    covariateSettings <- createDefaultCovariateSettings()
  }

  message("Creating combined cohort for feature extraction...")

  # Create temporary cohort table with both treatment and control
  combinedCohort <- createCombinedCohort(
    syTrialConnection,
    treatmentCohortId,
    controlWithIndex
  )

  message("Extracting high-dimensional covariates...")

  covariateData <- FeatureExtraction::getDbCovariateData(
    connection = syTrialConnection$connection,
    cdmDatabaseSchema = syTrialConnection$cdmDatabaseSchema,
    cohortDatabaseSchema = syTrialConnection$cohortDatabaseSchema,
    cohortTable = "sytrial_temp_cohort",
    cohortId = -1,  # Use all subjects in temp table
    covariateSettings = covariateSettings,
    aggregated = aggregated
  )

  message(sprintf("Extracted %d covariates for %d subjects",
                                  nrow(covariateData$covariateRef),
                                  length(unique(covariateData$covariates$rowId))))

  # Clean up temporary table
  cleanupTempCohort(syTrialConnection)

  return(covariateData)
}

#' Create Default Covariate Settings for SCA
#'
#' @description
#' Creates comprehensive covariate settings optimized for synthetic control
#' arm construction, including temporal features and risk scores.
#'
#' @param useConditionOccurrence Include condition occurrences (default: TRUE)
#' @param useDrugExposure Include drug exposures (default: TRUE)
#' @param useProcedureOccurrence Include procedures (default: TRUE)
#' @param useMeasurement Include measurements (default: TRUE)
#' @param useObservation Include observations (default: TRUE)
#' @param useCharlsonIndex Include Charlson Comorbidity Index (default: TRUE)
#' @param useDcsi Include Diabetes Comorbidity Severity Index (default: TRUE)
#' @param useChads2Vasc Include CHADS2VASc score (default: TRUE)
#' @param temporalStartDays Start of temporal window (default: -365)
#' @param temporalEndDays End of temporal window (default: 0)
#' @param includedCovariateConceptIds Specific concept IDs to include (optional)
#' @param excludedCovariateConceptIds Specific concept IDs to exclude (optional)
#'
#' @return A covariateSettings object
#' @export
createDefaultCovariateSettings <- function(useConditionOccurrence = TRUE,
                                           useDrugExposure = TRUE,
                                           useProcedureOccurrence = TRUE,
                                           useMeasurement = TRUE,
                                           useObservation = TRUE,
                                           useCharlsonIndex = TRUE,
                                           useDcsi = TRUE,
                                           useChads2Vasc = TRUE,
                                           temporalStartDays = -365,
                                           temporalEndDays = 0,
                                           includedCovariateConceptIds = c(),
                                           excludedCovariateConceptIds = c()) {

  covariateSettings <- FeatureExtraction::createDefaultCovariateSettings(
    useConditionOccurrenceAnyTimePrior = useConditionOccurrence,
    useConditionOccurrenceLongTerm = useConditionOccurrence,
    useConditionOccurrenceMediumTerm = useConditionOccurrence,
    useConditionOccurrenceShortTerm = useConditionOccurrence,
    useConditionEraAnyTimePrior = useConditionOccurrence,
    useConditionEraLongTerm = useConditionOccurrence,
    useConditionEraMediumTerm = useConditionOccurrence,
    useConditionEraShortTerm = useConditionOccurrence,
    useConditionEraOverlapping = useConditionOccurrence,
    useConditionGroupEraAnyTimePrior = useConditionOccurrence,
    useConditionGroupEraLongTerm = useConditionOccurrence,
    useConditionGroupEraMediumTerm = useConditionOccurrence,
    useConditionGroupEraShortTerm = useConditionOccurrence,
    useConditionGroupEraOverlapping = useConditionOccurrence,

    useDrugExposureAnyTimePrior = useDrugExposure,
    useDrugExposureLongTerm = useDrugExposure,
    useDrugExposureMediumTerm = useDrugExposure,
    useDrugExposureShortTerm = useDrugExposure,
    useDrugEraAnyTimePrior = useDrugExposure,
    useDrugEraLongTerm = useDrugExposure,
    useDrugEraMediumTerm = useDrugExposure,
    useDrugEraShortTerm = useDrugExposure,
    useDrugEraOverlapping = useDrugExposure,
    useDrugGroupEraAnyTimePrior = useDrugExposure,
    useDrugGroupEraLongTerm = useDrugExposure,
    useDrugGroupEraMediumTerm = useDrugExposure,
    useDrugGroupEraShortTerm = useDrugExposure,
    useDrugGroupEraOverlapping = useDrugExposure,

    useProcedureOccurrenceAnyTimePrior = useProcedureOccurrence,
    useProcedureOccurrenceLongTerm = useProcedureOccurrence,
    useProcedureOccurrenceMediumTerm = useProcedureOccurrence,
    useProcedureOccurrenceShortTerm = useProcedureOccurrence,

    useMeasurementAnyTimePrior = useMeasurement,
    useMeasurementLongTerm = useMeasurement,
    useMeasurementMediumTerm = useMeasurement,
    useMeasurementShortTerm = useMeasurement,
    useMeasurementValueAnyTimePrior = useMeasurement,
    useMeasurementValueLongTerm = useMeasurement,
    useMeasurementValueMediumTerm = useMeasurement,
    useMeasurementValueShortTerm = useMeasurement,

    useObservationAnyTimePrior = useObservation,
    useObservationLongTerm = useObservation,
    useObservationMediumTerm = useObservation,
    useObservationShortTerm = useObservation,

    useCharlsonIndex = useCharlsonIndex,
    useDcsi = useDcsi,
    useChads2Vasc = useChads2Vasc,

    useDistinctConditionCountLongTerm = TRUE,
    useDistinctConditionCountMediumTerm = TRUE,
    useDistinctConditionCountShortTerm = TRUE,
    useDistinctIngredientCountLongTerm = TRUE,
    useDistinctIngredientCountMediumTerm = TRUE,
    useDistinctIngredientCountShortTerm = TRUE,
    useDistinctProcedureCountLongTerm = TRUE,
    useDistinctProcedureCountMediumTerm = TRUE,
    useDistinctProcedureCountShortTerm = TRUE,
    useDistinctMeasurementCountLongTerm = TRUE,
    useDistinctMeasurementCountMediumTerm = TRUE,
    useDistinctMeasurementCountShortTerm = TRUE,

    useVisitCountLongTerm = TRUE,
    useVisitCountMediumTerm = TRUE,
    useVisitCountShortTerm = TRUE,
    useVisitConceptCountLongTerm = TRUE,
    useVisitConceptCountMediumTerm = TRUE,
    useVisitConceptCountShortTerm = TRUE,

    longTermStartDays = -365,
    mediumTermStartDays = -180,
    shortTermStartDays = -30,
    endDays = 0,

    includedCovariateConceptIds = includedCovariateConceptIds,
    excludedCovariateConceptIds = excludedCovariateConceptIds,
    addDescendantsToInclude = TRUE,
    addDescendantsToExclude = TRUE
  )

  return(covariateSettings)
}

# Internal function: Create combined cohort
createCombinedCohort <- function(syTrialConnection, treatmentCohortId, controlWithIndex) {

  # Extract treatment cohort
  sqlTreatment <- "
    SELECT
      subject_id,
      cohort_start_date,
      cohort_end_date,
      1 AS treatment
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
  treatmentData <- DatabaseConnector::querySql(syTrialConnection$connection, sqlTreatment)
  colnames(treatmentData) <- SqlRender::snakeCaseToCamelCase(colnames(treatmentData))

  # Prepare control data
  controlData <- data.frame(
    subjectId = controlWithIndex$subjectId,
    cohortStartDate = controlWithIndex$imputedIndexDate,
    cohortEndDate = controlWithIndex$latestEligibleDate,
    treatment = 0
  )

  # Combine
  combinedData <- rbind(treatmentData, controlData)

  # Create temporary table
  DatabaseConnector::insertTable(
    connection = syTrialConnection$connection,
    databaseSchema = syTrialConnection$cohortDatabaseSchema,
    tableName = "sytrial_temp_cohort",
    data = combinedData,
    dropTableIfExists = TRUE,
    createTable = TRUE,
    tempTable = FALSE,
    camelCaseToSnakeCase = TRUE
  )

  return(combinedData)
}

# Internal function: Cleanup temporary cohort
cleanupTempCohort <- function(syTrialConnection) {
  sql <- "DROP TABLE IF EXISTS @cohort_database_schema.sytrial_temp_cohort;"
  sql <- SqlRender::render(sql, cohort_database_schema = syTrialConnection$cohortDatabaseSchema)
  sql <- SqlRender::translate(sql, targetDialect = syTrialConnection$dbms)
  DatabaseConnector::executeSql(syTrialConnection$connection, sql, progressBar = FALSE, reportOverallTime = FALSE)
}
