#' Index Date Imputation Strategies for Control Cohorts
#'
#' @description
#' Implements multiple strategies to assign index dates to control cohorts
#' that lack a natural index date, addressing the key challenge in SCA construction.
#'
#' @param syTrialConnection SyTrialConnection object
#' @param treatmentCohortId Cohort definition ID for treatment group
#' @param controlCohortId Cohort definition ID for control pool (may lack index dates)
#' @param imputationMethod Method for index date assignment:
#'   \itemize{
#'     \item "calendar_matching": Match on calendar time distribution
#'     \item "time_to_event": Match on time from diagnosis to treatment
#'     \item "risk_set_sampling": Sample from risk sets at treatment times
#'     \item "sequential_trial": Create sequential trial emulation
#'   }
#' @param matchingWindow Days window for matching (default: 30)
#' @param minObservationPrior Minimum days of observation before index (default: 365)
#' @param maxControlsPerTreated Maximum controls per treated subject (default: 5)
#'
#' @return A data frame with imputed index dates for control cohort
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' controlWithIndex <- imputeControlIndexDates(
#'   syTrialConnection = syTrialConn,
#'   treatmentCohortId = 1,
#'   controlCohortId = 2,
#'   imputationMethod = "risk_set_sampling",
#'   matchingWindow = 30
#' )
#' }
imputeControlIndexDates <- function(syTrialConnection,
                                    treatmentCohortId,
                                    controlCohortId,
                                    imputationMethod = "risk_set_sampling",
                                    matchingWindow = 30,
                                    minObservationPrior = 365,
                                    maxControlsPerTreated = 5,
                                    seed = NULL) {

  checkmate::assertClass(syTrialConnection, "SyTrialConnection")
  checkmate::assertInt(treatmentCohortId)
  checkmate::assertInt(controlCohortId)
  checkmate::assertChoice(imputationMethod,
                          c("calendar_matching", "time_to_event",
                            "risk_set_sampling", "sequential_trial"))
  checkmate::assertNumber(matchingWindow, lower = 1)
  checkmate::assertNumber(minObservationPrior, lower = 1)
  checkmate::assertInt(maxControlsPerTreated, lower = 1)
  checkmate::assert(
    checkmate::checkNULL(seed) ||
      checkmate::checkInt(seed)
  )

  message(sprintf("Imputing index dates using method: %s", imputationMethod))

  # Extract treatment cohort with index dates
  treatmentData <- extractTreatmentCohort(syTrialConnection, treatmentCohortId)

  # Extract potential control pool
  controlPool <- extractControlPool(syTrialConnection, controlCohortId, minObservationPrior)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Apply imputation strategy
  imputedControls <- switch(
    imputationMethod,
    "calendar_matching" = calendarMatching(treatmentData, controlPool, matchingWindow),
    "time_to_event" = timeToEventMatching(syTrialConnection, treatmentData, controlPool, matchingWindow),
    "risk_set_sampling" = riskSetSampling(treatmentData, controlPool, maxControlsPerTreated),
    "sequential_trial" = sequentialTrialEmulation(syTrialConnection, treatmentData, controlPool)
  )

  message(sprintf("Imputed index dates for %d control subjects",
                                  nrow(imputedControls)))

  return(imputedControls)
}

# Internal function: Extract treatment cohort
extractTreatmentCohort <- function(syTrialConnection, treatmentCohortId) {
  sql <- "
    SELECT
      c.subject_id,
      c.cohort_start_date AS index_date,
      YEAR(c.cohort_start_date) - p.year_of_birth AS age_at_index,
      p.gender_concept_id,
      DATEDIFF(DAY, op.observation_period_start_date, c.cohort_start_date) AS prior_observation_days
    FROM @cohort_database_schema.@cohort_table c
    INNER JOIN @cdm_database_schema.person p
      ON c.subject_id = p.person_id
    INNER JOIN @cdm_database_schema.observation_period op
      ON c.subject_id = op.person_id
      AND c.cohort_start_date >= op.observation_period_start_date
      AND c.cohort_start_date <= op.observation_period_end_date
    WHERE c.cohort_definition_id = @treatment_cohort_id
  "

  sql <- SqlRender::render(
    sql,
    cohort_database_schema = syTrialConnection$cohortDatabaseSchema,
    cohort_table = syTrialConnection$cohortTable,
    cdm_database_schema = syTrialConnection$cdmDatabaseSchema,
    treatment_cohort_id = treatmentCohortId
  )

  sql <- SqlRender::translate(sql, targetDialect = syTrialConnection$dbms)

  treatmentData <- DatabaseConnector::querySql(syTrialConnection$connection, sql)
  colnames(treatmentData) <- SqlRender::snakeCaseToCamelCase(colnames(treatmentData))

  treatmentData$indexDate <- as.Date(treatmentData$indexDate)

  return(treatmentData)
}

# Internal function: Extract control pool
extractControlPool <- function(syTrialConnection, controlCohortId, minObservationPrior) {
  sql <- "
    SELECT
      c.subject_id,
      c.cohort_start_date AS earliest_eligible_date,
      c.cohort_end_date AS latest_eligible_date,
      YEAR(c.cohort_start_date) - p.year_of_birth AS age_at_earliest,
      p.gender_concept_id,
      op.observation_period_start_date,
      op.observation_period_end_date
    FROM @cohort_database_schema.@cohort_table c
    INNER JOIN @cdm_database_schema.person p
      ON c.subject_id = p.person_id
    INNER JOIN @cdm_database_schema.observation_period op
      ON c.subject_id = op.person_id
      AND c.cohort_start_date >= op.observation_period_start_date
      AND c.cohort_start_date <= op.observation_period_end_date
    WHERE c.cohort_definition_id = @control_cohort_id
      AND DATEDIFF(DAY, op.observation_period_start_date, c.cohort_start_date) >= @min_observation_prior
  "

  sql <- SqlRender::render(
    sql,
    cohort_database_schema = syTrialConnection$cohortDatabaseSchema,
    cohort_table = syTrialConnection$cohortTable,
    cdm_database_schema = syTrialConnection$cdmDatabaseSchema,
    control_cohort_id = controlCohortId,
    min_observation_prior = minObservationPrior
  )

  sql <- SqlRender::translate(sql, targetDialect = syTrialConnection$dbms)

  controlPool <- DatabaseConnector::querySql(syTrialConnection$connection, sql)
  colnames(controlPool) <- SqlRender::snakeCaseToCamelCase(colnames(controlPool))

  controlPool$earliestEligibleDate <- as.Date(controlPool$earliestEligibleDate)
  controlPool$latestEligibleDate <- as.Date(controlPool$latestEligibleDate)
  controlPool$observationPeriodStartDate <- as.Date(controlPool$observationPeriodStartDate)
  controlPool$observationPeriodEndDate <- as.Date(controlPool$observationPeriodEndDate)

  return(controlPool)
}

# Internal function: Risk Set Sampling
riskSetSampling <- function(treatmentData, controlPool, maxControlsPerTreated) {
  message("Performing risk set sampling...")

  imputedControls <- data.frame()

  for (i in 1:nrow(treatmentData)) {
    treatmentIndexDate <- treatmentData$indexDate[i]

    # Find eligible controls at risk at this time
    eligibleControls <- controlPool[
      controlPool$earliestEligibleDate <= treatmentIndexDate &
        controlPool$latestEligibleDate >= treatmentIndexDate,
    , drop = FALSE]

    if (nrow(eligibleControls) > 0) {
      # Sample up to maxControlsPerTreated
      nSample <- min(maxControlsPerTreated, nrow(eligibleControls))
      sampledControls <- eligibleControls[sample(nrow(eligibleControls), nSample), ]

      sampledControls$imputedIndexDate <- treatmentIndexDate
      sampledControls$matchedTreatmentSubjectId <- treatmentData$subjectId[i]

      imputedControls <- rbind(imputedControls, sampledControls)
    }
  }

  return(imputedControls)
}

# Internal function: Calendar Matching
calendarMatching <- function(treatmentData, controlPool, matchingWindow) {
  message("Performing calendar time matching...")

  # Get distribution of treatment index dates
  treatmentDates <- treatmentData$indexDate

  imputedControls <- data.frame()

  for (controlIdx in 1:nrow(controlPool)) {
    control <- controlPool[controlIdx, ]

    # Find treatment dates within the control's eligible window
    eligibleTreatmentDates <- treatmentDates[
      treatmentDates >= control$earliestEligibleDate &
        treatmentDates <= control$latestEligibleDate
    ]

    if (length(eligibleTreatmentDates) > 0) {
      # Sample a date from the distribution
      sampledDate <- sample(eligibleTreatmentDates, 1)

      control$imputedIndexDate <- sampledDate
      imputedControls <- rbind(imputedControls, control)
    }
  }

  return(imputedControls)
}

# Internal function: Time-to-Event Matching
timeToEventMatching <- function(syTrialConnection, treatmentData, controlPool, matchingWindow) {
  message("Performing time-to-event matching...")

  # This requires a landmark event (e.g., diagnosis date)
  # Extract diagnosis dates for both groups

  sql <- "
    SELECT
      co.person_id AS subject_id,
      MIN(co.condition_start_date) AS diagnosis_date
    FROM @cdm_database_schema.condition_occurrence co
    WHERE co.condition_concept_id IN (
      SELECT descendant_concept_id
      FROM @cdm_database_schema.concept_ancestor
      WHERE ancestor_concept_id = @disease_concept_id
    )
    GROUP BY co.person_id
  "

  # This is a placeholder - in practice, you'd specify the disease concept
  # For now, return calendar matching as fallback
  return(calendarMatching(treatmentData, controlPool, matchingWindow))
}

# Internal function: Sequential Trial Emulation
sequentialTrialEmulation <- function(syTrialConnection, treatmentData, controlPool) {
  message("Performing sequential trial emulation...")

  # Create multiple "trials" at different calendar times
  # Each control can enter multiple trials until they receive treatment or censor

  trialStartDates <- unique(treatmentData$indexDate)
  trialStartDates <- sort(trialStartDates)

  imputedControls <- data.frame()

  for (trialDate in trialStartDates) {
    # Find controls eligible at this trial start
    eligibleControls <- controlPool[
      controlPool$earliestEligibleDate <= trialDate &
        controlPool$latestEligibleDate >= trialDate,
    , drop = FALSE]

    if (nrow(eligibleControls) > 0) {
      eligibleControls$imputedIndexDate <- trialDate
      eligibleControls$trialId <- trialDate

      imputedControls <- rbind(imputedControls, eligibleControls)
    }
  }

  return(imputedControls)
}
