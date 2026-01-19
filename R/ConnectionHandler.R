#' Create Database Connection for SyTrial Analysis
#'
#' @description
#' Establishes connection to OMOP CDM database and validates schema structure.
#'
#' @param connectionDetails DatabaseConnector connection details object
#' @param cdmDatabaseSchema Schema name containing OMOP CDM tables
#' @param cohortDatabaseSchema Schema name for cohort tables
#' @param cohortTable Name of the cohort table
#' @param tempEmulationSchema Schema for temporary tables (Oracle/BigQuery)
#'
#' @return A SyTrialConnection object
#' @export
#'
#' @examples
#' \dontrun{
#' connectionDetails <- DatabaseConnector::createConnectionDetails(
#'   dbms = "postgresql",
#'   server = "localhost/ohdsi",
#'   user = "postgres",
#'   password = "password"
#' )
#'
#' syTrialConn <- createSyTrialConnection(
#'   connectionDetails = connectionDetails,
#'   cdmDatabaseSchema = "cdm_synpuf",
#'   cohortDatabaseSchema = "results",
#'   cohortTable = "cohort"
#' )
#' }
createSyTrialConnection <- function(connectionDetails,
                                    cdmDatabaseSchema,
                                    cohortDatabaseSchema,
                                    cohortTable = "cohort",
                                    tempEmulationSchema = NULL) {

  # Input validation
  checkmate::assertClass(connectionDetails, "ConnectionDetails")
  checkmate::assertCharacter(cdmDatabaseSchema, len = 1)
  checkmate::assertCharacter(cohortDatabaseSchema, len = 1)
  checkmate::assertCharacter(cohortTable, len = 1)

  message("Establishing database connection...")

  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)

  message("Validating OMOP CDM structure...")

  requiredTables <- c("person", "observation_period", "condition_occurrence",
                      "drug_exposure", "procedure_occurrence", "visit_occurrence")

  for (table in requiredTables) {
    tableExists <- DatabaseConnector::existsTable(
      connection = connection,
      databaseSchema = cdmDatabaseSchema,
      tableName = table
    )

    if (!tableExists) {
      DatabaseConnector::disconnect(connection)
      stop(sprintf("Required OMOP CDM table '%s' not found in schema '%s'",
                   table, cdmDatabaseSchema))
    }
  }

  message("Connection established and validated successfully.")

  conn <- structure(
    list(
      connection = connection,
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      cohortDatabaseSchema = cohortDatabaseSchema,
      cohortTable = cohortTable,
      tempEmulationSchema = tempEmulationSchema,
      dbms = connectionDetails$dbms
    ),
    class = "SyTrialConnection"
  )

  on.exit(NULL, add = FALSE)
  conn
}

#' Disconnect SyTrial Connection
#'
#' @param syTrialConnection A SyTrialConnection object
#' @export
disconnectSyTrial <- function(syTrialConnection) {
  checkmate::assertClass(syTrialConnection, "SyTrialConnection")

  try(DatabaseConnector::disconnect(syTrialConnection$connection), silent = TRUE)
  message("Database connection closed.")
  invisible(TRUE)
}
