#' Generate Covariate Balance Diagnostics (Love Plot)
#'
#' @description
#' Creates Love plot showing standardized mean differences before and after
#' weighting/matching for covariate balance assessment.
#'
#' @param covariateData CovariateData object
#' @param treatment Binary treatment indicator
#' @param weights Optional weights for balance assessment
#' @param threshold SMD threshold for balance (default: 0.1)
#' @param topN Number of top covariates to display (default: 50)
#' @param fileName Optional file name to save plot
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' lovePlot <- createLovePlot(
#'   covariateData = covariateData,
#'   treatment = treatment,
#'   weights = overlapWeights,
#'   threshold = 0.1
#' )
#' print(lovePlot)
#' }
createLovePlot <- function(covariateData,
                           treatment,
                           weights = NULL,
                           threshold = 0.1,
                           topN = 50,
                           fileName = NULL) {

  checkmate::assertClass(covariateData, "CovariateData")
  checkmate::assertNumeric(treatment, lower = 0, upper = 1)

  message("Calculating covariate balance...")

  if (!requireNamespace("FeatureExtraction", quietly = TRUE)) {
    stop("Package 'FeatureExtraction' is required to convert covariates to a sparse matrix. Please install it.")
  }

  # Convert covariate data to sparse matrix (FeatureExtraction long format -> dgCMatrix)
  covMatrixSparse <- FeatureExtraction::covariateDataToSparseMatrix(covariateData)

  # For diagnostics/plotting we operate on a dense matrix. This can be memory intensive for large data.
  covMatrix <- as.matrix(covMatrixSparse)

  # Calculate SMD before weighting
  smdBefore <- calculateSMD(covMatrix, treatment, weights = NULL)

  # Calculate SMD after weighting
  if (!is.null(weights)) {
    smdAfter <- calculateSMD(covMatrix, treatment, weights = weights)
  } else {
    smdAfter <- NULL
  }

  # Prepare data for plotting: map matrix columns (covariateId) to covariate names
  covariateIds <- colnames(covMatrixSparse)
  if (is.null(covariateIds)) {
    stop("Covariate matrix has no column names (covariateIds). Cannot map covariates to names.")
  }
  covariateIds <- as.integer(covariateIds)

  covariateNameMap <- unique(covariateData$covariateRef[, c("covariateId", "covariateName")])
  covariateNameMap$covariateId <- as.integer(covariateNameMap$covariateId)

  covariateNames <- covariateNameMap$covariateName[match(covariateIds, covariateNameMap$covariateId)]
  covariateNames[is.na(covariateNames)] <- as.character(covariateIds[is.na(covariateNames)])

  balanceData <- data.frame(
    covariate = covariateNames,
    smdBefore = abs(smdBefore),
    stringsAsFactors = FALSE
  )

  if (!is.null(smdAfter)) {
    balanceData$smdAfter <- abs(smdAfter)
  }

  # Select top N covariates by SMD before
  balanceData <- balanceData[order(balanceData$smdBefore, decreasing = TRUE), ]
  balanceData <- head(balanceData, topN)

  # Reshape for plotting
  if (!is.null(smdAfter)) {
    balanceDataLong <- tidyr::pivot_longer(
      balanceData,
      cols = c("smdBefore", "smdAfter"),
      names_to = "timing",
      values_to = "smd"
    )
    balanceDataLong$timing <- factor(
      balanceDataLong$timing,
      levels = c("smdBefore", "smdAfter"),
      labels = c("Before Weighting", "After Weighting")
    )
  } else {
    balanceDataLong <- data.frame(
      covariate = balanceData$covariate,
      timing = "Before Weighting",
      smd = balanceData$smdBefore
    )
  }

  # Create Love plot
  p <- ggplot2::ggplot(balanceDataLong, ggplot2::aes(x = smd, y = reorder(covariate, smd))) +
    ggplot2::geom_point(ggplot2::aes(color = timing), size = 2) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = -threshold, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Covariate Balance Assessment (Love Plot)",
      subtitle = sprintf("Threshold: %.2f", threshold),
      x = "Absolute Standardized Mean Difference",
      y = "Covariate",
      color = "Timing"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 8)
    )

  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot = p, width = 10, height = 12, dpi = 300)
    message(sprintf("Love plot saved to: %s", fileName))
  }

  # Calculate balance summary
  if (!is.null(smdAfter)) {
    nBalancedBefore <- sum(abs(smdBefore) < threshold)
    nBalancedAfter <- sum(abs(smdAfter) < threshold)

    message(sprintf(
      "Covariate balance: Before = %d/%d (%.1f%%), After = %d/%d (%.1f%%)",
      nBalancedBefore, length(smdBefore), 100 * nBalancedBefore / length(smdBefore),
      nBalancedAfter, length(smdAfter), 100 * nBalancedAfter / length(smdAfter)
    ))
  }

  return(p)
}

#' Create Kaplan-Meier Survival Curves
#'
#' @description
#' Generates Kaplan-Meier curves comparing treatment and synthetic control arms.
#'
#' @param data Data frame with survival data
#' @param time Character; name of time variable
#' @param event Character; name of event indicator variable
#' @param treatment Character; name of treatment variable
#' @param weights Optional weights for weighted survival curves
#' @param title Plot title
#' @param fileName Optional file name to save plot
#'
#' @return A survminer ggsurvplot object
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' simData <- data.frame(
#'   time = rexp(n, rate = 0.05),
#'   event = rbinom(n, size = 1, prob = 0.8),
#'   treatment = rbinom(n, size = 1, prob = 0.5)
#' )
#' km <- createKaplanMeierPlot(
#'   data = simData,
#'   time = "time",
#'   event = "event",
#'   treatment = "treatment",
#'   title = "Simulated Kaplan-Meier Curves"
#' )
#' print(km)
createKaplanMeierPlot <- function(data,
                                  time,
                                  event,
                                  treatment,
                                  weights = NULL,
                                  title = "Kaplan-Meier Survival Curves",
                                  fileName = NULL) {

  checkmate::assertDataFrame(data)
  checkmate::assertChoice(time, colnames(data))
  checkmate::assertChoice(event, colnames(data))
  checkmate::assertChoice(treatment, colnames(data))

  message("Creating Kaplan-Meier survival curves...")

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Please install it.")
  }

  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is required. Please install it.")
  }

  # Ensure treatment is a 2-level factor so legend labels map deterministically
  data[[treatment]] <- factor(data[[treatment]], levels = c(0, 1), labels = c("Control", "Treatment"))

  # Create survival object
  survFormula <- as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time, event, treatment))

  if (!is.null(weights)) {
    fit <- survival::survfit(survFormula, data = data, weights = weights)
  } else {
    fit <- survival::survfit(survFormula, data = data)
  }

  # Create plot
  p <- survminer::ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    ggtheme = ggplot2::theme_minimal(),
    palette = c("#E7B800", "#2E9FDF"),
    title = title,
    xlab = "Time (days)",
    ylab = "Survival Probability",
    legend.title = "Group",
    legend.labs = c("Control", "Treatment")
  )

  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot = p$plot, width = 10, height = 8, dpi = 300)
    message(sprintf("Kaplan-Meier plot saved to: %s", fileName))
  }

  # Log-rank test
  logRankTest <- survival::survdiff(survFormula, data = data)
  pValue <- 1 - pchisq(logRankTest$chisq, length(logRankTest$n) - 1)

  message(sprintf("Log-rank test p-value: %.4f", pValue))

  return(p)
}

#' Calculate Standardized Mean Difference (SMD)
#'
#' @param covariates Matrix of covariates
#' @param treatment Binary treatment indicator
#' @param weights Optional weights
#'
#' @return Vector of SMDs
#' @keywords internal
calculateSMD <- function(covariates, treatment, weights = NULL) {

  if (is.null(weights)) {
    weights <- rep(1, length(treatment))
  }

  smd <- numeric(ncol(covariates))

  for (i in 1:ncol(covariates)) {
    x <- covariates[, i]

    # Weighted means
    mean1 <- sum(x[treatment == 1] * weights[treatment == 1]) / sum(weights[treatment == 1])
    mean0 <- sum(x[treatment == 0] * weights[treatment == 0]) / sum(weights[treatment == 0])

    # Weighted variances
    var1 <- sum(weights[treatment == 1] * (x[treatment == 1] - mean1)^2) / sum(weights[treatment == 1])
    var0 <- sum(weights[treatment == 0] * (x[treatment == 0] - mean0)^2) / sum(weights[treatment == 0])

    # Pooled standard deviation
    pooledSD <- sqrt((var1 + var0) / 2)

    # SMD
    if (pooledSD > 0) {
      smd[i] <- (mean1 - mean0) / pooledSD
    } else {
      smd[i] <- 0
    }
  }

  return(smd)
}

#' Generate Comprehensive Diagnostic Report
#'
#' @param syTrialResult SyTrial analysis result object
#' @param outputDir Directory to save diagnostic outputs
#'
#' @return List of diagnostic objects
#' @export
generateDiagnosticReport <- function(syTrialResult, outputDir = ".") {

  checkmate::assertClass(syTrialResult, "SyTrialResult")

  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  }
  checkmate::assertDirectoryExists(outputDir)

  message("Generating comprehensive diagnostic report...")

  diagnostics <- list()

  # 1. Love plot
  lovePlot <- createLovePlot(
    covariateData = syTrialResult$covariateData,
    treatment = syTrialResult$treatment,
    weights = syTrialResult$weights,
    fileName = file.path(outputDir, "love_plot.png")
  )
  diagnostics$lovePlot <- lovePlot

  # 2. Propensity score distribution
  psPlot <- createPropensityScorePlot(
    propensityScores = syTrialResult$propensityScores,
    treatment = syTrialResult$treatment,
    fileName = file.path(outputDir, "ps_distribution.png")
  )
  diagnostics$psPlot <- psPlot

  # 3. Kaplan-Meier curves (if survival outcome)
  if (!is.null(syTrialResult$survivalData)) {
    kmPlot <- createKaplanMeierPlot(
      data = syTrialResult$survivalData,
      time = "time",
      event = "event",
      treatment = "treatment",
      weights = syTrialResult$weights,
      fileName = file.path(outputDir, "kaplan_meier.png")
    )
    diagnostics$kmPlot <- kmPlot
  }

  # 4. Summary statistics table
  summaryTable <- createSummaryTable(syTrialResult)
  write.csv(summaryTable, file.path(outputDir, "summary_statistics.csv"), row.names = FALSE)
  diagnostics$summaryTable <- summaryTable

  message(sprintf("Diagnostic report saved to: %s", outputDir))

  return(diagnostics)
}

# Internal function: Create propensity score distribution plot
createPropensityScorePlot <- function(propensityScores, treatment, fileName = NULL) {

  plotData <- data.frame(
    ps = propensityScores,
    treatment = factor(treatment, levels = c(0, 1), labels = c("Control", "Treatment"))
  )

  p <- ggplot2::ggplot(plotData, ggplot2::aes(x = ps, fill = treatment)) +
    ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
    ggplot2::labs(
      title = "Propensity Score Distribution",
      x = "Propensity Score",
      y = "Count",
      fill = "Group"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot = p, width = 10, height = 6, dpi = 300)
  }

  return(p)
}

# Internal function: Create summary statistics table
createSummaryTable <- function(syTrialResult) {

  summaryTable <- data.frame(
    Metric = c(
      "Total Subjects",
      "Treatment Group",
      "Control Group",
      "Number of Covariates",
      "Propensity Score AUC",
      "Effective Sample Size (Weighted)",
      "Average Treatment Effect",
      "95% CI Lower",
      "95% CI Upper",
      "P-value"
    ),
    Value = c(
      syTrialResult$n,
      syTrialResult$nTreated,
      syTrialResult$nControl,
      length(syTrialResult$covariates),
      round(syTrialResult$psAUC, 3),
      round(sum(syTrialResult$weights)^2 / sum(syTrialResult$weights^2), 1),
      round(syTrialResult$ate, 4),
      round(syTrialResult$ci[1], 4),
      round(syTrialResult$ci[2], 4),
      format.pval(syTrialResult$pValue, digits = 4)
    )
  )

  return(summaryTable)
}
