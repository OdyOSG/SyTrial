# SyTrial: Synthetic Control Arm Construction Using OMOP CDM

## Overview

**SyTrial** is an R package for constructing synthetic control arms from Real-World Data (RWD) stored in OMOP Common Data Model (CDM) databases. It implements advanced causal inference methods including Targeted Maximum Likelihood Estimation (TMLE), G-computation, and Inverse Probability of Treatment Weighting (IPTW).

### Key Innovation: Index Date Imputation

Unlike traditional propensity score methods that require both treatment and control groups to have defined index dates, **SyTrial** implements multiple strategies to intelligently assign index dates to control cohorts that lack a natural treatment initiation date.

## Features

- **Index Date Imputation**: Risk set sampling, calendar matching, time-to-event matching, sequential trial emulation
- **High-Dimensional Covariate Extraction**: Leverages OHDSI HADES FeatureExtraction for thousands of covariates
- **Advanced Causal Inference**: TMLE with Super Learner, G-computation, IPTW
- **Flexible Weighting**: Overlap weights, IPTW, propensity score matching
- **Comprehensive Diagnostics**: Love plots, propensity score distributions, Kaplan-Meier curves
- **CRAN-Ready**: Full documentation, testing, and vignettes

## Installation

```r
# SyTrial is not (yet) published on CRAN.
# Install the development version from GitHub:
install.packages("remotes")
remotes::install_github("OdyOSG/SyTrial")
```

## Quick example (no database required)

The example below simulates a treated vs control cohort, computes propensity-score-based weights using SyTrial helpers, and runs TMLE and G-computation.

```r
set.seed(1)

# Simulate cohort-like data
n <- 400
analysisData <- data.frame(
  subjectId = seq_len(n),
  treatment = rbinom(n, 1, 0.5),
  age       = rnorm(n, mean = 60, sd = 10),
  male      = rbinom(n, 1, 0.48),
  x1        = rnorm(n),
  x2        = rbinom(n, 1, 0.3)
)

# Simulate a binary outcome with treatment effect + confounding
linpred <- with(
  analysisData,
  -2 + 0.5 * treatment + 0.03 * (age - 60) + 0.4 * male + 0.6 * x1 - 0.3 * x2
)
analysisData$outcome <- rbinom(n, 1, plogis(linpred))

covariates <- c("age", "male", "x1", "x2")

# Estimate propensity scores with a simple logistic regression
psModel <- glm(
  treatment ~ age + male + x1 + x2,
  data = analysisData,
  family = binomial()
)
analysisData$propensityScore <- as.numeric(stats::predict(psModel, type = "response"))

# Calculate weights (helpers provided by SyTrial)
analysisData$iptw <- calculateIPTW(
  propensityScores = analysisData$propensityScore,
  treatment = analysisData$treatment
)
analysisData$overlapW <- calculateOverlapWeights(
  propensityScores = analysisData$propensityScore,
  treatment = analysisData$treatment
)

# Run causal estimators
tmleResult <- performTMLE(
  data = analysisData,
  treatment = "treatment",
  outcome = "outcome",
  covariates = covariates
)

gcompResult <- performGComputation(
  data = analysisData,
  treatment = "treatment",
  outcome = "outcome",
  covariates = covariates
)

tmleResult
gcompResult
```

## Full pipeline requires OMOP CDM + cohort IDs + outcomeData

SyTrial’s full synthetic control arm workflow expects cohorts extracted from an OMOP CDM. Outcomes are commonly passed as an `outcomeData` object/data frame keyed by `subjectId` (and, depending on your design, including time-at-risk information).

Example `outcomeData` shape:

```r
# Minimal example structure (adapt as required by your design/functions)
outcomeData <- data.frame(
  subjectId = c(1, 2, 3),
  outcome   = c(0, 1, 0),
  timeAtRiskStart = as.Date(c("2020-01-01", "2020-01-01", "2020-01-01")),
  timeAtRiskEnd   = as.Date(c("2020-12-31", "2020-12-31", "2020-12-31"))
)
```

```r
## Not run:
# Connect to your OMOP CDM (details depend on your environment)
syTrialConnection <- createSyTrialConnection(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable
)

result <- runSyTrialAnalysis(
  syTrialConnection = syTrialConnection,
  treatmentCohortId = treatmentCohortId,
  controlCohortId = controlCohortId,
  outcomeData = outcomeData,
  outcomeVariable = outcomeVariable
)
## End(Not run)
```
