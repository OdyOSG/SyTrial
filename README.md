## 5. README.md

# SyTrial: Synthetic Control Arm Construction Using OMOP CDM

[![CRAN status](https://www.r-pkg.org/badges/version/SyTrial)](https://CRAN.R-project.org/package=SyTrial)
[![R build status](https://github.com/OHDSI/SyTrial/workflows/R-CMD-check/badge.svg)](https://github.com/OHDSI/SyTrial/actions)

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
# From CRAN (once published)
install.packages("SyTrial")

# Development version from GitHub
devtools::install_github("OHDSI/SyTrial")
```
