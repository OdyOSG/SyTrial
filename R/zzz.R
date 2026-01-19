#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    sprintf(
      "SyTrial v%s - Synthetic Control Arm Construction\n",
      packageVersion("SyTrial")
    ),
    "Developed for OHDSI HADES Framework\n",
    "Author: Alexander Aleksiayuk (alexander_aleksiayuk@epam.com)\n",
    "\n",
    "Key Features:\n",
    "  - Index date imputation for control cohorts\n",
    "  - High-dimensional covariate extraction\n",
    "  - TMLE and G-computation for causal inference\n",
    "  - Comprehensive diagnostic tools\n",
    "\n",
    "Get started: vignette('SyTrial_Tutorial')\n",
    "Report issues: https://github.com/OHDSI/SyTrial/issues"
  )
}
