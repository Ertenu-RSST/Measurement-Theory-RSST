# ==============================================================================
# R BACKEND SETUP AND VERIFICATION FOR RSST / SFORR GRM PATHWAYS
# ==============================================================================
# Purpose:
#   This helper script installs or verifies the R package required by the
#   GRM-based RSST/SFORR pathways used in Supplementary Files S2 and S3.
#
# Tested environment:
#   R 4.5.3
#
# Required R package:
#   mirt
#
# Note:
#   This script is a setup and verification helper. It does not create a locked
#   R package snapshot and does not replace the archived Python/R source code.
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("mirt", quietly = TRUE)) {
    message("The R package 'mirt' was not detected. Installing from CRAN...")
    install.packages("mirt")
} else {
    message("The R package 'mirt' is already installed.")
}

suppressPackageStartupMessages(library(mirt))

message("R mirt backend verification completed successfully.")
message("The R environment is ready for the RSST/SFORR S2 and S3 GRM pathways.")
