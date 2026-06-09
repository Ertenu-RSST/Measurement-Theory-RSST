#
# ====================================================================================
# PROGRAM NAME: Parametric GRM Refit-Bootstrap Ordering-Stability Diagnostic
# FILENAME:     RSST_GRM_refit_bootstrap_ordering_stability_diagnostic_program_v1.py
# VERSION:      v1.0 (Empirical Resampling Guard Production Suite)
# AUTHOR:       Ali Nihat Ertenu
# EMAIL:        ertenuan@msn.com
# PROFESSIONAL ROLE:  Data Analyst, IT Project Manager, and Independent Researcher
# CREDENTIALS:  B.S. in Mathematics, Middle East Technical University (1981)
#               Specialization: Applied Mathematics (Numerical Analysis Option)
#               IT Project Manager
# CORE THEORY:  RSST theory, including the Density–Granularity Law,
#               Metric Silence framework, scale-faithful Omega operator, and
#               boundary-preserving ordinal refinement architecture, developed
#               by Ali Nihat Ertenu. GRM/IRT methods and the R mirt package are
#               used as established psychometric auxiliary-ordering tools.
# EMPIRICAL DATA SOURCE: Deniz Dilara Ertenu 
#               B.A. in Psychology, Johns Hopkins University (2016)
#               M.A. in Guidance and Counseling, Bahcesehir University (2020)
# COPYRIGHT:    Copyright (c) 2026 Ali Nihat Ertenu. All Rights Reserved.
# LICENSE:      MIT Open Science License (Archived for Peer-Review Verification)
# ====================================================================================

from __future__ import annotations

import math
import os
import subprocess
import tempfile
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kendalltau

# ====================================================================================
# USER CONFIGURATION — CHANGE ONLY THESE TWO LINES FOR ROUTINE USE
# ====================================================================================
INPUT_INDICATOR = "PANAS_PA_SeverityAligned"
BOOTSTRAP_ITERATIONS = 1000
# ====================================================================================

# RSST_GRM_refit_bootstrap_ordering_stability_diagnostic_program_v1.py
# -----------------------------------------------------------------------------
# Purpose:
#   Supplementary File S3 — the RSST_GRM_refit_bootstrap_ordering_stability
#   diagnostic program — is the supplementary GRM refit-bootstrap ordering-
#   stability diagnostic for the RSST/SFORR ITEM_LEVEL_GRM pathway.
#
# Why this diagnostic program exists:
#   In the RSST/SFORR article, conditional GRM-assisted ordinal refinement
#   depends on a GRM-derived auxiliary ordering vector, A_aux. Item-target
#   consistency is an entry condition for ITEM_LEVEL_GRM, but item-target
#   consistency by itself does not establish that the GRM-derived auxiliary
#   ordering remains stable under sampling/refit variation.
#
#   This diagnostic program therefore evaluates the sampling/refit sensitivity
#   of A_aux inside native observed-score tied cells. For each bootstrap
#   iteration, cases are resampled with replacement, the graded response model
#   is refit using R mirt, and the refit model scores the original case set.
#   Within each native observed-score tied cell, pairwise ordering directions
#   are compared between the reference A_aux vector and the bootstrap/refit
#   A_aux vector.
#
# Scientific role:
#   This diagnostic program evaluates model-relative auxiliary-ordering
#   stability. It is not population-level GRM validation, not latent-truth
#   reconstruction, not scale recalibration, and not score modification.
#
# Governance boundary:
#   This Supplementary File S3 diagnostic program does NOT generate or modify GRM_RSST_Score. It does
#   NOT apply Omega (Ω), select or change K′, modify native observed scores,
#   change thresholds, alter clinical/severity categories, or override SILENT
#   governance from the production GRM-only RSST/SFORR program.
#
# Expected primary computational input:
#   The same clean target-specific item-level GRM input dataset used by
#   Supplementary File S2 — the production GRM-only RSST/SFORR program for
#   conditional Ω-based ordinal refinement and GRM_RSST_Score reporting:
#
#       GRM-Target_with_Ordinal_Items_<INPUT_INDICATOR>_Input.xlsx
#
#   Sheet: GR-Target_with_Ordinal_Items
#   Column 1: case ID
#   Column 2: bounded observed target total score
#   Columns 3+: target-specific ordinal item responses
#
# Output:
#   A diagnostic Excel workbook reporting run configuration, item-target
#   consistency, R mirt refit status, iteration-level stability, cell-level
#   stability, unstable/sensitivity-limited cells, and method notes.
# -----------------------------------------------------------------------------

# =============================================================================
# ROLE IN THE RSST/SFORR ARTICLE PIPELINE AND SUPPLEMENTARY FILE SET
# =============================================================================
# Supplementary File S3 — the RSST_GRM_refit_bootstrap_ordering_stability_diagnostic program — is used
# after the production GRM-only RSST/SFORR program has been run for the same
# target-score core and the same clean target-specific item-level input dataset.
#
# Pipeline relation:
#   Supplementary File S1 — RSST_Stage1_Diagnostics program:
#       Structural saturation / bounded-region diagnostic program. It evaluates
#       density-granularity structure, threshold/manual-region occupancy,
#       tied-cell burden, and regional heterogeneity. It does not apply Omega
#       (Ω), does not derive GRM_RSST_Score, and does not decide ITEM_LEVEL_GRM
#       refinement eligibility.
#
#   Supplementary File S2 — production GRM-only RSST/SFORR program:
#       The RSST-Scale-Faithful_Ordinal_Resolution_Refinement_SFORR GRM-only
#       program for conditional Ω-based ordinal refinement and GRM_RSST_Score
#       reporting. It uses the clean target-specific item-level input dataset
#       to evaluate ITEM_LEVEL_GRM eligibility, construct the GRM-derived
#       auxiliary ordering vector A_aux through R mirt, apply or withhold Ω,
#       apply K′ preservation governance, and report either GRM_RSST_Score or
#       the SILENT-mode preserved observed score.
#
#   Supplementary File S3 — RSST_GRM_refit_bootstrap_ordering_stability_diagnostic program:
#       This Supplementary File S3 diagnostic program uses the same clean target-specific
#       item-level input dataset as the production GRM-only RSST/SFORR program,
#       but it does not perform production ordinal refinement. Instead, it
#       evaluates whether the GRM-derived A_aux ordering is stable under
#       repeated bootstrap/refit estimation inside native observed-score tied
#       cells.
#
# Scientific question answered by this diagnostic program:
#   If the GRM is refit repeatedly on bootstrap-resampled cases, and each refit
#   model is then used to score the original case set, does the within-cell
#   ordering implied by the refit A_aux vectors materially agree with the
#   reference GRM-derived A_aux ordering?
#
# Scientific question NOT answered by this diagnostic program:
#   This diagnostic program does not prove latent-truth ordering, does not
#   validate the GRM at the population level, does not recalibrate the native
#   score, and does not determine clinical interpretation. It is a
#   model-relative auxiliary-ordering stability diagnostic only.
# =============================================================================

# =============================================================================
# PRIMARY COMPUTATIONAL INPUT — SHARED TARGET-SPECIFIC ITEM-LEVEL GRM DATASET
# =============================================================================
# Example:
#   GRM-Target_with_Ordinal_Items_MMCGI_PSB_Input.xlsx
#
# This dataset contains BOTH:
#   - the native observed bounded target total score, for example
#       MMCGI_PSB_Score
#   - the target-specific ordinal item-level responses that reproduce that
#     total score under the declared RAW_SUM scoring rule.
#
# The same clean target-specific item-level dataset is used by the production
# GRM-only RSST/SFORR program and by this diagnostic program.
#
# Role in the production GRM-only RSST/SFORR program:
#   The target score plus its item-level response structure are used for
#   conditional Ω-based ordinal refinement and GRM_RSST_Score reporting. If all
#   ITEM_LEVEL_GRM and preservation-governance conditions are satisfied, the
#   production program reports GRM_RSST_Score. If those conditions are not
#   satisfied, the production program withholds Ω and preserves the observed
#   score under SILENT mode.
#
# Role in this diagnostic program:
#   The same target score plus item-level response structure are used to refit
#   the GRM under bootstrap resampling and evaluate the refit/sampling stability
#   of the GRM-derived A_aux ordering inside native observed-score tied cells.
#
# Corresponding production output artifact:
#   Example:
#     MMCGI_PSB_Score_RSST_output_GRM-only_v1.xlsx
#
#   File-name logic:
#     MMCGI_PSB_Score
#       = native observed target-score column / target-score core.
#
#     RSST_output_GRM-only_v1
#       = production GRM-only RSST/SFORR output tag.
#
#   This output artifact contains the production output after the production
#   GRM-only RSST/SFORR run, including the observed target score,
#   GRM_RSST_Score reporting context or SILENT-mode preserved score, K′/Ω
#   governance outputs, and related audit sheets. It is not the computational
#   source for refit estimation in this diagnostic program. This diagnostic
#   program computes from the shared target-specific item-level input dataset
#   and is interpreted alongside the corresponding production output artifact
#   for traceability.
# =============================================================================

# =============================================================================
# IMPORTANT GOVERNANCE BOUNDARIES
# =============================================================================
# This Supplementary File S3 diagnostic program does NOT:
#   - apply Ω;
#   - select, optimize, or modify K′;
#   - modify native observed scores;
#   - modify thresholds, cutoffs, severity bands, or interpretive categories;
#   - modify the production GRM-only RSST/SFORR output artifact;
#   - report a new GRM_RSST_Score;
#   - override SILENT governance from the production GRM-only RSST/SFORR program;
#   - transform a sensitivity-limited tied cell into a production failure.
#
# A sensitivity-limited tied cell means reduced local confidence in the
# GRM-derived A_aux ordering for that native observed-score cell. It does not
# imply violation of native-score preservation because this diagnostic program
# does not apply Ω and does not alter the K′-bounded offsets produced by the
# production GRM-only RSST/SFORR program.
# =============================================================================

# =============================================================================
# USER CONFIGURATION NOTES
# =============================================================================
# -----------------------------------------------------------------------------
# CONFIGURATION PHILOSOPHY
# -----------------------------------------------------------------------------
# Routine use should require changing only the target indicator label and the
# requested number of bootstrap/refit iterations. All other controls are kept
# explicit for auditability and reproducibility.
#
# INPUT_INDICATOR:
#   Controls the indicator-specific input and output file names. It should match
#   the target-score core used in the shared target-specific item-level input
#   dataset, without the final "_Score" suffix when the project naming convention
#   uses files such as:
#       GRM-Target_with_Ordinal_Items_MMCGI_PSB_Input.xlsx
#
# BOOTSTRAP_ITERATIONS:
#   Major diagnostic-intensity control for this program. It changes how many GRM refit-
#   bootstrap estimations are requested and therefore how much sampling/refit
#   variation the diagnostic can reveal. It is NOT a parameter of the production GRM-only RSST/SFORR
#   conditional Ω-based ordinal refinement and GRM_RSST_Score reporting pathway and does not change Ω, K′, GRM_RSST_Score reporting,
#   native observed scores, thresholds, or SILENT governance.
#
# Practical interpretation:
#   20       smoke test / fast code-path verification
#   50       development-level diagnostic
#   100      stronger diagnostic screen
#   500      manuscript-grade final diagnostic when runtime and convergence allow
#   1000     stronger final diagnostic when runtime, convergence, and reporting
#            burden permit
#
# Larger values may reveal local ordering-sensitivity concerns that are not
# visible in smaller pilot runs. Therefore B20, B50, and B100 outputs should be
# interpreted as increasing diagnostic intensity, not as different Supplementary File S2
# production GRM-only RSST/SFORR runs.
#
# -----------------------------------------------------------------------------
# INPUT_INDICATOR drives all indicator-specific names below.
# Expected input filename convention:
#     GRM-Target_with_Ordinal_Items_<INPUT_INDICATOR>_Input.xlsx
# Output folder and filename convention:
#     <INPUT_INDICATOR>/<INPUT_INDICATOR>_GRM_refit_bootstrap_ordering_stability_B<BOOTSTRAP_ITERATIONS>.xlsx
#
# Example:
#     INPUT_INDICATOR = "PANAS_PA_SeverityAligned"
#     BOOTSTRAP_ITERATIONS = 50
# Start with 20/50/100/500. Use 1000 for final if runtime is acceptable.
#
# =============================================================================
# DERIVED FILE NAMES — DO NOT EDIT FOR ROUTINE INDICATOR CHANGES
# =============================================================================
# INPUT_FILE is derived from INPUT_INDICATOR and points to the same clean
# target-specific item-level dataset used by Supplementary File S2. The file must contain the
# GR-Target_with_Ordinal_Items sheet with:
#   column 1 = case identifier;
#   column 2 = native observed bounded target total score;
#   columns 3+ = target-specific ordinal item responses.
#
# OUTPUT_FILE is this diagnostic program's output artifact. Its B-number records
# the requested diagnostic intensity:
#   B20 = 20 requested refit-bootstrap iterations;
#   B50 = 50 requested refit-bootstrap iterations;
#   B100 = 100 requested refit-bootstrap iterations;
#   B500 = 500 requested refit-bootstrap iterations;
#   B1000 = 1000 requested refit-bootstrap iterations.
#
# The output artifact is a diagnostic record only. It is not the Supplementary File S2
# production GRM-only RSST/SFORR output artifact and it is not a replacement for the native observed score or
# the conditionally reported GRM_RSST_Score.
# -----------------------------------------------------------------------------
INPUT_FILE = fr"GRM-Target_with_Ordinal_Items_{INPUT_INDICATOR}_Input.xlsx"
INPUT_SHEET = "GR-Target_with_Ordinal_Items"

# Output workbooks are written into an indicator-specific folder so that
# B20/B50/B100/B500/B1000 diagnostic runs for the same target indicator remain
# grouped together and do not overwrite or clutter the main project directory.
OUTPUT_DIR = Path(INPUT_INDICATOR)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / f"{INPUT_INDICATOR}_GRM_refit_bootstrap_ordering_stability_B{BOOTSTRAP_ITERATIONS}.xlsx"

# -----------------------------------------------------------------------------
# R mirt / GRM REFIT MODULE ROLE
# -----------------------------------------------------------------------------
# The R mirt call is used here for the same auxiliary-ordering reason it is used in Supplementary File S2: to obtain a
# defensible GRM-derived auxiliary ordering vector from target-specific ordinal
# item responses. In this diagnostic program, however, this ordering vector is used only for stability
# diagnostics by Supplementary File S3. The theta/factor-score values are not reported as a replacement
# measurement scale and are not converted into GRM_RSST_Score by this script.
# -----------------------------------------------------------------------------
# Use the same Rscript path used in the production RSST/SFORR GRM-only program.
R_MIRT_EXECUTABLE = r"C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe"

# Refit-bootstrap controls.
BOOTSTRAP_SEED = 20260529
R_MIRT_FSCORES_METHOD = "EAP"
R_MIRT_MAX_EM_CYCLES = 1000

# -----------------------------------------------------------------------------
# STABILITY INTERPRETATION GUIDES
# -----------------------------------------------------------------------------
# These values are reporting guides for Supplementary File S3. They
# are not RSST theory-defining constants, do not determine structural saturation,
# and do not modify production GRM_RSST_Score reporting.
#
# A cell below the moderate guide is interpreted as sensitivity-limited for
# local auxiliary ordering. This is not the same as failure of native-score
# preservation. Boundary preservation remains governed by production Ω/K′/α constraints.
#
# A cell can also be labelled not_evaluable if strict finite pairwise ordering
# comparisons cannot be computed for that cell under the diagnostic rule, even
# when n_cell >= 2. Not_evaluable is not an instability label; it means that Supplementary File S3
# does not assert a pairwise ordering-stability conclusion for that cell.
# -----------------------------------------------------------------------------
# Tied-cell and stability interpretation controls.
MIN_TIED_CELL_SIZE = 2
LOW_STABILITY_FLAG = 0.70
MODERATE_STABILITY_FLAG = 0.80
STRONG_STABILITY_FLAG = 0.90

# RAW_SUM item-target consistency tolerance. Must match the production logic.
ITEM_SUM_TOLERANCE = 1e-8


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
# -----------------------------------------------------------------------------
# HELPER-FUNCTION ROLE
# -----------------------------------------------------------------------------
# The helper functions below protect the diagnostic from over-interpreting
# degenerate comparisons. Correlations are returned only when there is enough
# finite nonconstant information. Pairwise ordering agreement is computed only
# for evaluable pairs with strict nonzero ordering differences in both the
# reference and refit A_aux vectors.
# -----------------------------------------------------------------------------
def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    if len(np.unique(x[mask])) < 2 or len(np.unique(y[mask])) < 2:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).correlation)


def safe_kendall(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    if len(np.unique(x[mask])) < 2 or len(np.unique(y[mask])) < 2:
        return np.nan
    return float(kendalltau(x[mask], y[mask]).correlation)


# -----------------------------------------------------------------------------
# CORE DIAGNOSTIC METRIC: PAIRWISE ORDERING AGREEMENT
# -----------------------------------------------------------------------------
# This function is the central diagnostic comparison used by this diagnostic program.
#
# Within each native observed-score tied cell, every possible pair of cases is
# compared. The diagnostic asks whether the sign of the reference GRM-derived
# A_aux difference is preserved after bootstrap/refit estimation:
#
#   sign(reference_A_aux_i - reference_A_aux_j)
#       ==
#   sign(refit_A_aux_i - refit_A_aux_j)
#
# Agreement is therefore ordinal-directional. It does not require numeric
# equality of theta values, does not compare GRM_RSST_Score values, does not
# evaluate K′ or Ω displacement, and does not claim model-free latent-truth
# ordering. It evaluates model-relative A_aux ordering stability only.
#
# Pairs with missing or exactly equal differences are counted as ties_or_missing
# because no strict pairwise ordering direction can be evaluated for that pair.
# -----------------------------------------------------------------------------

def pairwise_agreement(reference_theta: np.ndarray, bootstrap_theta: np.ndarray, indices: list[int]) -> tuple[float, int, int, int]:
    """Return pairwise ordering agreement within one observed-score tied cell.

    A pair is evaluable when both reference and bootstrap ordering differences
    are finite and nonzero. Agreement means the pairwise order sign is identical
    in the bootstrap/refit ordering and the original/reference ordering.
    """
    stable = 0
    discordant = 0
    ties_or_missing = 0
    total_evaluable = 0

    for i, j in combinations(indices, 2):
        ref_diff = reference_theta[i] - reference_theta[j]
        boot_diff = bootstrap_theta[i] - bootstrap_theta[j]
        if not np.isfinite(ref_diff) or not np.isfinite(boot_diff):
            ties_or_missing += 1
            continue
        if ref_diff == 0 or boot_diff == 0:
            ties_or_missing += 1
            continue
        total_evaluable += 1
        if np.sign(ref_diff) == np.sign(boot_diff):
            stable += 1
        else:
            discordant += 1

    agreement = stable / total_evaluable if total_evaluable else np.nan
    return float(agreement) if np.isfinite(agreement) else np.nan, stable, discordant, ties_or_missing


def stability_label(value: float) -> str:
    if not np.isfinite(value):
        return "not_evaluable"
    if value >= STRONG_STABILITY_FLAG:
        return "strong_ordering_stability"
    if value >= MODERATE_STABILITY_FLAG:
        return "moderate_ordering_stability"
    if value >= LOW_STABILITY_FLAG:
        return "sensitivity_concern"
    return "unstable_ordering_flag"

def status_value_is_true(value) -> bool:
    """Return True for explicit truth-like convergence values."""
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


# -----------------------------------------------------------------------------
# TEMPORARY R SCRIPT WRITER
# -----------------------------------------------------------------------------
# The Python program writes a temporary R script so that all refit-bootstrap GRM
# estimations are performed through the same R mirt engine. The R script:
#   1. fits the original/reference GRM using all original cases;
#   2. computes the reference A_aux theta/factor-score vector for the original
#      case set;
#   3. loops over bootstrap-resampled case matrices;
#   4. refits the GRM on each bootstrap sample;
#   5. scores the original case set under each refit model;
#   6. writes theta_matrix.csv and mirt_refit_status.csv back to Python.
#
# This R script does not apply Ω and does not produce GRM_RSST_Score.
# -----------------------------------------------------------------------------
def write_mirt_bootstrap_r_script(path: Path) -> None:
    """Write one R script that fits the original GRM and all bootstrap/refit GRMs."""
    r_code = r'''
suppressPackageStartupMessages(library(mirt))

args <- commandArgs(trailingOnly = TRUE)
items_csv <- args[[1]]
boot_indices_csv <- args[[2]]
out_dir <- args[[3]]
fs_method <- args[[4]]
max_cycles <- as.integer(args[[5]])

dat <- read.csv(items_csv, check.names = FALSE)
boot_indices <- as.matrix(read.csv(boot_indices_csv, check.names = FALSE))
B <- nrow(boot_indices)
N <- nrow(dat)

theta_matrix <- matrix(NA_real_, nrow = N, ncol = B + 1)
colnames(theta_matrix) <- c("theta_reference", paste0("theta_boot_", seq_len(B)))

status <- data.frame(
  iteration = c(0, seq_len(B)),
  status = NA_character_,
  converged = NA,
  logLik = NA_real_,
  error_message = NA_character_,
  stringsAsFactors = FALSE
)

fit_and_score_original <- function(train_dat, original_dat) {
  fit <- mirt(
    train_dat,
    1,
    itemtype = "graded",
    method = "EM",
    verbose = FALSE,
    technical = list(NCYCLES = max_cycles)
  )
  fs <- as.data.frame(fscores(
    fit,
    method = fs_method,
    full.scores = TRUE,
    full.scores.SE = TRUE,
    response.pattern = original_dat
  ))
  score_col <- if ("F1" %in% names(fs)) "F1" else names(fs)[[1]]
  list(
    theta = as.numeric(fs[[score_col]]),
    converged = fit@OptimInfo$converged,
    logLik = as.numeric(logLik(fit))
  )
}

# Original/reference fit.
ref_result <- tryCatch(
  fit_and_score_original(dat, dat),
  error = function(e) e
)

if (inherits(ref_result, "error")) {
  status$status[1] <- "failed"
  status$error_message[1] <- conditionMessage(ref_result)
} else {
  theta_matrix[, 1] <- ref_result$theta
  reference_converged <- isTRUE(ref_result$converged)
  status$status[1] <- if (reference_converged) "reference_fit_ok" else "reference_nonconverged"
  status$converged[1] <- reference_converged
  status$logLik[1] <- ref_result$logLik
}

# Bootstrap/refit loop.
for (b in seq_len(B)) {
  idx <- as.integer(boot_indices[b, ])
  train_dat <- dat[idx, , drop = FALSE]

  result <- tryCatch(
    fit_and_score_original(train_dat, dat),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    status$status[b + 1] <- "failed"
    status$error_message[b + 1] <- conditionMessage(result)
    next
  }

  theta_matrix[, b + 1] <- result$theta
  iteration_converged <- isTRUE(result$converged)
  status$status[b + 1] <- if (iteration_converged) "ok" else "nonconverged"
  status$converged[b + 1] <- iteration_converged
  status$logLik[b + 1] <- result$logLik
}

write.csv(theta_matrix, file.path(out_dir, "theta_matrix.csv"), row.names = FALSE)
write.csv(status, file.path(out_dir, "mirt_refit_status.csv"), row.names = FALSE)
'''
    path.write_text(r_code, encoding="utf-8")


# =============================================================================
# MAIN DIAGNOSTIC WORKFLOW
# =============================================================================
# -----------------------------------------------------------------------------
# RUNTIME WORKFLOW SUMMARY
# -----------------------------------------------------------------------------
# 1. Load the shared clean target-specific item-level GRM input dataset.
# 2. Detect the case ID column, native observed target-score column, and item
#    response columns.
# 3. Verify RAW_SUM item-target consistency. If item totals do not reproduce the
#    target score, the diagnostic stops because the item-level evidence is not a
#    valid target-specific basis for GRM ordering-stability evaluation.
# 4. Identify native observed-score tied cells from the target total score.
# 5. Generate bootstrap case-resampling indices at the original sample size.
# 6. Write and run the temporary R mirt refit script.
# 7. Read the reference and bootstrap/refit A_aux theta matrices back into
#    Python.
# 8. Align refit latent direction to the reference direction when needed.
# 9. Within each tied cell, compare reference and refit pairwise ordering signs.
# 10. Export run configuration, summary, cell-level stability, iteration-level
#     stability, item-target consistency, mirt refit status, and method notes.
# -----------------------------------------------------------------------------
def main() -> None:
    input_path = Path(INPUT_FILE)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path.resolve()}")
    if not Path(R_MIRT_EXECUTABLE).exists():
        raise FileNotFoundError(
            "Rscript executable not found. Update R_MIRT_EXECUTABLE at the top of this script.\n"
            f"Current path: {R_MIRT_EXECUTABLE}"
        )

    # -------------------------------------------------------------------------
    # Load the shared target-specific item-level GRM input dataset.
    # -------------------------------------------------------------------------
    # This is the same clean item-level dataset used by the production GRM-only RSST/SFORR program. It contains the
    # native observed target total score plus the item responses that reproduce
    # that score under the declared RAW_SUM scoring rule. This diagnostic program reads it for
    # refit-bootstrap stability diagnostics, not for production ordinal
    # refinement.
    # -------------------------------------------------------------------------
    df = pd.read_excel(input_path, sheet_name=INPUT_SHEET)
    if df.shape[1] < 3:
        raise ValueError("Input sheet must contain ID, target score, and at least one item column.")

    id_col = df.columns[0]
    target_col = df.columns[1]
    item_cols = list(df.columns[2:])

    # -------------------------------------------------------------------------
    # RAW_SUM item-target consistency gate.
    # -------------------------------------------------------------------------
    # The item columns must reproduce the bounded target score under the same
    # declared scoring rule used by the production GRM-only RSST/SFORR program.
    # This protects the diagnostic from fitting a GRM to variables that are
    # merely ordinal-looking but not the actual target-specific item responses
    # for the observed score.
    #
    # If this check fails, this diagnostic program stops. That behavior is
    # intentional and consistent with the RSST/SFORR governance rule that
    # item-level role must be supported by item-target consistency before
    # GRM-derived A_aux ordering is evaluated.
    # -------------------------------------------------------------------------

    target = pd.to_numeric(df[target_col], errors="coerce")
    items = df[item_cols].apply(pd.to_numeric, errors="coerce")

    # Missingness gate:
    # Supplementary File S3 uses the same clean target-specific item-level
    # evidence required by Supplementary File S2. R mirt can estimate factor
    # scores with missing item responses, but this diagnostic does not submit
    # incomplete item-total evidence to the refit-bootstrap GRM ordering-stability
    # procedure.
    missing_item_cells = int(items.isna().sum().sum())
    rows_with_missing_items = int(items.isna().any(axis=1).sum())
    missing_target_scores = int(target.isna().sum())

    if missing_item_cells > 0 or missing_target_scores > 0:
        if missing_target_scores > 0 and missing_item_cells > 0:
            missing_reason = (
                "Supplementary File S3 diagnostic was not evaluated because "
                "missing bounded target total score values and missing item-level "
                "score components were detected."
            )
        elif missing_target_scores > 0:
            missing_reason = (
                "Supplementary File S3 diagnostic was not evaluated because "
                "one or more bounded target total score values were missing."
            )
        else:
            missing_reason = (
                "Supplementary File S3 diagnostic was not evaluated because "
                "one or more item-level score components required to reproduce "
                "the bounded target total score were missing."
            )

        raise ValueError(
            missing_reason
            + f" missing_target_scores={missing_target_scores}; "
            + f"missing_item_cells={missing_item_cells}; "
            + f"rows_with_missing_items={rows_with_missing_items}. "
            + "Incomplete item-total evidence is not passed to R mirt for the "
            + "S3 refit-bootstrap ordering-stability diagnostic."
        )

    item_sum = items.sum(axis=1, skipna=False)
    diff = item_sum - target
    valid_rows = target.notna() & item_sum.notna()
    max_abs_diff = float(np.nanmax(np.abs(diff[valid_rows]))) if valid_rows.any() else np.nan
    mean_abs_diff = float(np.nanmean(np.abs(diff[valid_rows]))) if valid_rows.any() else np.nan
    item_sum_matches_target = bool(valid_rows.any() and max_abs_diff <= ITEM_SUM_TOLERANCE)

    item_target_consistency = pd.DataFrame([
        {"field": "id_column", "value": id_col, "meaning": "Case identifier column."},
        {"field": "target_score_column", "value": target_col, "meaning": "Bounded observed target score."},
        {"field": "n_cases", "value": len(df), "meaning": "Number of input cases."},
        {"field": "n_item_columns", "value": len(item_cols), "meaning": "Number of candidate item-response columns."},
        {"field": "target_score_derivation_rule", "value": "RAW_SUM", "meaning": "Declared item-target consistency rule."},
        {"field": "missing_item_cells", "value": missing_item_cells, "meaning": "Total missing cells across item-level score component columns."},
        {"field": "rows_with_missing_items", "value": rows_with_missing_items, "meaning": "Number of cases with at least one missing item-level score component."},
        {"field": "missing_target_scores", "value": missing_target_scores, "meaning": "Number of cases with missing bounded target total score values."},
        {"field": "max_abs_difference_item_sum_minus_target", "value": max_abs_diff, "meaning": "Maximum absolute item-sum discrepancy."},
        {"field": "mean_abs_difference_item_sum_minus_target", "value": mean_abs_diff, "meaning": "Mean absolute item-sum discrepancy."},
        {"field": "item_sum_matches_target", "value": item_sum_matches_target, "meaning": "Whether item total reproduces target score within tolerance."},
    ])

    if not item_sum_matches_target:
        raise ValueError(
            "Item-target RAW_SUM consistency failed. Use the clean GRM input file for stability diagnostics. "
            f"max_abs_diff={max_abs_diff}"
        )

    # -------------------------------------------------------------------------
    # Define native observed-score tied cells from the bounded target score.
    # -------------------------------------------------------------------------
    # This diagnostic program evaluates ordering stability only inside cells where two or more cases
    # share the same native observed target score. Singleton cells do not contain
    # unresolved within-cell ordering and are therefore excluded from pairwise
    # stability evaluation.
    # -------------------------------------------------------------------------
    score_values = target.to_numpy(dtype=float)
    tied_cells = {}
    for score, idx in pd.Series(score_values).groupby(score_values).groups.items():
        indices = list(idx)
        if len(indices) >= MIN_TIED_CELL_SIZE:
            tied_cells[float(score)] = indices

    # -------------------------------------------------------------------------
    # Generate bootstrap case-resampling indices.
    # -------------------------------------------------------------------------
    # Each bootstrap sample has the original sample size and is drawn with
    # replacement. These indices control only the refit-bootstrap diagnostic.
    # They do not alter the original dataset or any production GRM-only RSST/SFORR output.
    # -------------------------------------------------------------------------
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    boot_indices = rng.integers(0, len(df), size=(BOOTSTRAP_ITERATIONS, len(df))) + 1  # R is 1-based.

    # -------------------------------------------------------------------------
    # Temporary R/mirt execution workspace.
    # -------------------------------------------------------------------------
    # Python writes the item matrix, bootstrap indices, and temporary R script to
    # a disposable folder. R mirt performs the reference and bootstrap/refit GRM
    # estimations and returns theta matrices plus refit status diagnostics.
    # -------------------------------------------------------------------------
    with tempfile.TemporaryDirectory(prefix="rsst_grm_refit_bootstrap_") as tmp:
        tmp_path = Path(tmp)
        items_csv = tmp_path / "items.csv"
        boot_csv = tmp_path / "boot_indices.csv"
        r_script = tmp_path / "run_mirt_refit_bootstrap.R"

        items.to_csv(items_csv, index=False)
        pd.DataFrame(boot_indices).to_csv(boot_csv, index=False)
        write_mirt_bootstrap_r_script(r_script)

        cmd = [
            R_MIRT_EXECUTABLE,
            str(r_script),
            str(items_csv),
            str(boot_csv),
            str(tmp_path),
            R_MIRT_FSCORES_METHOD,
            str(R_MIRT_MAX_EM_CYCLES),
        ]

        print("Running R/mirt refit-bootstrap diagnostic...")
        print(f"Input indicator: {INPUT_INDICATOR}")
        print(f"Bootstrap iterations: {BOOTSTRAP_ITERATIONS}")
        print(f"Input file: {INPUT_FILE}")
        print(f"Output file: {OUTPUT_FILE}")
        print(" ".join(cmd))
        r_env = os.environ.copy()
        r_env.setdefault("TMPDIR", r"C:\Rtmp")
        r_env.setdefault("TEMP", r"C:\Rtmp")
        r_env.setdefault("TMP", r"C:\Rtmp")
        completed = subprocess.run(cmd, capture_output=True, text=True, env=r_env)
        print(completed.stdout)
        if completed.returncode != 0:
            print(completed.stderr)
            raise RuntimeError("R/mirt refit-bootstrap diagnostic failed. See R stderr above.")

        theta_matrix = pd.read_csv(tmp_path / "theta_matrix.csv")
        refit_status = pd.read_csv(tmp_path / "mirt_refit_status.csv")

        reference_theta = theta_matrix["theta_reference"].to_numpy(dtype=float)

        refit_status["iteration"] = pd.to_numeric(refit_status["iteration"], errors="coerce").astype("Int64")
        refit_status_indexed = refit_status.set_index("iteration", drop=False)

        reference_status_row = refit_status_indexed.loc[0] if 0 in refit_status_indexed.index else None
        if reference_status_row is None:
            raise RuntimeError(
                "Reference R mirt GRM fit status was not reported; Supplementary File S3 diagnostic was not evaluated."
            )

        reference_status_text = str(reference_status_row.get("status", "")).strip().lower()
        reference_converged = status_value_is_true(reference_status_row.get("converged", False))

        if reference_status_text != "reference_fit_ok" or not reference_converged:
            raise RuntimeError(
                "Reference R mirt GRM fit did not explicitly report convergence; "
                "Supplementary File S3 diagnostic was not evaluated."
            )

        if not np.isfinite(reference_theta).all():
            raise RuntimeError(
                "Reference R mirt GRM theta vector contained non-finite values; "
                "Supplementary File S3 diagnostic was not evaluated."
            )


    iteration_rows = []
    cell_iteration_rows = []

    # -------------------------------------------------------------------------
    # Iteration-level stability evaluation.
    # -------------------------------------------------------------------------
    # For each bootstrap/refit theta vector, the script optionally reverses the
    # sign when the refit latent direction is opposite to the reference direction.
    # This protects the diagnostic from arbitrary factor-score sign indeterminacy
    # and preserves the comparison of ordinal direction rather than raw theta
    # orientation.
    # -------------------------------------------------------------------------

    for b in range(1, BOOTSTRAP_ITERATIONS + 1):
        col = f"theta_boot_{b}"

        status_row = refit_status_indexed.loc[b] if b in refit_status_indexed.index else None
        if status_row is None:
            iteration_status = "missing_refit_status"
            iteration_converged = False
        else:
            iteration_status = str(status_row.get("status", "")).strip().lower()
            iteration_converged = status_value_is_true(status_row.get("converged", False))

        boot_theta = theta_matrix[col].to_numpy(dtype=float)

        if iteration_status != "ok" or not iteration_converged or not np.isfinite(boot_theta).all():
            if iteration_status != "ok":
                exclusion_reason = f"mirt_status_{iteration_status}"
            elif not iteration_converged:
                exclusion_reason = "mirt_refit_nonconverged"
            else:
                exclusion_reason = "nonfinite_bootstrap_theta"

            iteration_rows.append({
                "iteration": b,
                "status": "excluded_nonconverged_or_failed",
                "mirt_iteration_status": iteration_status,
                "mirt_converged": iteration_converged,
                "iteration_included_in_stability_summary": False,
                "exclusion_reason": exclusion_reason,
                "direction_sign_flipped_to_match_reference": False,
                "reference_bootstrap_spearman": np.nan,
                "evaluated_tied_cells": len(tied_cells),
                "n_evaluable_pairs": 0,
                "stable_pairs": 0,
                "discordant_pairs": 0,
                "ties_or_missing_pairs": 0,
                "overall_pairwise_agreement": np.nan,
                "mean_cell_pairwise_agreement": np.nan,
                "median_cell_pairwise_agreement": np.nan,
            })
            continue

        # Align latent direction to the reference fit when needed.
        
        rho = safe_spearman(reference_theta, boot_theta)
        sign_flipped = False
        if np.isfinite(rho) and rho < 0:
            boot_theta = -boot_theta
            sign_flipped = True
            rho = safe_spearman(reference_theta, boot_theta)

        total_stable = 0
        total_discordant = 0
        total_tie_missing = 0
        total_pairs = 0
        agreements = []

        for score, indices in tied_cells.items():
            agreement, stable, discordant, tie_missing = pairwise_agreement(reference_theta, boot_theta, indices)
            n_pairs = stable + discordant
            if np.isfinite(agreement):
                agreements.append(agreement)
            total_stable += stable
            total_discordant += discordant
            total_tie_missing += tie_missing
            total_pairs += n_pairs
            cell_iteration_rows.append({
                "iteration": b,
                "observed_score_cell": score,
                "n_cell": len(indices),
                "n_evaluable_pairs": n_pairs,
                "stable_pairs": stable,
                "discordant_pairs": discordant,
                "ties_or_missing_pairs": tie_missing,
                "pairwise_agreement": agreement,
                "spearman_within_cell": safe_spearman(reference_theta[indices], boot_theta[indices]),
                "kendall_within_cell": safe_kendall(reference_theta[indices], boot_theta[indices]),
            })

        iteration_rows.append({
            "iteration": b,
            "status": "evaluated" if total_pairs > 0 else "not_evaluable",
            "mirt_iteration_status": iteration_status,
            "mirt_converged": iteration_converged,
            "iteration_included_in_stability_summary": True,
            "exclusion_reason": "not_applicable",
            "direction_sign_flipped_to_match_reference": sign_flipped,
            "reference_bootstrap_spearman": rho,
            "evaluated_tied_cells": len(tied_cells),
            "n_evaluable_pairs": total_pairs,
            "stable_pairs": total_stable,
            "discordant_pairs": total_discordant,
            "ties_or_missing_pairs": total_tie_missing,
            "overall_pairwise_agreement": total_stable / total_pairs if total_pairs else np.nan,
            "mean_cell_pairwise_agreement": float(np.nanmean(agreements)) if agreements else np.nan,
            "median_cell_pairwise_agreement": float(np.nanmedian(agreements)) if agreements else np.nan,
        })

    iteration_df = pd.DataFrame(iteration_rows)
    cell_iter_df = pd.DataFrame(cell_iteration_rows)


    if cell_iter_df.empty:
        cell_iter_df = pd.DataFrame(columns=[
            "iteration",
            "observed_score_cell",
            "n_cell",
            "n_evaluable_pairs",
            "stable_pairs",
            "discordant_pairs",
            "ties_or_missing_pairs",
            "pairwise_agreement",
            "spearman_within_cell",
            "kendall_within_cell",
        ])


    # -------------------------------------------------------------------------
    # Cell-level summary aggregated over bootstrap/refit iterations.
    # -------------------------------------------------------------------------
    # Each row summarizes one native observed-score tied cell. The median and
    # mean pairwise agreement values describe how consistently the refit GRM
    # preserves the reference A_aux ordering inside that cell.
    #
    # Interpretation:
    #   strong/moderate labels support stronger local ordering stability;
    #   sensitivity_concern indicates reduced local confidence in A_aux ordering;
    #   not_evaluable means strict pairwise comparison was not available under
    #   the diagnostic rule.
    # -------------------------------------------------------------------------
    cell_summary_rows = []
    for score, indices in tied_cells.items():
        tmp = cell_iter_df[cell_iter_df["observed_score_cell"] == score]
        vals = tmp["pairwise_agreement"].to_numpy(dtype=float)
        cell_summary_rows.append({
            "observed_score_cell": score,
            "n_cell": len(indices),
            "n_pairs_per_iteration_nominal": len(indices) * (len(indices) - 1) // 2,
            "mean_pairwise_agreement": float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan,
            "median_pairwise_agreement": float(np.nanmedian(vals)) if np.isfinite(vals).any() else np.nan,
            "min_pairwise_agreement": float(np.nanmin(vals)) if np.isfinite(vals).any() else np.nan,
            "max_pairwise_agreement": float(np.nanmax(vals)) if np.isfinite(vals).any() else np.nan,
            "mean_spearman_within_cell": float(np.nanmean(tmp["spearman_within_cell"])) if np.isfinite(tmp["spearman_within_cell"]).any() else np.nan,
            "mean_kendall_within_cell": float(np.nanmean(tmp["kendall_within_cell"])) if np.isfinite(tmp["kendall_within_cell"]).any() else np.nan,
            "stability_label_by_median": stability_label(float(np.nanmedian(vals)) if np.isfinite(vals).any() else np.nan),
        })
    cell_summary_df = pd.DataFrame(cell_summary_rows).sort_values(["median_pairwise_agreement", "observed_score_cell"], ascending=[True, True])

    unstable_cells = cell_summary_df[cell_summary_df["median_pairwise_agreement"] < MODERATE_STABILITY_FLAG].copy()

    ok_status = refit_status[refit_status["iteration"] > 0].copy()

    if "status" in ok_status.columns:
        status_lower = ok_status["status"].astype(str).str.lower()
        n_ok = int((status_lower == "ok").sum())
        n_failed = int((status_lower == "failed").sum())
        n_nonconverged = int((status_lower == "nonconverged").sum())
    else:
        n_ok = 0
        n_failed = 0
        n_nonconverged = 0

    if "converged" in ok_status.columns:
        n_converged = int(ok_status["converged"].apply(status_value_is_true).sum())
    else:
        n_converged = np.nan

    included_iteration_df = iteration_df[
        iteration_df["iteration_included_in_stability_summary"] == True
    ].copy()

    n_included_iterations = int(included_iteration_df.shape[0])
    n_excluded_iterations = int(BOOTSTRAP_ITERATIONS - n_included_iterations)

    overall_values = included_iteration_df["overall_pairwise_agreement"].to_numpy(dtype=float)
    overall_values = overall_values[np.isfinite(overall_values)]

    median_overall_pairwise_agreement = float(np.nanmedian(overall_values)) if overall_values.size else np.nan
    mean_overall_pairwise_agreement = float(np.nanmean(overall_values)) if overall_values.size else np.nan
    min_overall_pairwise_agreement = float(np.nanmin(overall_values)) if overall_values.size else np.nan
    # -------------------------------------------------------------------------
    # Run-level bootstrap summary.
    # -------------------------------------------------------------------------
    # This sheet summarizes requested and completed refits, convergence status,
    # evaluated tied cells, overall pairwise agreement, and cells below the
    # reporting guides. It is reviewer-facing diagnostic documentation, not a
    # production score sheet.
    # -------------------------------------------------------------------------
    summary_rows = [
        {"field": "indicator", "value": INPUT_INDICATOR, "meaning": "Target indicator evaluated."},
        {"field": "n_cases", "value": len(df), "meaning": "Number of original cases scored in each refit."},
        {"field": "n_items", "value": len(item_cols), "meaning": "Number of item-level ordinal response columns."},
        {"field": "bootstrap_iterations_requested", "value": BOOTSTRAP_ITERATIONS, "meaning": "Number of refit-bootstrap iterations requested."},
        {"field": "bootstrap_iterations_ok", "value": n_ok, "meaning": "Number of bootstrap/refit iterations completed successfully and explicitly reported as converged by mirt."},
        {"field": "bootstrap_iterations_failed", "value": n_failed, "meaning": "Number of bootstrap/refit iterations that failed."},
        {"field": "bootstrap_iterations_nonconverged", "value": n_nonconverged, "meaning": "Number of bootstrap/refit iterations that completed but did not explicitly report convergence."},
        {"field": "bootstrap_iterations_converged", "value": n_converged, "meaning": "Number of bootstrap/refit iterations reported as converged by mirt."},
        {"field": "bootstrap_iterations_included_in_stability_summary", "value": n_included_iterations, "meaning": "Number of refit-bootstrap iterations included in ordering-stability calculations."},
        {"field": "bootstrap_iterations_excluded_from_stability_summary", "value": n_excluded_iterations, "meaning": "Number of failed, nonconverged, or nonfinite-theta iterations excluded from ordering-stability calculations."},
        {"field": "GRM_engine", "value": "R mirt", "meaning": "GRM refit engine."},
        {"field": "factor_score_method", "value": R_MIRT_FSCORES_METHOD, "meaning": "Factor-score method used for A_aux ordering."},
        {"field": "evaluated_tied_cells", "value": len(tied_cells), "meaning": "Observed-score cells with n >= 2."},
        {"field": "median_overall_pairwise_agreement", "value": median_overall_pairwise_agreement, "meaning": "Median pairwise ordering agreement across included bootstrap/refit iterations."},
        {"field": "mean_overall_pairwise_agreement", "value": mean_overall_pairwise_agreement, "meaning": "Mean pairwise ordering agreement across included bootstrap/refit iterations."},
        {"field": "min_overall_pairwise_agreement", "value": min_overall_pairwise_agreement, "meaning": "Minimum overall pairwise agreement across included iterations."},
        {"field": "cells_below_0.90_median", "value": int((cell_summary_df["median_pairwise_agreement"] < STRONG_STABILITY_FLAG).sum()), "meaning": "Cells below strong-stability guide."},
        {"field": "cells_below_0.80_median", "value": int((cell_summary_df["median_pairwise_agreement"] < MODERATE_STABILITY_FLAG).sum()), "meaning": "Cells below moderate-stability guide."},
        {"field": "cells_below_0.70_median", "value": int((cell_summary_df["median_pairwise_agreement"] < LOW_STABILITY_FLAG).sum()), "meaning": "Cells below low-stability guide."},
        {"field": "diagnostic_role", "value": "A_aux ordering stability sensitivity check only", "meaning": "This diagnostic does not alter Omega, K-prime, observed scores, or GRM_RSST_Score generation."},
    ]
    bootstrap_summary = pd.DataFrame(summary_rows)

    # -------------------------------------------------------------------------
    # Full reproducibility record for the diagnostic run.
    # -------------------------------------------------------------------------
    # The run_configuration sheet records the actual input file, input sheet,
    # target-score column, item columns, bootstrap iteration count, seed, R mirt
    # configuration, and governance role of the diagnostic.
    #
    # Note: existing output-note wording is preserved in this comment-only
    # revision. Conceptually, the diagnostic is outside the production GRM-only conditional
    # Ω-based ordinal refinement / GRM_RSST_Score reporting pathway.
    # -------------------------------------------------------------------------
    # Full reproducibility record for the diagnostic run. This sheet is intentionally
    # separate from bootstrap_summary so reviewers can see the actual input and
    # configuration values that generated the workbook.
    run_configuration = pd.DataFrame([
        {"field": "diagnostic_name", "value": "GRM refit-bootstrap ordering-stability diagnostic", "meaning": "Supplementary diagnostic type."},
        {"field": "input_indicator", "value": INPUT_INDICATOR, "meaning": "Single user-facing indicator label controlling derived filenames and output labels."},
        {"field": "input_file", "value": INPUT_FILE, "meaning": "Source clean target-specific item-level input dataset used by Supplementary File S3."},
        {"field": "input_sheet", "value": INPUT_SHEET, "meaning": "Input sheet read from the clean target-specific item-level input dataset."},
        {"field": "output_file", "value": OUTPUT_FILE, "meaning": "Diagnostic output artifact generated by Supplementary File S3."},
        {"field": "id_column", "value": id_col, "meaning": "Case identifier column detected from the input dataset."},
        {"field": "target_score_column", "value": target_col, "meaning": "Bounded observed target-score column detected from the input dataset."},
        {"field": "item_columns_used", "value": ", ".join(map(str, item_cols)), "meaning": "Ordinal item-response columns submitted to R mirt refit diagnostics."},
        {"field": "n_cases", "value": len(df), "meaning": "Number of original cases scored after each refit."},
        {"field": "n_items", "value": len(item_cols), "meaning": "Number of item-level ordinal response columns."},
        {"field": "item_target_rule", "value": "RAW_SUM", "meaning": "Declared item-target consistency rule for this diagnostic."},
        {"field": "item_sum_matches_target", "value": bool(max_abs_diff <= ITEM_SUM_TOLERANCE), "meaning": "Whether item sums reproduced the target score within tolerance."},
        {"field": "item_sum_tolerance", "value": ITEM_SUM_TOLERANCE, "meaning": "Tolerance used for RAW_SUM item-target consistency."},
        {"field": "bootstrap_iterations_requested", "value": BOOTSTRAP_ITERATIONS, "meaning": "Number of refit-bootstrap iterations requested."},
        {"field": "bootstrap_seed", "value": BOOTSTRAP_SEED, "meaning": "Random seed used to generate bootstrap case-resampling indices."},
        {"field": "bootstrap_sample_size", "value": len(df), "meaning": "Number of cases sampled with replacement for each refit."},
        {"field": "bootstrap_sampling_rule", "value": "case resampling with replacement", "meaning": "Bootstrap sampling design."},
        {"field": "GRM_engine", "value": "R mirt", "meaning": "GRM refit engine."},
        {"field": "R_MIRT_EXECUTABLE", "value": R_MIRT_EXECUTABLE, "meaning": "Rscript executable path used to run mirt."},
        {"field": "factor_score_method", "value": R_MIRT_FSCORES_METHOD, "meaning": "Factor-score method used for A_aux ordering."},
        {"field": "max_EM_cycles", "value": R_MIRT_MAX_EM_CYCLES, "meaning": "Maximum EM cycles supplied to mirt technical controls."},
        {"field": "scoring_target", "value": "original cases", "meaning": "Each refit model is used to score the original case set for ordering comparison."},
        {"field": "diagnostic_role", "value": "A_aux ordering stability sensitivity check only", "meaning": "This diagnostic is not part of the Supplementary File S2 conditional Omega-based ordinal refinement and GRM_RSST_Score reporting pathway."},
        {"field": "modifies_observed_score", "value": False, "meaning": "The diagnostic does not modify native observed scores."},
        {"field": "modifies_K_prime", "value": False, "meaning": "The diagnostic does not modify K-prime."},
        {"field": "modifies_Omega", "value": False, "meaning": "The diagnostic does not modify Omega."},
        {"field": "generates_GRM_RSST_Score", "value": False, "meaning": "The diagnostic does not generate GRM_RSST_Score."},
    ])

    # -------------------------------------------------------------------------
    # Method notes exported into the diagnostic workbook.
    # -------------------------------------------------------------------------
    # These notes provide a compact interpretation record for spreadsheet users.
    # The full methodological explanation is carried in the comments above and
    # in the RSST/SFORR manuscript's Optional Supplement B. This comment-only
    # revision does not change the existing executable output strings.
    # -------------------------------------------------------------------------
    methods_note = pd.DataFrame([
        {"topic": "purpose", "note": "This Supplementary File S3 diagnostic output reports a supplementary GRM refit-bootstrap diagnostic for within-cell auxiliary-ordering stability."},
        {"topic": "not_score_generation", "note": "The diagnostic is not part of the RSST_Score generation pathway and does not change Omega, K-prime, observed scores, or GRM_RSST_Score."},
        {"topic": "reference_ordering", "note": "The original/reference GRM fit supplies the baseline A_aux ordering."},
        {"topic": "bootstrap_refit", "note": "Each bootstrap iteration fits a new R mirt graded response model on a case-resampled item-response matrix and scores the original cases."},
        {"topic": "pairwise_agreement", "note": "Within each native observed-score tied cell, all case pairs are compared. Agreement means the bootstrap/refit A_aux order matches the reference A_aux order."},
        {"topic": "interpretation", "note": "Higher agreement supports stronger model-relative stability of the auxiliary ordering. Low agreement flags sensitivity concerns rather than invalidating the native observed score."},
        {"topic": "governance", "note": "Cells with unstable ordering should be flagged for sensitivity review rather than interpreted as model-free latent-truth reconstruction."},
    ])

    # -------------------------------------------------------------------------
    # Export this diagnostic program's workbook.
    # -------------------------------------------------------------------------
    # Sheet roles:
    #   run_configuration       = reproducibility and diagnostic settings.
    #   bootstrap_summary       = run-level stability summary.
    #   cell_level_stability    = tied-cell-specific ordering stability.
    #   iteration_level_stability = bootstrap/refit iteration summaries.
    #   unstable_cells          = cells below the moderate-stability guide.
    #   item_target_consistency = RAW_SUM target-item consistency audit.
    #   mirt_refit_status       = R mirt fit/convergence/error status.
    #   methods_note            = compact user-facing interpretation notes.
    #
    # None of these sheets is a production GRM_RSST_Score sheet.
    # -------------------------------------------------------------------------
    with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
        run_configuration.to_excel(writer, sheet_name="run_configuration", index=False)
        bootstrap_summary.to_excel(writer, sheet_name="bootstrap_summary", index=False)
        cell_summary_df.to_excel(writer, sheet_name="cell_level_stability", index=False)
        iteration_df.to_excel(writer, sheet_name="iteration_level_stability", index=False)
        unstable_cells.to_excel(writer, sheet_name="unstable_cells", index=False)
        item_target_consistency.to_excel(writer, sheet_name="item_target_consistency", index=False)
        refit_status.to_excel(writer, sheet_name="mirt_refit_status", index=False)
        methods_note.to_excel(writer, sheet_name="methods_note", index=False)

    print("Done.")
    print(f"Output workbook: {Path(OUTPUT_FILE).resolve()}")
    print(f"Bootstrap iterations included in stability summary: {n_included_iterations}/{BOOTSTRAP_ITERATIONS}")

    if np.isfinite(median_overall_pairwise_agreement):
        print(f"Median overall pairwise agreement: {median_overall_pairwise_agreement:.4f}")
    else:
        print("Median overall pairwise agreement: not_evaluable")

    if np.isfinite(mean_overall_pairwise_agreement):
        print(f"Mean overall pairwise agreement: {mean_overall_pairwise_agreement:.4f}")
    else:
        print("Mean overall pairwise agreement: not_evaluable")

    print(f"Cells below 0.80 median agreement: {int((cell_summary_df['median_pairwise_agreement'] < MODERATE_STABILITY_FLAG).sum())}")

    print("=" * 72)
    print("************ End of Program. ************")
    print("=" * 72)

if __name__ == "__main__":
    main()
