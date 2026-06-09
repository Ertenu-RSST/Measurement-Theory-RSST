#
#######################################################################################
# PROGRAM NAME: RSST Structural Saturation Diagnostics Pipeline
# FILENAME:     RSST_Stage1_Diagnostics_v1.py
# VERSION:      v1.0 (Official Production Release)
# ARTICLE TITLE: Regional Score Saturation Theory (RSST): Resolution Loss and
#                Scale-Faithful Ordinal Refinement in Bounded Measurement Systems
# ARTICLE SUBTITLE: A focused item-level implementation using
#                   graded-response-model-assisted within-cell ordinal refinement
# AUTHOR:       Ali Nihat Ertenu
# EMAIL:        ertenuan@msn.com
# PROFESSIONAL ROLE:  Data Analyst, IT Project Manager, and Independent Researcher
# CREDENTIALS:  B.S. in Mathematics, Middle East Technical University (1981)
#               Specialization: Applied Mathematics (Numerical Analysis Option)
#               IT Project Manager
# CORE ROLE:    Supplementary File S1 for the RSST article. Implements Stage-1
#               structural saturation diagnostics, including density–granularity
#               assessment, Metric Silence screening, bounded-region diagnostics,
#               and descriptive regional heterogeneity analysis.
# EMPIRICAL DATA SOURCE: Deniz Dilara Ertenu 
#               B.A. in Psychology, Johns Hopkins University (2016)
#               M.A. in Guidance and Counseling, Bahcesehir University (2020)
# COPYRIGHT:    Copyright (c) 2026 Ali Nihat Ertenu. All Rights Reserved.
# LICENSE:      MIT Open Science License (Archived for Peer-Review Verification)
#######################################################################################

import matplotlib
matplotlib.use("Agg")

from pathlib import Path
import textwrap
import argparse
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, shapiro, normaltest, skew, kurtosis

# ====================================================================================
INPUT_FILE = Path("RSST-Data_Cutoff_Filled_SYN.xlsx")
OUTPUT_DIR = Path("RSST_Stage1_Diagnostics-v1")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
# ====================================================================================
#
# RSST Stage‑1 Diagnostics Pipeline:  
# This program performs a complete Stage‑1 RSST diagnostic analysis using a
# manual‑threshold‑first approach. It loads all cutoff/threshold metadata
# from the workbook, constructs manual bands, and runs indicator‑level Data‑
# Health Diagnostics, Distribution Analysis, Manual Region Structure Analysis,
# Structural Saturation Analysis (granularity, density, saturation indices,
# compression evidence), and Regional Heterogeneity Analysis (median–peak
# displacement, asymmetry, multimodality). The pipeline generates all Stage‑1
# outputs, including analytical tables, diagnostic figures, and summary notes.
# All analyses are descriptive/diagnostic only; no new clinical cutoffs are
# created, and Stage‑1 does not compute Omega (Ω) or RSST_Score.
#
# ================================================================================
# MANUAL-THRESHOLD-FIRST, RSST STRUCTURAL SATURATION REPORT
# WITH Δmin MODE AND K′ PRESERVATION-SETTING REPORTING
# ================================================================================
# Core rules:
# 1. Cutoff values and Threshold/Reference values are read from separate rows
#    in the Cutoff sheet.
# 2. No new clinical cutoffs are created.
# 3. CI, TDR, KDE peak, tie burden, outlier flags, regional summaries,
#    and RSST structural indicators are descriptive/diagnostic only.
# 4. Clinically meaningful observations are retained.
# 5. Region n >= 5 permits descriptive structural classification.
#    Region n < 5 is still reported, but saturation status is marked
#    indeterminate because the manual region is too sparse for a stable
#    structural decision.
# 6. This program classifies STRUCTURAL saturation at the indicator-region
#    level. It does not apply Omega (Ω), does not derive or report
#    GRM_RSST_Score, and does not determine refinement eligibility.
# ===============================================================================
#
# -------------------------------------------------------------------------------
# Figures
# -------------------------------------------------------------------------------
FIG1_PNG = OUTPUT_DIR / "Figure_1_Manual_Threshold_Occupancy.png"
FIG1_PDF = OUTPUT_DIR / "Figure_1_Manual_Threshold_Occupancy.pdf"

FIG2_PNG = OUTPUT_DIR / "Figure_2_Data_Health_Heatmap.png"
FIG2_PDF = OUTPUT_DIR / "Figure_2_Data_Health_Heatmap.pdf"

FIG3_PNG = OUTPUT_DIR / "Figure_3_Manual_Region_Structure_Map.png"
FIG3_PDF = OUTPUT_DIR / "Figure_3_Manual_Region_Structure_Map.pdf"

FIG4_PNG = OUTPUT_DIR / "Figure_4_Regional_Heterogeneity_Map.png"
FIG4_PDF = OUTPUT_DIR / "Figure_4_Regional_Heterogeneity_Map.pdf"

# ----------------------------------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------------------------------
TABLE1_FULL_XLSX = OUTPUT_DIR / "Table_I_Manual_Threshold_Status_Full.xlsx"
TABLE1_FULL_CSV = OUTPUT_DIR / "Table_I_Manual_Threshold_Status_Full.csv"
TABLE1_DISPLAY_XLSX = OUTPUT_DIR / "Table_I_Manual_Threshold_Status_Display.xlsx"
TABLE1_DISPLAY_CSV = OUTPUT_DIR / "Table_I_Manual_Threshold_Status_Display.csv"

TABLE2_FULL_XLSX = OUTPUT_DIR / "Table_II_Data_Health_and_Assumption_Audit_Full.xlsx"
TABLE2_FULL_CSV = OUTPUT_DIR / "Table_II_Data_Health_and_Assumption_Audit_Full.csv"
TABLE2_DISPLAY_XLSX = OUTPUT_DIR / "Table_II_Data_Health_and_Assumption_Audit_Display.xlsx"
TABLE2_DISPLAY_CSV = OUTPUT_DIR / "Table_II_Data_Health_and_Assumption_Audit_Display.csv"

TABLE3_FULL_XLSX = OUTPUT_DIR / "Table_III_Manual_Region_Structure_Full.xlsx"
TABLE3_FULL_CSV = OUTPUT_DIR / "Table_III_Manual_Region_Structure_Full.csv"
TABLE3_DISPLAY_XLSX = OUTPUT_DIR / "Table_III_Manual_Region_Structure_Display.xlsx"
TABLE3_DISPLAY_CSV = OUTPUT_DIR / "Table_III_Manual_Region_Structure_Display.csv"

TABLE4_FULL_XLSX = OUTPUT_DIR / "Table_IV_RSST_Structural_Saturation_Classification_Full.xlsx"
TABLE4_FULL_CSV = OUTPUT_DIR / "Table_IV_RSST_Structural_Saturation_Classification_Full.csv"
TABLE4_DISPLAY_XLSX = OUTPUT_DIR / "Table_IV_RSST_Structural_Saturation_Classification_Display.xlsx"
TABLE4_DISPLAY_CSV = OUTPUT_DIR / "Table_IV_RSST_Structural_Saturation_Classification_Display.csv"

TABLE5_FULL_XLSX = OUTPUT_DIR / "Table_V_Regional_Structural_Heterogeneity_Profile_Full.xlsx"
TABLE5_FULL_CSV = OUTPUT_DIR / "Table_V_Regional_Structural_Heterogeneity_Profile_Full.csv"
TABLE5_DISPLAY_XLSX = OUTPUT_DIR / "Table_V_Regional_Structural_Heterogeneity_Profile_Display.xlsx"
TABLE5_DISPLAY_CSV = OUTPUT_DIR / "Table_V_Regional_Structural_Heterogeneity_Profile_Display.csv"

# ------------------------------------------------------------
# Text outputs
# ------------------------------------------------------------
PIPELINE_NOTE_TXT = OUTPUT_DIR / "Pipeline_Methodology_Note.txt"
TABLE1_NOTE_TXT = OUTPUT_DIR / "Table_I_Note.txt"
TABLE2_NOTE_TXT = OUTPUT_DIR / "Table_II_Note.txt"
TABLE3_NOTE_TXT = OUTPUT_DIR / "Table_III_Note.txt"
TABLE4_NOTE_TXT = OUTPUT_DIR / "Table_IV_Note.txt"
EXEC_SUMMARY_TXT = OUTPUT_DIR / "Executive_Summary.txt"
RECOMMENDATIONS_TXT = OUTPUT_DIR / "Indicator_Recommendations.txt"
USER_GUIDE_TXT = OUTPUT_DIR / "RSST_Stage1_Diagnostics_User_Guide_v1.txt"
TABLE5_NOTE_TXT = OUTPUT_DIR / "Table_V_Note.txt"
HETEROGENEITY_SUMMARY_TXT = OUTPUT_DIR / "Regional_Heterogeneity_Summary.txt"

# =============================================================================================
# Manual threshold metadata
# =============================================================================================
# v5 dynamic metadata layer
# ---------------------------------------------------------------------------------------------
# The previous version used a hard-coded indicator_specs list.
# This version keeps the downstream calculations and decision logic
# intact, but builds indicator_specs dynamically from the workbook's
# Cutoff sheet.
#
# Required workbook structure:
# - A Data sheet containing one ID column plus score columns.
# - A Cutoff sheet where columns correspond to the same score columns.
# - The first Cutoff-sheet column contains row labels such as:
#       Cutoff, Threshold/Reference, Min/Max
#
# Supported Cutoff and Threshold/Reference formats:
# - Single numeric value: 65
# - Multiple values in one cell: 60/70, 60,70, 60;70, 60|70
# - Multiple rows: Cutoff 1, Cutoff 2, Threshold/Reference 1, etc.
#
# Supported Min/Max formats:
# - 20/120
# - 0/3.00
# - 0,3
# - min=0 max=3
#
# Rules preserved:
# - Cutoff values remain primary when supplied.
# - Threshold/Reference values are used only when Cutoff is empty.
# - No new clinical cutoffs are created.
# - If both Cutoff and Threshold/Reference are empty, the indicator is skipped
#   and reported in the console as: No cutoff or Threshold/Reference defined.
# - For k supplied values, k+1 bands are generated:
#       <c1, c1 to <c2, ..., >=ck
# - Elevated/reference regions are all bands at or above the first supplied value.
# =============================================================================================

DATA_SHEET_NAME = "Data"
CUTOFF_SHEET_NAME = "Cutoff"

indicator_specs = []
short_labels = {}
instrument_colors = {}
band_colors = {}


def _clean_text(x):
    """Return a compact string representation, preserving readable labels."""
    if pd.isna(x):
        return ""
    return str(x).strip()


def _normalize_label(x):
    """Normalize row labels from the Cutoff sheet for robust matching."""
    return re.sub(r"[^a-z0-9]+", "", _clean_text(x).lower())


def _is_cutoff_row(label):
    """Return True for Cutoff, Cutoff 1, Cutoff_2, etc. only."""
    norm = _normalize_label(label)
    return norm.startswith("cutoff")


def _is_threshold_reference_row(label):
    """Return True for Threshold/Reference rows, excluding Cutoff rows."""
    norm = _normalize_label(label)
    return (
        norm.startswith("thresholdreference") or
        norm.startswith("threshold") or
        norm.startswith("reference")
    ) and not norm.startswith("cutoff")


def _is_score_metric_type_row(label):
    """Optional metadata row for native score metric type, e.g., integer or decimal."""
    norm = _normalize_label(label)
    return norm in {"scoremetrictype", "metrictype", "scoretype", "metric"}


def _is_minmax_row(label):
    """Return True for Min/Max, Min Max, Range, Scale Range, etc."""
    norm = _normalize_label(label)
    return norm in {"minmax", "minimummaximum", "rangeminmax", "scalerange", "range"}


def _is_instrument_row(label):
    """Optional metadata row for instrument names."""
    norm = _normalize_label(label)
    return norm in {"instrument", "scale", "measure"}


def _is_indicator_label_row(label):
    """Optional metadata row for human-readable indicator labels."""
    norm = _normalize_label(label)
    return norm in {"indicatorlabel", "indicatorname", "label", "domain", "subscale"}


def _is_manual_note_row(label):
    """Optional metadata row for analyst-supplied manual notes."""
    norm = _normalize_label(label)
    return norm in {"manualnote", "note", "manualrule", "rule"}


def _parse_numbers_from_cell(cell):
    """
    Parse one or more numbers from a cutoff/min-max cell.

    Examples
    --------
    65              -> [65.0]
    "60/70"        -> [60.0, 70.0]
    "1.44; 1.78"   -> [1.44, 1.78]
    "min=0 max=3"  -> [0.0, 3.0]
    """
    if pd.isna(cell):
        return []
    if isinstance(cell, (int, float, np.integer, np.floating)):
        if np.isfinite(cell):
            return [float(cell)]
        return []

    text = str(cell).strip()
    if text == "":
        return []

    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    return [float(x) for x in nums]


def _unique_sorted(values):
    """Return sorted unique numeric values while protecting against tiny float artifacts."""
    out = []
    for v in values:
        if pd.notna(v) and np.isfinite(v):
            rv = round(float(v), 12)
            if rv not in out:
                out.append(rv)
    return sorted(out)


def _format_threshold_value(x):
    """Format thresholds cleanly for labels, preserving decimal cutoffs."""
    if pd.isna(x):
        return "NA"
    x = float(x)
    if abs(x - round(x)) < 1e-10:
        return str(int(round(x)))
    return (f"{x:.12g}")


def _parse_min_max(cell, observed_values=None):
    """
    Parse scale minimum and maximum.

    If Min/Max is not available, the program falls back to observed minimum
    and maximum so the code can still run, but this fallback is explicitly
    less ideal because manual scale bounds are preferred.
    """
    nums = _parse_numbers_from_cell(cell)
    if len(nums) >= 2:
        lo, hi = float(nums[0]), float(nums[1])
        if hi < lo:
            lo, hi = hi, lo
        return lo, hi

    if observed_values is not None:
        arr = pd.to_numeric(pd.Series(observed_values), errors="coerce").dropna().to_numpy(dtype=float)
        if arr.size > 0:
            return float(np.min(arr)), float(np.max(arr))

    raise ValueError(
        "Scale Min/Max could not be parsed. Add a Min/Max row to the Cutoff sheet, e.g., 20/120 or 0/3.00."
    )


def _infer_instrument_from_column(column_name):
    """Infer instrument name from the score-column prefix when no metadata row is supplied."""
    col = str(column_name)

    # Generic fallback: use the first token before _Score when possible.
    cleaned = re.sub(r"_?Score$", "", col, flags=re.IGNORECASE)
    token = cleaned.split("_")[0] if "_" in cleaned else cleaned
    return token.replace("-", " ").strip() or "Unknown instrument"


def _infer_indicator_from_column(column_name, instrument_name=None):
    """Infer a readable indicator label from the score-column name."""
    col = str(column_name)
    cleaned = re.sub(r"_?Score$", "", col, flags=re.IGNORECASE)

    # If an instrument name is available, remove a matching leading token
    # from the score-column label before generating a readable indicator name.
    if instrument_name:
        inst_token = re.sub(r"[^A-Za-z0-9]+", "", str(instrument_name)).upper()
        cleaned_token = re.sub(r"[^A-Za-z0-9]+", "", cleaned).upper()
        if inst_token and cleaned_token.startswith(inst_token):
            cleaned = re.sub(r"^[^_]+_?", "", cleaned, count=1)

    cleaned = cleaned.replace("_AVERAGE", "")
    cleaned = cleaned.replace("_", " ")
    cleaned = re.sub(r"\s+", " ", cleaned).strip()

    # Preserve common uppercase abbreviations while keeping labels readable.
    replacements = {
        "Gec": "GEC",
        "Mci": "MCI",
        "Iv": "IV",
    }
    title = cleaned.title()
    for old, new in replacements.items():
        title = re.sub(rf"\b{old}\b", new, title)
    return title or str(column_name)


def _make_short_label(instrument, indicator, column_name):
    """Create compact labels for figures without requiring a hard-coded dictionary."""
    inst = re.sub(r"[^A-Za-z0-9]+", "", str(instrument)).upper()
    inst = inst[:6] if inst else "IND"
    words = re.findall(r"[A-Za-z0-9]+", str(indicator))
    if not words:
        words = re.findall(r"[A-Za-z0-9]+", str(column_name))
    abbrev = "".join(w[0].upper() for w in words[:4])
    return f"{inst}-{abbrev}" if abbrev else str(column_name)[:12]


def _make_manual_bands(cutoffs):
    """Create manual bands from one or more official cutoffs."""
    cutoffs = _unique_sorted(cutoffs)
    if len(cutoffs) == 0:
        raise ValueError("At least one value value is required for each analyzed indicator.")

    bands = []
    first = cutoffs[0]
    bands.append({"label": f"<{_format_threshold_value(first)}", "lower": None, "upper": first})

    for lo, hi in zip(cutoffs[:-1], cutoffs[1:]):
        bands.append({
            "label": f"{_format_threshold_value(lo)} to <{_format_threshold_value(hi)}",
            "lower": lo,
            "upper": hi,
        })

    last = cutoffs[-1]
    bands.append({"label": f">={_format_threshold_value(last)}", "lower": last, "upper": None})
    return bands


def _cutoff_note_from_bands(instrument, indicator, bands):
    """Create a generic manual-note line from dynamically generated bands."""
    labels = ", ".join([b["label"] for b in bands])
    return f"Manual thresholds supplied in Cutoff sheet for {instrument} | {indicator}: {labels}."


def _get_cutoff_sheet_value(cutoff_df, row_mask_func, column):
    """Return the first non-empty value for a metadata row type and indicator column."""
    label_col = cutoff_df.columns[0]
    matches = cutoff_df[cutoff_df[label_col].apply(row_mask_func)]
    if matches.empty:
        return np.nan
    for _, r in matches.iterrows():
        val = r.get(column, np.nan)
        if _clean_text(val) != "":
            return val
    return np.nan


def _get_all_cutoffs_for_column(cutoff_df, column):
    """
    Collect boundary values for one indicator.

    Minimal publication-safety update:
    - If Cutoff row is filled, use Cutoff.
    - Else if Threshold/Reference row is filled, use Threshold/Reference.
    - Else return no values and source = "None".
    """
    label_col = cutoff_df.columns[0]

    cutoff_rows = cutoff_df[cutoff_df[label_col].apply(_is_cutoff_row)]
    cutoffs = []
    for _, r in cutoff_rows.iterrows():
        cutoffs.extend(_parse_numbers_from_cell(r.get(column, np.nan)))
    cutoffs = _unique_sorted(cutoffs)
    if len(cutoffs) > 0:
        return cutoffs, "Cutoff"

    threshold_rows = cutoff_df[cutoff_df[label_col].apply(_is_threshold_reference_row)]
    thresholds = []
    for _, r in threshold_rows.iterrows():
        thresholds.extend(_parse_numbers_from_cell(r.get(column, np.nan)))
    thresholds = _unique_sorted(thresholds)
    if len(thresholds) > 0:
        return thresholds, "Threshold/Reference"

    return [], "None"


def load_indicator_specs_from_workbook(input_file, data_sheet=DATA_SHEET_NAME, cutoff_sheet=CUTOFF_SHEET_NAME):
    """
    Load the Data sheet and dynamically construct indicator_specs from Cutoff sheet.

    Returns
    -------
    df : pandas.DataFrame
        Data sheet used for analysis.

    specs : list[dict]
        Dynamic replacement for the former hard-coded indicator_specs list.
    """
    xls = pd.ExcelFile(input_file)

    if data_sheet in xls.sheet_names:
        df = pd.read_excel(input_file, sheet_name=data_sheet)
    else:
        df = pd.read_excel(input_file, sheet_name=xls.sheet_names[0])

    if cutoff_sheet not in xls.sheet_names:
        raise ValueError(
            f"Cutoff sheet '{cutoff_sheet}' was not found. Available sheets: {xls.sheet_names}"
        )

    cutoff_df = pd.read_excel(input_file, sheet_name=cutoff_sheet)
    if cutoff_df.shape[1] < 2:
        raise ValueError("Cutoff sheet must contain one row-label column plus at least one indicator column.")

    label_col = cutoff_df.columns[0]
    score_columns = [c for c in cutoff_df.columns if c != label_col and c in df.columns]

    missing_from_data = [c for c in cutoff_df.columns if c != label_col and c not in df.columns]
    if missing_from_data:
        print("Warning: these Cutoff-sheet columns were not found in Data and were skipped:")
        for c in missing_from_data:
            print(f"- {c}")

    if len(score_columns) == 0:
        raise ValueError("No Cutoff-sheet indicator columns matched columns in the Data sheet.")

    specs = []
    for order, column in enumerate(score_columns):
        cutoffs, boundary_source = _get_all_cutoffs_for_column(cutoff_df, column)
        if len(cutoffs) == 0:
            print(f"No cutoff or Threshold/Reference defined for {column}; this column was skipped.")
            continue

        values = pd.to_numeric(df[column], errors="coerce")
        minmax_cell = _get_cutoff_sheet_value(cutoff_df, _is_minmax_row, column)
        scale_min, scale_max = _parse_min_max(minmax_cell, observed_values=values)

        instrument_cell = _get_cutoff_sheet_value(cutoff_df, _is_instrument_row, column)
        indicator_cell = _get_cutoff_sheet_value(cutoff_df, _is_indicator_label_row, column)
        manual_note_cell = _get_cutoff_sheet_value(cutoff_df, _is_manual_note_row, column)
        score_metric_type_cell = _get_cutoff_sheet_value(cutoff_df, _is_score_metric_type_row, column)

        instrument = _clean_text(instrument_cell) or _infer_instrument_from_column(column)
        indicator = _clean_text(indicator_cell) or _infer_indicator_from_column(column, instrument)

        bands = _make_manual_bands(cutoffs)
        elevated_regions = bands[1:]
        primary_cutoff = cutoffs[-1]
        secondary_cutoff = cutoffs[-2] if len(cutoffs) >= 2 else np.nan
        manual_note = _clean_text(manual_note_cell) or _cutoff_note_from_bands(instrument, indicator, bands)

        explicit_metric_type = _clean_text(score_metric_type_cell).lower()
        if explicit_metric_type in {"integer", "int", "unit", "unit-spaced", "unit_spaced"}:
            score_metric_type = "integer"
        elif explicit_metric_type in {"decimal", "fixed-precision", "fixed_precision", "average", "continuous"}:
            score_metric_type = "decimal"
        else:
            observed_numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
            has_decimal_observed = observed_numeric.size > 0 and np.any(np.abs(observed_numeric - np.round(observed_numeric)) > 1e-10)
            has_decimal_scale_bounds = any(
                pd.notna(v) and abs(float(v) - round(float(v))) > 1e-10
                for v in [scale_min, scale_max]
            )
            score_metric_type = "decimal" if (has_decimal_observed or has_decimal_scale_bounds) else "integer"

        specs.append({
            "_order": order,
            "instrument": instrument,
            "indicator": indicator,
            "column": column,
            "scale_min": float(scale_min),
            "scale_max": float(scale_max),
            "manual_bands": bands,
            "elevated_regions": elevated_regions,
            "primary_cutoff": float(primary_cutoff),
            "secondary_cutoff": float(secondary_cutoff) if pd.notna(secondary_cutoff) else np.nan,
            "manual_note": manual_note,
            "cutoff_values": cutoffs,
            "boundary_source": boundary_source,
            "score_metric_type": score_metric_type,
        })

    if len(specs) == 0:
        raise ValueError("No usable indicator specifications could be constructed from the Cutoff sheet.")

    return df, specs


def initialize_dynamic_plot_labels(specs):
    """Populate short_labels dynamically after specs are loaded."""
    global short_labels
    short_labels = {}
    for spec in specs:
        short_labels[spec["indicator"]] = _make_short_label(
            spec["instrument"],
            spec["indicator"],
            spec["column"],
        )


def get_band_color(label):
    """Return a stable color for any dynamically generated manual band label."""
    if label not in band_colors:
        cmap = plt.get_cmap("tab20")
        band_colors[label] = cmap(len(band_colors) % 20)
    return band_colors[label]


def get_instrument_color(instrument):
    """Return a stable color for any instrument name."""
    if instrument not in instrument_colors:
        cmap = plt.get_cmap("tab10")
        instrument_colors[instrument] = cmap(len(instrument_colors) % 10)
    return instrument_colors[instrument]


def short_label_for(indicator):
    """Return a compact plotting label, with safe fallback."""
    return short_labels.get(indicator, str(indicator)[:16])

# ============================================================
# Δmin / K′ operational settings
# ============================================================
# Stage 1 uses Δmin to compute regional granularity (G) and saturation
# indices (S = D/G). K′ is NOT used to classify structural saturation
# because Omega (Ω) is not applied in this Stage 1 program.
#
# K′ is reported here only as a recommended preservation setting for
# possible later Omega / RSST_Score refinement. It is not used in
# Stage 1 structural saturation classification.
#
# Implementation rule:
# - Integer-valued indicators under Δmin AUTO mode:
#   K′ = 10000
# - Decimal / fixed-precision / average-score indicators under Δmin AUTO mode:
#   K′ = 1000000
# - RSS comparison mode, when used in the refinement module:
#   K′ = 1
# ============================================================

DELTA_MIN_MODE = "AUTO"        # Options: "AUTO", "USER", "HYBRID"
DEFAULT_USER_DELTA = None      # Used only for USER or HYBRID modes

# Optional per-column override for USER / HYBRID modes.
# Leave empty for AUTO mode.
USER_DELTA_BY_COLUMN = {
    # Example:
    # "Example_Decimal_Average_Score": 0.01,
}

# ============================================================
# Helper functions
# ============================================================
def pct(n, total):
    return np.nan if total == 0 else (float(n) / float(total)) * 100.0


def format_n_pct(n, total):
    return "" if total == 0 else f"{int(n)} ({pct(n, total):.1f}%)"


def compute_kde(values, grid_min, grid_max, grid_points=4000):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return np.array([np.nan]), np.array([np.nan])

    if len(values) == 1 or np.allclose(values, values[0]):
        x = np.linspace(grid_min, grid_max, grid_points)
        y = np.zeros_like(x)
        y[np.argmin(np.abs(x - values[0]))] = 1.0
        return x, y

    kde = gaussian_kde(values)
    x = np.linspace(grid_min, grid_max, grid_points)
    y = kde(x)
    return x, y


def compute_kde_peak(values, scale_min, scale_max):
    x, y = compute_kde(values, scale_min, scale_max)
    if np.isnan(x).all():
        return np.nan
    return float(x[np.nanargmax(y)])


def tie_proportion(values):
    values = np.asarray(values)
    n = len(values)
    if n <= 1:
        return 0.0
    counts = pd.Series(values).value_counts(dropna=False)
    tied_n = counts[counts > 1].sum()
    return (tied_n / n) * 100.0


def iqr_outlier_counts(values):
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1

    mild_lo = q1 - 1.5 * iqr
    mild_hi = q3 + 1.5 * iqr
    ext_lo = q1 - 3.0 * iqr
    ext_hi = q3 + 3.0 * iqr

    mild_n = int(np.sum((values < mild_lo) | (values > mild_hi)))
    extreme_n = int(np.sum((values < ext_lo) | (values > ext_hi)))
    return mild_n, extreme_n, q1, q3, iqr


def normality_flag(shapiro_p, dagostino_p):
    shapiro_available = pd.notna(shapiro_p)
    dagostino_available = pd.notna(dagostino_p)

    shapiro_normal = shapiro_available and shapiro_p > 0.05
    dagostino_normal = dagostino_available and dagostino_p > 0.05

    shapiro_nonnormal = shapiro_available and shapiro_p <= 0.05
    dagostino_nonnormal = dagostino_available and dagostino_p <= 0.05

    if shapiro_available and dagostino_available:
        if shapiro_normal and dagostino_normal:
            return "Approximately normal (Shapiro-Wilk, D'Agostino-Pearson)"
        if shapiro_nonnormal and dagostino_nonnormal:
            return "Non-normal (Shapiro-Wilk, D'Agostino-Pearson)"
        if shapiro_nonnormal and dagostino_normal:
            return "Mixed normality results (Shapiro-Wilk non-normal, D'Agostino-Pearson approximately normal)"
        if shapiro_normal and dagostino_nonnormal:
            return "Mixed normality results (Shapiro-Wilk approximately normal, D'Agostino-Pearson non-normal)"

    if shapiro_available and not dagostino_available:
        if shapiro_normal:
            return "Approximately normal (Shapiro-Wilk only)"
        return "Non-normal (Shapiro-Wilk only)"

    if dagostino_available and not shapiro_available:
        if dagostino_normal:
            return "Approximately normal (D'Agostino-Pearson only)"
        return "Non-normal (D'Agostino-Pearson only)"

    return "Normality not assessed"


def mask_region(values, lower, upper):
    values = np.asarray(values, dtype=float)
    mask = np.ones(len(values), dtype=bool)

    if lower is not None:
        mask &= values >= lower
    if upper is not None:
        mask &= values < upper

    return mask


def density_ratio_in_region(values, scale_min, scale_max, lower, upper):
    x, y = compute_kde(values, scale_min, scale_max)

    if np.isnan(x).all():
        return np.nan

    total = np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)
    if total <= 0:
        return np.nan

    lo = scale_min if lower is None else max(scale_min, lower)
    hi = scale_max if upper is None else min(scale_max, upper)

    mask = (x >= lo) & ((x < hi) if upper is not None else (x <= hi))
    if not np.any(mask):
        return np.nan

    area = np.trapezoid(y[mask], x[mask]) if hasattr(np, "trapezoid") else np.trapz(y[mask], x[mask])
    return float(area / total)


def band_counts_for_spec(values, spec):
    counts = {}
    for band in spec["manual_bands"]:
        counts[band["label"]] = int(np.sum(mask_region(values, band["lower"], band["upper"])))
    return counts


def threshold_status_label(spec, band_counts):
    """
    Generic cutoff/reference status label.

    Cutoff and Threshold/Reference values are read from separate input rows.
    The downstream band logic is unchanged.
    """
    bands = spec["manual_bands"]
    elevated = bands[1:]
    source = spec.get("boundary_source", "Cutoff")

    occupied_elevated = [b["label"] for b in elevated if band_counts.get(b["label"], 0) > 0]
    if occupied_elevated:
        highest = occupied_elevated[-1]
        if source == "Threshold/Reference":
            return f"Threshold/Reference observations present ({highest})"
        return f"Cutoff observations present ({highest})"

    if source == "Threshold/Reference":
        return f"No Threshold/Reference elevation ({bands[0]['label']} only)"
    return f"No cutoff elevation ({bands[0]['label']} only)"


def region_pattern_label(region_n, tie_pct, occupied_width, scale_start, scale_end, unique_values):
    if region_n == 0:
        return "Empty manual region"
    if region_n < 5:
        return "Sparse manual region"

    full_width = scale_end - scale_start
    width_ratio = np.nan if full_width <= 0 else occupied_width / full_width

    if tie_pct >= 50 and unique_values <= max(3, int(region_n / 3)):
        return "Densely repeated region"

    if pd.notna(width_ratio) and width_ratio <= 0.25 and region_n >= 10:
        return "Narrow but heavily occupied region"

    if pd.notna(width_ratio) and width_ratio >= 0.60:
        return "Broadly occupied region"

    return "Structured manual region"


def region_structure(values, lower, upper, scale_min, scale_max):
    vals = np.asarray(values, dtype=float)
    vals = vals[mask_region(vals, lower, upper)]
    n = len(vals)

    if n == 0:
        return {
            "Region n": 0,
            "Observed start": np.nan,
            "Observed end": np.nan,
            "Observed width": np.nan,
            "Unique values": 0,
            "Tie proportion (%)": np.nan,
            "IQR": np.nan,
            "KDE peak": np.nan,
            "Interpretation note": "No observations in this manual region.",
        }

    q1 = np.percentile(vals, 25)
    q3 = np.percentile(vals, 75)
    reg_iqr = q3 - q1
    uniq = int(pd.Series(vals).nunique())
    ties = tie_proportion(vals)
    obs_start = float(np.min(vals))
    obs_end = float(np.max(vals))
    width = obs_end - obs_start

    if obs_end > obs_start:
        peak = compute_kde_peak(vals, obs_start, obs_end)
    else:
        peak = float(obs_start)

    full_region_start = scale_min if lower is None else lower
    full_region_end = scale_max if upper is None else upper
    full_region_width = full_region_end - full_region_start

    if n < 5:
        note = "Region contains fewer than 5 observations; descriptive interpretation is limited."
    else:
        width_ratio = np.nan if full_region_width <= 0 else width / full_region_width

        if ties >= 50 and uniq <= max(3, int(n / 3)):
            note = "Possible saturation signal: the region is heavily repeated and may display apparent homogeneity."
        elif pd.notna(width_ratio) and width_ratio <= 0.25 and n >= 10:
            note = "Possible compression signal: the region is narrow relative to its full manual width."
        else:
            note = "Region size is sufficient for descriptive reporting; interpret with manual thresholds as primary."

    return {
        "Region n": n,
        "Observed start": obs_start,
        "Observed end": obs_end,
        "Observed width": width,
        "Unique values": uniq,
        "Tie proportion (%)": ties,
        "IQR": reg_iqr,
        "KDE peak": peak,
        "Interpretation note": note,
    }


# ============================================================
# RSST structural-saturation helper functions
# ============================================================
# These helpers operationalize the RSST density-granularity logic
# at the indicator-region level. They are intentionally conservative:
# S > 1 is treated as necessary structural information but is not treated
# as sufficient evidence of saturation by itself. Structural saturation
# requires density-granularity imbalance together with observable
# compression evidence such as high tie burden, reduced unique-score
# representation, restricted occupied width, or low regional IQR.
#
# IMPORTANT:
# - These helpers do NOT apply Omega (Ω).
# - They do NOT derive or report RSST_Score.
# - They identify whether a region should move forward to a separate
#   auxiliary-structure / refinement-eligibility assessment.
# ============================================================

def safe_div(numer, denom):
    """Return numer / denom, with NaN protection for zero or missing denominators."""
    if pd.isna(numer) or pd.isna(denom) or denom == 0:
        return np.nan
    return float(numer) / float(denom)


def compute_delta_min(values, mode="AUTO", user_delta=None, tol=1e-12, rounding_decimals=12):
    """
    Compute Δmin under AUTO, USER, or HYBRID mode.

    AUTO:
        Δmin is computed empirically as the smallest positive difference
        between distinct observed values.

    USER:
        Δmin is supplied directly by the analyst.

    HYBRID:
        A user-specified Δmin is accepted, but diagnostic warnings are
        generated when it is much smaller or much larger than the empirical
        Δmin. If no positive user_delta is supplied, HYBRID falls back to
        empirical AUTO behavior and records a diagnostic note.

    Returns
    -------
    final_delta : float
        The Δmin value used for regional granularity computation.

    empirical_delta : float
        The observed empirical Δmin.

    warning_text : str
        Diagnostic warning text, if applicable.
    """

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]

    if arr.size <= 1:
        empirical_delta = 1.0
    else:
        uniq = np.sort(np.unique(arr))
        if uniq.size <= 1:
            empirical_delta = 1.0
        else:
            diffs = np.diff(uniq)
            positive = diffs[diffs > tol]
            empirical_delta = float(np.min(positive)) if positive.size > 0 else 1.0

    # Remove floating-point artifacts, e.g., 0.00999999999999979 -> 0.01.
    empirical_delta = round(float(empirical_delta), rounding_decimals)

    mode = str(mode).upper()
    warning_text = ""

    if mode == "AUTO":
        final_delta = empirical_delta

    elif mode == "USER":
        if user_delta is None or user_delta <= 0:
            raise ValueError("USER mode requires a positive user_delta.")
        final_delta = float(user_delta)

    elif mode == "HYBRID":
        if user_delta is None or user_delta <= 0:
            final_delta = empirical_delta
            warning_text = "HYBRID mode used empirical Δmin because no positive user_delta was provided."
        else:
            final_delta = float(user_delta)

            warnings = []
            if final_delta < empirical_delta / 10:
                warnings.append(
                    f"user_delta={final_delta} is much smaller than empirical Δmin={empirical_delta}"
                )
            if final_delta > empirical_delta * 10:
                warnings.append(
                    f"user_delta={final_delta} is much larger than empirical Δmin={empirical_delta}"
                )

            warning_text = " | ".join(warnings)

    else:
        raise ValueError(f"Unknown Δmin mode: {mode}")

    final_delta = round(float(final_delta), rounding_decimals)
    return final_delta, empirical_delta, warning_text


def select_k_prime(spec, delta_mode="AUTO", refinement_mode="RSST"):
    """
    Select the preservation parameter K′ for the later Ω / RSST_Score module.

    IMPORTANT:
    This Stage 1 program reports K′ but does not use K′ to classify
    structural saturation. K′ becomes operational only when Ω is applied.

    Generalized v5 rules:
    - RSS mode uses K′ = 1.
    - Integer-valued clinical score metrics under Δmin AUTO mode use K′ = 10000.
    - Decimal / fixed-precision / average-score metrics under Δmin AUTO mode
      use K′ = 1000000.

    This generalizes the earlier scale-specific decimal rule so that the
    same program can be used for integer, decimal, or fixed-precision datasets
    with the same workbook format.
    """

    delta_mode = str(delta_mode).upper()
    refinement_mode = str(refinement_mode).upper()

    if refinement_mode == "RSS":
        return 1

    if refinement_mode != "RSST":
        raise ValueError(f"Unknown refinement_mode: {refinement_mode}")

    metric_type = str(spec.get("score_metric_type", "integer")).lower()

    if metric_type == "decimal" and delta_mode == "AUTO":
        return 1000000

    return 10000


def k_prime_rule_note(spec, delta_mode="AUTO"):
    """
    Return a readable note explaining why K′ was selected.
    """

    metric_type = str(spec.get("score_metric_type", "integer")).lower()

    if metric_type == "decimal" and str(delta_mode).upper() == "AUTO":
        return "Decimal/fixed-precision score metric under Δmin AUTO mode; K′ = 1000000 for later Ω/RSST_Score preservation."

    return "Integer-valued score metric under current settings; K′ = 10000 for later Ω/RSST_Score preservation."


def count_metric_positions(start, end, delta_min, upper_exclusive=False):
    """
    Count the number of available metric positions in a bounded interval.

    Parameters
    ----------
    start : float
        Lower boundary of the region.
    end : float
        Upper boundary of the region.
    delta_min : float
        Smallest admissible metric step.
    upper_exclusive : bool
        True when the manual region is [start, end), e.g., 60 to <70.
        False when the region includes the endpoint, e.g., >=70 up to scale max.

    Notes
    -----
    For [60, 70) with Δmin = 1, this returns 10 score positions: 60-69.
    For [70, 120] with Δmin = 1, this returns 51 score positions: 70-120.
    """
    if pd.isna(start) or pd.isna(end) or pd.isna(delta_min) or delta_min <= 0:
        return np.nan

    width = float(end) - float(start)
    if width < 0:
        return np.nan

    eps = delta_min * 1e-7

    if upper_exclusive:
        # Number of k values for start + k*Δmin < end.
        # The small epsilon protects exact integer/decimal boundaries.
        return int(max(0, np.floor((width - eps) / delta_min) + 1))

    # Inclusive endpoint count.
    return int(max(1, np.floor((width + eps) / delta_min) + 1))


def classify_rsst_structural_status(region_n, s_manual, s_occupied, compression_count):
    """
    Classify RSST structural saturation status for an indicator-region.

    Decision states
    ---------------
    - Empty manual region:
      no empirical observations occupy the target manual region.
    - Indeterminate:
      region n < 5. The region is reported but too sparse for stable
      structural classification.
    - Structurally saturated:
      density-granularity imbalance is present AND at least two observable
      compression markers are present.
    - Compression present; structural saturation not fully confirmed:
      either imbalance or compression is present, but the full convergence
      condition is not met.
    - Saturation absent:
      no density-granularity imbalance and no meaningful compression pattern.
    """
    if region_n == 0:
        return "Absent - empty manual region"

    if region_n < 5:
        return "Indeterminate - sparse manual region"

    density_imbalance = (
        (pd.notna(s_manual) and s_manual > 1.0) or
        (pd.notna(s_occupied) and s_occupied > 1.0)
    )

    if density_imbalance and compression_count >= 2:
        return "Structurally saturated"

    if density_imbalance or compression_count >= 2:
        return "Compression present; structural saturation not fully confirmed"

    return "Saturation absent"


def omega_pathway_status(rsst_status):
    """
    Provide the correct Stage 1 statement about Ω.

    Stage 1 never applies Ω and never creates RSST_Score. It only states
    whether a region should proceed to auxiliary-structure assessment.
    """
    if rsst_status == "Structurally saturated":
        return "Auxiliary-structure assessment required; Ω not applied at Stage 1"

    if rsst_status == "Compression present; structural saturation not fully confirmed":
        return "Optional auxiliary assessment if substantively justified; Ω not applied at Stage 1"

    if "Indeterminate" in rsst_status:
        return "Not assessed for Ω; manual region is too sparse at Stage 1"

    return "Ω not indicated at Stage 1"

# ============================================================
# Table I: Manual threshold status
# ============================================================
def build_table1(df):
    rows = []

    for idx, spec in enumerate(indicator_specs):
        values = pd.to_numeric(df[spec["column"]], errors="coerce").dropna().to_numpy(dtype=float)
        n = len(values)
        counts = band_counts_for_spec(values, spec)
        bands = spec["manual_bands"]

        row = {
            "_order": idx,
            "Instrument": spec["instrument"],
            "Indicator": spec["indicator"],
            "Data column": spec["column"],
            "N non-missing": n,
            "Cutoff/Threshold status": threshold_status_label(spec, counts),
            "Number of supplied cutoff/threshold values": len(spec.get("cutoff_values", [])),
            "Cutoff/Threshold values": "; ".join(_format_threshold_value(v) for v in spec.get("cutoff_values", [])),
            "Manual rule": spec["manual_note"],
        }

        # Dynamic band columns allow 1, 2, 3, 4, or more manual cutoffs.
        for band_idx, band in enumerate(bands, start=1):
            label = band["label"]
            row[f"Band {band_idx} label"] = label
            row[f"Band {band_idx} n (%)"] = format_n_pct(counts[label], n)

        rows.append(row)

    out = pd.DataFrame(rows).sort_values("_order").drop(columns=["_order"]).reset_index(drop=True)
    return out


# ============================================================
# Table II: Data health and assumption audit
# ============================================================
def build_table2(df):
    rows = []

    for idx, spec in enumerate(indicator_specs):
        values_all = pd.to_numeric(df[spec["column"]], errors="coerce")
        values = values_all.dropna().to_numpy(dtype=float)

        n_total = len(values_all)
        n = len(values)
        n_missing = n_total - n

        scale_min = float(spec["scale_min"])
        scale_max = float(spec["scale_max"])
        scale_range = scale_max - scale_min
        primary = float(spec["primary_cutoff"])
        secondary = spec["secondary_cutoff"]

        mean_val = float(np.mean(values))
        median_val = float(np.median(values))
        sd_val = float(np.std(values, ddof=1)) if n > 1 else 0.0
        min_obs = float(np.min(values))
        max_obs = float(np.max(values))

        mild_n, extreme_n, q1, q3, iqr = iqr_outlier_counts(values)

        skew_val = float(skew(values, bias=False)) if n >= 3 else np.nan
        kurt_val = float(kurtosis(values, fisher=True, bias=False)) if n >= 4 else np.nan
        shapiro_p = float(shapiro(values)[1]) if 3 <= n <= 5000 else np.nan
        dagostino_p = float(normaltest(values)[1]) if n >= 8 else np.nan

        ties_pct = tie_proportion(values)
        unique_n = int(pd.Series(values).nunique())

        kde_peak = compute_kde_peak(values, scale_min, scale_max)
        ci = abs(mean_val - kde_peak) / scale_range if scale_range > 0 else np.nan

        tdr_primary = density_ratio_in_region(values, scale_min, scale_max, primary, None)
        tdr_secondary = (
            density_ratio_in_region(values, scale_min, scale_max, secondary, primary)
            if pd.notna(secondary) else np.nan
        )

        floor_near = pct(np.sum(values <= scale_min + 0.05 * scale_range), n)
        ceiling_near = pct(np.sum(values >= scale_max - 0.05 * scale_range), n)

        rows.append({
            "_order": idx,
            "Instrument": spec["instrument"],
            "Indicator": spec["indicator"],
            "N total": n_total,
            "N non-missing": n,
            "N missing": n_missing,
            "% missing": pct(n_missing, n_total),

            "Scale minimum": scale_min,
            "Scale maximum": scale_max,
            "Observed minimum": min_obs,
            "Observed maximum": max_obs,

            "Mean": mean_val,
            "Median": median_val,
            "SD": sd_val,
            "Q1": q1,
            "Q3": q3,
            "IQR": iqr,

            "Skewness": skew_val,
            "Kurtosis": kurt_val,
            "Shapiro-Wilk p": shapiro_p,
            "D'Agostino-Pearson p": dagostino_p,
            "Normality flag": normality_flag(shapiro_p, dagostino_p),

            "Unique values": unique_n,
            "Tie proportion (%)": ties_pct,

            "IQR mild outliers n": mild_n,
            "IQR mild outliers %": pct(mild_n, n),
            "IQR extreme outliers n": extreme_n,
            "IQR extreme outliers %": pct(extreme_n, n),

            "Primary cutoff": primary,
            "Secondary cutoff": secondary,
            "% >= primary cutoff": pct(np.sum(values >= primary), n),
            "% in secondary manual region": (
                pct(np.sum(mask_region(values, secondary, primary)), n)
                if pd.notna(secondary) else np.nan
            ),

            "Floor-near occupancy (%)": floor_near,
            "Ceiling-near occupancy (%)": ceiling_near,

            "KDE peak": kde_peak,
            "Compression Index (CI)": ci,
            "Threshold Density Ratio (TDR) primary": tdr_primary,
            "Threshold Density Ratio (TDR) secondary": tdr_secondary,

            "Variance homogeneity check": "Not applicable at indicator-only stage",
            "Manual note": spec["manual_note"],
        })

    out = pd.DataFrame(rows).sort_values("_order").drop(columns=["_order"]).reset_index(drop=True)
    return out

# ============================================================
# Table III: Manual-region structure
# ============================================================
def build_table3(df):
    rows = []

    for idx, spec in enumerate(indicator_specs):
        values = pd.to_numeric(df[spec["column"]], errors="coerce").dropna().to_numpy(dtype=float)
        n = len(values)

        for region in spec["elevated_regions"]:
            region_info = region_structure(
                values,
                region["lower"],
                region["upper"],
                spec["scale_min"],
                spec["scale_max"],
            )

            region_start = region["lower"] if region["lower"] is not None else spec["scale_min"]
            region_end = region["upper"] if region["upper"] is not None else spec["scale_max"]

            pattern = region_pattern_label(
                region_info["Region n"],
                region_info["Tie proportion (%)"],
                region_info["Observed width"],
                region_start,
                region_end,
                region_info["Unique values"],
            )

            rows.append({
                "_order": idx,
                "Instrument": spec["instrument"],
                "Indicator": spec["indicator"],
                "Manual region": region["label"],
                "Manual region start": region_start,
                "Manual region end": region_end,

                "Region n": region_info["Region n"],
                "Region %": pct(region_info["Region n"], n),
                "Observed occupied start": region_info["Observed start"],
                "Observed occupied end": region_info["Observed end"],
                "Observed occupied width": region_info["Observed width"],

                "Regional unique values": region_info["Unique values"],
                "Regional tie proportion (%)": region_info["Tie proportion (%)"],
                "Regional IQR": region_info["IQR"],
                "Regional KDE peak": region_info["KDE peak"],
                "Region pattern": pattern,

                "Interpretation note": region_info["Interpretation note"],
                "Manual note": spec["manual_note"],
            })

    out = pd.DataFrame(rows).sort_values(["_order", "Manual region"]).drop(columns=["_order"]).reset_index(drop=True)
    return out


# ============================================================
# Table IV: RSST structural saturation classification
# ============================================================
def build_table4(df, table3):
    """
    Build an RSST structural saturation classification table.

    This table links the Stage 1 manual-region structure directly to the
    RSST formal architecture:

        S = D / G
        R = 1 / S

    where D is empirical density in the indicator-region and G is available
    regional granularity. The table also reports realized density per observed
    unique score value and observable compression flags.

    Important interpretation rule:
    S > 1 alone is not treated as sufficient proof of structural saturation.
    In bounded integer or low-granularity clinical metrics, repeated values can
    occur naturally. Therefore, structural saturation is classified as present
    only when density-granularity imbalance is accompanied by observable
    compression evidence.
    """
    rows = []

    for idx, spec in enumerate(indicator_specs):
        values = pd.to_numeric(df[spec["column"]], errors="coerce").dropna().to_numpy(dtype=float)
        user_delta = USER_DELTA_BY_COLUMN.get(spec["column"], DEFAULT_USER_DELTA)

        delta_min, delta_min_empirical, delta_min_warning = compute_delta_min(
            values,
            mode=DELTA_MIN_MODE,
            user_delta=user_delta,
        )

        k_prime = select_k_prime(
            spec,
            delta_mode=DELTA_MIN_MODE,
            refinement_mode="RSST",
        )

        k_prime_note = k_prime_rule_note(
            spec,
            delta_mode=DELTA_MIN_MODE,
        )

        for region in spec["elevated_regions"]:
            region_start = region["lower"] if region["lower"] is not None else spec["scale_min"]
            region_end = region["upper"] if region["upper"] is not None else spec["scale_max"]
            upper_exclusive = region["upper"] is not None

            # Reuse Table III information to keep the diagnostic pipeline aligned.
            t3_match = table3[
                (table3["Instrument"] == spec["instrument"]) &
                (table3["Indicator"] == spec["indicator"]) &
                (table3["Manual region"] == region["label"])
            ]
            if len(t3_match) != 1:
                raise ValueError(
                    f"Table III lookup failed for {spec['instrument']} | {spec['indicator']} | {region['label']}"
                )
            t3 = t3_match.iloc[0]

            d_region = int(t3["Region n"])
            obs_start = t3["Observed occupied start"]
            obs_end = t3["Observed occupied end"]
            obs_width = t3["Observed occupied width"]
            u_observed = int(t3["Regional unique values"])
            tie_pct = t3["Regional tie proportion (%)"]
            reg_iqr = t3["Regional IQR"]
            reg_kde_peak = t3["Regional KDE peak"]
            pattern = t3["Region pattern"]

            full_region_width = float(region_end) - float(region_start)
            width_ratio = safe_div(obs_width, full_region_width)
            iqr_ratio = safe_div(reg_iqr, full_region_width)

            g_manual = count_metric_positions(region_start, region_end, delta_min, upper_exclusive=upper_exclusive)

            if d_region > 0:
                g_occupied = count_metric_positions(obs_start, obs_end, delta_min, upper_exclusive=False)
            else:
                g_occupied = np.nan

            s_manual = safe_div(d_region, g_manual)
            s_occupied = safe_div(d_region, g_occupied)
            s_realized = safe_div(d_region, u_observed)

            r_manual = safe_div(1.0, s_manual)
            r_occupied = safe_div(1.0, s_occupied)
            r_realized = safe_div(1.0, s_realized)

            high_tie_flag = bool(pd.notna(tie_pct) and tie_pct >= 50.0)
            reduced_unique_flag = bool(d_region >= 5 and u_observed <= max(3, int(d_region / 3)))
            restricted_width_flag = bool(pd.notna(width_ratio) and d_region >= 10 and width_ratio <= 0.25)
            low_iqr_flag = bool(pd.notna(iqr_ratio) and d_region >= 5 and iqr_ratio <= 0.25)
            dense_pattern_flag = pattern in ["Densely repeated region", "Narrow but heavily occupied region"]

            compression_flags = {
                "High tie burden": high_tie_flag,
                "Reduced unique-score representation": reduced_unique_flag,
                "Restricted occupied width": restricted_width_flag,
                "Low regional IQR": low_iqr_flag,
                "Dense/restricted region pattern": dense_pattern_flag,
            }
            compression_count = int(sum(compression_flags.values()))
            compression_flag_text = "; ".join([k for k, v in compression_flags.items() if v])
            if compression_flag_text == "":
                compression_flag_text = "None"

            rsst_status = classify_rsst_structural_status(
                d_region,
                s_manual,
                s_occupied,
                compression_count,
            )

            rows.append({
                "_order": idx,
                "Instrument": spec["instrument"],
                "Indicator": spec["indicator"],
                "Manual region": region["label"],
                "Manual region start": region_start,
                "Manual region end": region_end,
                "Region upper boundary rule": "Upper exclusive" if upper_exclusive else "Upper inclusive to scale maximum",

                "D_region": d_region,
                "Delta_min_mode": DELTA_MIN_MODE,
                "Delta_min_user": user_delta,
                "Delta_min_empirical": delta_min_empirical,
                "Delta_min_final": delta_min,
                "Delta_min_warning": delta_min_warning,
                "K_prime_preservation_for_later_Omega": k_prime,
                "K_prime_rule_note": k_prime_note,
                "G_manual": g_manual,
                "G_occupied": g_occupied,
                "U_observed": u_observed,

                "S_manual = D/G_manual": s_manual,
                "S_occupied = D/G_occupied": s_occupied,
                "S_realized = D/U_observed": s_realized,
                "R_manual = 1/S_manual": r_manual,
                "R_occupied = 1/S_occupied": r_occupied,
                "R_realized = 1/S_realized": r_realized,

                "Observed occupied width": obs_width,
                "Manual width": full_region_width,
                "Occupied width ratio": width_ratio,
                "Regional IQR": reg_iqr,
                "Regional IQR ratio": iqr_ratio,
                "Regional tie proportion (%)": tie_pct,
                "Regional KDE peak": reg_kde_peak,
                "Region pattern": pattern,

                "High tie burden flag": high_tie_flag,
                "Reduced unique-score representation flag": reduced_unique_flag,
                "Restricted occupied width flag": restricted_width_flag,
                "Low regional IQR flag": low_iqr_flag,
                "Dense/restricted region pattern flag": dense_pattern_flag,
                "Compression evidence count": compression_count,
                "Compression evidence flags": compression_flag_text,

                "RSST structural saturation status": rsst_status,
                "Omega pathway status": omega_pathway_status(rsst_status),
                "Refinement eligibility status": "Not assessed at Stage 1",
                "RSST Stage 1 interpretation": (
                    "Structural saturation classification only; Ω/RSST_Score requires a separate item-level or proxy-indicator refinement-eligibility module."
                ),
                "Manual note": spec["manual_note"],
            })

    out = pd.DataFrame(rows).sort_values(["_order", "Manual region"]).drop(columns=["_order"]).reset_index(drop=True)
    return out


# ============================================================
# Table V: Regional structural heterogeneity profile
# ============================================================
# This module adds a heterogeneity-oriented interpretation layer without
# changing the existing Tables I-IV calculations or saturation decisions.
# It evaluates whether cutoff-defined/manual regions contain internally
# dispersed or asymmetric structures that may be clinically meaningful.
#
# Important distinction:
# - Table IV is the primary RSST structural saturation table.
#   It asks whether a region is structurally saturated/compressed.
# - Table V is a complementary regional heterogeneity layer.
#   It does not replace or modify Table IV. It asks whether the same
#   cutoff/reference region also shows internal dispersion, peak-median
#   displacement, lower/upper asymmetry, regional-normality mismatch,
#   or multimodal shape signals.
#
# A region can therefore be:
# - predominantly compressed/homogeneous,
# - structurally heterogeneous,
# - mixed (saturation + heterogeneity), or
# - indeterminate if too sparse.
#
# Reporting rule:
# Table IV remains intact and primary. Table V complements Table IV;
# it does not replace it.
# ============================================================

def _safe_float(x):
    try:
        if pd.isna(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def _position_ratio(value, start, end):
    """Return position of a value inside a region, where 0=start and 1=end."""
    value = _safe_float(value)
    start = _safe_float(start)
    end = _safe_float(end)
    if pd.isna(value) or pd.isna(start) or pd.isna(end) or end <= start:
        return np.nan
    return float((value - start) / (end - start))


def _region_values(values, lower, upper):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return arr[mask_region(arr, lower, upper)]


def _try_skew(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3 or np.allclose(values, values[0]):
        return np.nan
    try:
        return float(skew(values, bias=False))
    except Exception:
        return np.nan


def _try_kurtosis(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 4 or np.allclose(values, values[0]):
        return np.nan
    try:
        return float(kurtosis(values, fisher=True, bias=False))
    except Exception:
        return np.nan


def _try_shapiro_p(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3 or values.size > 5000 or np.allclose(values, values[0]):
        return np.nan
    try:
        return float(shapiro(values)[1])
    except Exception:
        return np.nan


def _try_dagostino_p(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    # Regional normality diagnostic rule:
    # SciPy's D'Agostino-Pearson normaltest can technically run for n >= 8,
    # but it emits a reliability warning when n < 20 because the kurtosis
    # component is not considered valid for smaller regional samples.
    # For the regional heterogeneity module, we therefore use the more
    # conservative rule: compute D'Agostino-Pearson only when n >= 20.
    # Shapiro-Wilk remains available for 3 <= n <= 5000.
    if values.size < 20 or np.allclose(values, values[0]):
        return np.nan
    try:
        return float(normaltest(values)[1])
    except Exception:
        return np.nan


def _normality_category(flag):
    text = str(flag).lower()
    if "non-normal" in text:
        return "Non-normal"
    if "approximately normal" in text:
        return "Approximately normal"
    if "mixed" in text:
        return "Mixed"
    return "Not assessed"


def _kde_peak_count(values, grid_min, grid_max, min_relative_height=0.25):
    """
    Count exploratory KDE peaks inside a region.

    This is a descriptive multimodality signal, not a formal Hartigan dip test.
    Peaks below min_relative_height of the maximum density are ignored to avoid
    counting numerical ripples as clinically meaningful modes.
    """
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 5 or grid_max <= grid_min or np.allclose(vals, vals[0]):
        return np.nan, "Not assessed"

    try:
        x, y = compute_kde(vals, grid_min, grid_max, grid_points=800)
        if np.isnan(x).all() or np.nanmax(y) <= 0:
            return np.nan, "Not assessed"

        ymax = float(np.nanmax(y))
        peaks = []
        for i in range(1, len(y) - 1):
            if y[i] > y[i - 1] and y[i] >= y[i + 1] and y[i] >= min_relative_height * ymax:
                peaks.append(float(x[i]))

        # Merge very close peaks that can occur from numerical wiggles.
        merged = []
        min_sep = (grid_max - grid_min) * 0.08
        for p in peaks:
            if not merged or abs(p - merged[-1]) >= min_sep:
                merged.append(p)

        count = len(merged)
        if count >= 2:
            return count, "Possible multimodality"
        if count == 1:
            return count, "Single dominant peak"
        return 0, "No clear peak"
    except Exception:
        return np.nan, "Not assessed"


def _heterogeneity_status(region_n, heterogeneity_count, saturation_status):
    if region_n == 0:
        return "Absent - empty manual/reference region"
    if region_n < 5:
        return "Indeterminate - sparse region"

    saturated = saturation_status == "Structurally saturated"
    compression = str(saturation_status).startswith("Compression present")

    if saturated and heterogeneity_count >= 3:
        return "Mixed structural pattern - saturation plus heterogeneity"
    if heterogeneity_count >= 3:
        return "Structurally heterogeneous"
    if heterogeneity_count == 2:
        return "Possible structural heterogeneity"
    if saturated:
        return "Structurally homogeneous/compressed dominant"
    if compression:
        return "Compression dominant; heterogeneity not fully supported"
    return "No structural heterogeneity signal"


def _heterogeneity_interpretation(status):
    if status == "Mixed structural pattern - saturation plus heterogeneity":
        return "The region shows both compression/saturation and meaningful internal dispersion; visible homogeneity should not be interpreted as true uniformity."
    if status == "Structurally heterogeneous":
        return "The region shows multiple dispersion/asymmetry signals compatible with internal structural heterogeneity."
    if status == "Possible structural heterogeneity":
        return "The region shows limited but notable dispersion/asymmetry signals; interpret as possible heterogeneity."
    if status == "Structurally homogeneous/compressed dominant":
        return "The region is dominated by compression/saturation indicators with insufficient heterogeneity evidence."
    if status == "Compression dominant; heterogeneity not fully supported":
        return "The region shows compression evidence, but heterogeneity evidence is insufficient for a heterogeneity classification."
    if "Sparse" in status or "sparse" in status:
        return "The region is too sparse for stable heterogeneity interpretation."
    if "empty" in status:
        return "No observations occupy this region."
    return "No meaningful regional heterogeneity signal was detected by the current descriptive criteria."


def build_table5(df, table2, table3, table4):
    """
    Build a regional structural heterogeneity profile table.

    The table keeps the existing Table IV saturation classifications intact and
    adds region-level heterogeneity diagnostics based on median/peak position,
    lower-vs-upper median split, regional normality, peak shape, and dispersion.
    """
    rows = []

    for idx, spec in enumerate(indicator_specs):
        values = pd.to_numeric(df[spec["column"]], errors="coerce").dropna().to_numpy(dtype=float)
        t2_match = table2[(table2["Instrument"] == spec["instrument"]) & (table2["Indicator"] == spec["indicator"])]
        whole_normality = t2_match.iloc[0]["Normality flag"] if len(t2_match) == 1 else "Not assessed"
        whole_norm_cat = _normality_category(whole_normality)

        for region in spec["elevated_regions"]:
            region_start = region["lower"] if region["lower"] is not None else spec["scale_min"]
            region_end = region["upper"] if region["upper"] is not None else spec["scale_max"]
            region_width = float(region_end) - float(region_start)

            vals = _region_values(values, region["lower"], region["upper"])
            n = int(vals.size)

            t3_match = table3[
                (table3["Instrument"] == spec["instrument"]) &
                (table3["Indicator"] == spec["indicator"]) &
                (table3["Manual region"] == region["label"])
            ]
            t4_match = table4[
                (table4["Instrument"] == spec["instrument"]) &
                (table4["Indicator"] == spec["indicator"]) &
                (table4["Manual region"] == region["label"])
            ]
            if len(t3_match) != 1 or len(t4_match) != 1:
                raise ValueError(f"Table V lookup failed for {spec['instrument']} | {spec['indicator']} | {region['label']}")

            t3 = t3_match.iloc[0]
            t4 = t4_match.iloc[0]

            obs_start = _safe_float(t3["Observed occupied start"])
            obs_end = _safe_float(t3["Observed occupied end"])
            obs_width = _safe_float(t3["Observed occupied width"])
            u_observed = int(t3["Regional unique values"])
            tie_pct = _safe_float(t3["Regional tie proportion (%)"])
            reg_iqr = _safe_float(t3["Regional IQR"])
            reg_peak = _safe_float(t3["Regional KDE peak"])
            saturation_status = str(t4["RSST structural saturation status"])

            if n > 0:
                reg_median = float(np.median(vals))
                lower_half = vals[vals <= reg_median]
                upper_half = vals[vals > reg_median]
                lower_n = int(lower_half.size)
                upper_n = int(upper_half.size)
                lower_unique = int(pd.Series(lower_half).nunique()) if lower_n > 0 else 0
                upper_unique = int(pd.Series(upper_half).nunique()) if upper_n > 0 else 0
                lower_tie = tie_proportion(lower_half) if lower_n > 1 else 0.0
                upper_tie = tie_proportion(upper_half) if upper_n > 1 else 0.0
            else:
                reg_median = np.nan
                lower_n = upper_n = lower_unique = upper_unique = 0
                lower_tie = upper_tie = np.nan

            median_position = _position_ratio(reg_median, region_start, region_end)
            peak_position = _position_ratio(reg_peak, region_start, region_end)
            peak_to_cutoff = safe_div(reg_peak - region_start, region_width) if pd.notna(reg_peak) else np.nan
            median_to_cutoff = safe_div(reg_median - region_start, region_width) if pd.notna(reg_median) else np.nan
            peak_to_median_abs = abs(reg_peak - reg_median) if pd.notna(reg_peak) and pd.notna(reg_median) else np.nan
            peak_to_median_ratio = safe_div(peak_to_median_abs, region_width)
            occupied_width_ratio = safe_div(obs_width, region_width)
            regional_iqr_ratio = safe_div(reg_iqr, region_width)
            unique_value_ratio = safe_div(u_observed, n)
            lower_upper_n_ratio = safe_div(lower_n, upper_n)
            lower_upper_unique_ratio = safe_div(lower_unique, upper_unique)
            lower_upper_tie_difference = lower_tie - upper_tie if pd.notna(lower_tie) and pd.notna(upper_tie) else np.nan

            reg_skew = _try_skew(vals)
            reg_kurt = _try_kurtosis(vals)
            reg_shapiro_p = _try_shapiro_p(vals)
            reg_dagostino_p = _try_dagostino_p(vals)
            reg_normality = normality_flag(reg_shapiro_p, reg_dagostino_p)
            reg_norm_cat = _normality_category(reg_normality)
            norm_mismatch = (
                reg_norm_cat != "Not assessed" and
                whole_norm_cat != "Not assessed" and
                reg_norm_cat != whole_norm_cat
            )

            peak_count, multimodality_signal = _kde_peak_count(vals, region_start, region_end)

            # Heterogeneity evidence flags. These criteria are descriptive and
            # intentionally require convergence; no single flag is sufficient.
            broad_width_flag = bool(pd.notna(occupied_width_ratio) and n >= 5 and occupied_width_ratio >= 0.50)
            high_iqr_flag = bool(pd.notna(regional_iqr_ratio) and n >= 5 and regional_iqr_ratio >= 0.25)
            preserved_unique_flag = bool(pd.notna(unique_value_ratio) and n >= 5 and unique_value_ratio >= 0.40 and u_observed >= 5)
            median_displaced_flag = bool(pd.notna(median_position) and n >= 5 and median_position >= 0.25)
            peak_median_separation_flag = bool(pd.notna(peak_to_median_ratio) and n >= 5 and peak_to_median_ratio >= 0.15)
            lower_upper_asymmetry_flag = bool(
                n >= 5 and pd.notna(lower_upper_n_ratio) and
                (lower_upper_n_ratio <= 0.67 or lower_upper_n_ratio >= 1.50)
            )
            normality_mismatch_flag = bool(norm_mismatch)
            multimodality_flag = bool(pd.notna(peak_count) and peak_count >= 2)

            # Reporting label is source-aware only for wording.
            # This does not change the median-displacement calculation,
            # heterogeneity evidence count, or structural heterogeneity status.
            boundary_source = str(spec.get("boundary_source", "Cutoff"))
            if boundary_source == "Threshold/Reference":
                median_displacement_label = "Median displaced from Threshold/Reference"
            else:
                median_displacement_label = "Median displaced from cutoff"

            heterogeneity_flags = {
                "Broad occupied width": broad_width_flag,
                "High regional IQR ratio": high_iqr_flag,
                "Preserved unique-score representation": preserved_unique_flag,
                median_displacement_label: median_displaced_flag,
                "Peak-median separation": peak_median_separation_flag,
                "Lower/upper median-split asymmetry": lower_upper_asymmetry_flag,
                "Regional-vs-indicator normality mismatch": normality_mismatch_flag,
                "Possible multimodality": multimodality_flag,
            }
            heterogeneity_count = int(sum(heterogeneity_flags.values()))
            heterogeneity_flag_text = "; ".join([k for k, v in heterogeneity_flags.items() if v]) or "None"
            structural_heterogeneity_status = _heterogeneity_status(n, heterogeneity_count, saturation_status)

            rows.append({
                "_order": idx,
                "Instrument": spec["instrument"],
                "Indicator": spec["indicator"],
                "Manual/reference region": region["label"],
                "Region start/cutoff boundary": region_start,
                "Region end": region_end,
                "Region n": n,
                "Observed occupied start": obs_start,
                "Observed occupied end": obs_end,
                "Observed occupied width": obs_width,
                "Occupied width ratio": occupied_width_ratio,
                "Regional median": reg_median,
                "Regional KDE peak": reg_peak,
                "Median position ratio": median_position,
                "Peak position ratio": peak_position,
                "Median-to-cutoff distance ratio": median_to_cutoff,
                "Peak-to-cutoff distance ratio": peak_to_cutoff,
                "Peak-to-median absolute distance": peak_to_median_abs,
                "Peak-to-median distance ratio": peak_to_median_ratio,
                "Regional unique values": u_observed,
                "Unique value ratio": unique_value_ratio,
                "Regional tie proportion (%)": tie_pct,
                "Regional IQR": reg_iqr,
                "Regional IQR ratio": regional_iqr_ratio,
                "Lower-half n (<= regional median)": lower_n,
                "Upper-half n (> regional median)": upper_n,
                "Lower-half unique values": lower_unique,
                "Upper-half unique values": upper_unique,
                "Lower-half tie proportion (%)": lower_tie,
                "Upper-half tie proportion (%)": upper_tie,
                "Lower/upper n ratio": lower_upper_n_ratio,
                "Lower/upper unique ratio": lower_upper_unique_ratio,
                "Lower-upper tie difference": lower_upper_tie_difference,
                "Regional skewness": reg_skew,
                "Regional kurtosis": reg_kurt,
                "Regional Shapiro-Wilk p": reg_shapiro_p,
                "Regional D'Agostino-Pearson p": reg_dagostino_p,
                "Regional normality flag": reg_normality,
                "Whole-indicator normality flag": whole_normality,
                "Regional-vs-indicator normality mismatch": normality_mismatch_flag,
                "KDE peak count": peak_count,
                "KDE multimodality signal": multimodality_signal,
                "Table IV saturation status": saturation_status,
                "Heterogeneity evidence count": heterogeneity_count,
                "Heterogeneity evidence flags present": heterogeneity_flag_text,
                "Structural heterogeneity status": structural_heterogeneity_status,
                "Heterogeneity interpretation": _heterogeneity_interpretation(structural_heterogeneity_status),
            })

    out = pd.DataFrame(rows).sort_values(["_order", "Manual/reference region"]).drop(columns=["_order"]).reset_index(drop=True)
    return out


def build_table5_display(table5):
    cols = [
        "Instrument",
        "Indicator",
        "Manual/reference region",
        "Region n",
        "Regional median",
        "Regional KDE peak",
        "Median position ratio",
        "Peak position ratio",
        "Occupied width ratio",
        "Regional unique values",
        "Unique value ratio",
        "Regional tie proportion (%)",
        "Regional IQR ratio",
        "Regional normality flag",
        "Whole-indicator normality flag",
        "Regional-vs-indicator normality mismatch",
        "KDE peak count",
        "KDE multimodality signal",
        "Table IV saturation status",
        "Heterogeneity evidence count",
        "Heterogeneity evidence flags present",
        "Structural heterogeneity status",
        "Heterogeneity interpretation",
    ]
    return table5[cols].copy()


def plot_figure4(table5):
    if table5.empty:
        return
    fig, ax = plt.subplots(figsize=(11, 8), dpi=300)

    for _, r in table5.iterrows():
        x = 0.0 if pd.isna(r["Occupied width ratio"]) else float(r["Occupied width ratio"])
        y = 0.0 if pd.isna(r["Unique value ratio"]) else float(r["Unique value ratio"])
        size_base = 50 if pd.isna(r["Heterogeneity evidence count"]) else 50 + float(r["Heterogeneity evidence count"]) * 35
        color = get_instrument_color(r["Instrument"])
        ax.scatter(x, y, s=size_base, color=color, alpha=0.75, edgecolor="black", linewidth=0.8)
        ax.annotate(
            f"{short_label_for(r['Indicator'])} | {r['Manual/reference region']}",
            (x, y),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
        )

    ax.set_xlabel("Occupied width ratio within region")
    ax.set_ylabel("Unique value ratio within region")
    ax.set_title("Figure 4. Regional Structural Heterogeneity Map")
    ax.grid(alpha=0.2)
    plt.tight_layout()
    fig.savefig(FIG4_PNG, bbox_inches="tight")
    fig.savefig(FIG4_PDF, bbox_inches="tight")
    plt.close(fig)


def write_table5_note_and_summary(table5):
    TABLE5_NOTE_TXT.write_text(
        textwrap.dedent(
            """
            Table V note.
            Regional structural heterogeneity is evaluated inside the same manual/reference regions used in Tables III and IV. Table IV remains the primary RSST structural saturation classification. Table V does not replace Table IV and does not change any Table IV saturation classifications. Instead, Table V complements Table IV by examining regional median position, KDE peak position, peak-to-median separation, occupied width, unique-score representation, median-split lower/upper asymmetry, regional normality, regional-versus-whole-indicator normality mismatch, and exploratory KDE peak count.

            Structural heterogeneity is not inferred from any single statistic. It is classified only when multiple dispersion/asymmetry signals converge. A region may be structurally saturated and structurally heterogeneous at the same time; such cases are labeled as mixed structural patterns. This means that apparent homogeneity from score compression may coexist with clinically meaningful internal dispersion.

            The “Heterogeneity evidence flags present” column lists only positive heterogeneity evidence flags; criteria not listed were not met for that region.

            Regional normality diagnostic rule: Regional D'Agostino-Pearson normality testing is computed only for regions with n >= 20. For regions with n < 20, the D'Agostino-Pearson p-value is reported as missing/not assessed to avoid the SciPy kurtosis-test warning and to preserve conservative regional interpretation. Shapiro-Wilk remains the primary regional normality check when 3 <= n <= 5000.
            """
        ).strip(),
        encoding="utf-8",
    )

    lines = ["Regional heterogeneity summary", "==============================", ""]
    counts = table5["Structural heterogeneity status"].value_counts().to_dict()
    lines.append("Structural heterogeneity status counts:")
    for status, count in counts.items():
        lines.append(f"- {status}: {int(count)} region(s)")
    lines.append("")

    mixed = table5[table5["Structural heterogeneity status"] == "Mixed structural pattern - saturation plus heterogeneity"]
    heter = table5[table5["Structural heterogeneity status"] == "Structurally heterogeneous"]
    possible = table5[table5["Structural heterogeneity status"] == "Possible structural heterogeneity"]

    def add_section(title, sub):
        lines.append(title)
        if len(sub) == 0:
            lines.append("- None")
        else:
            for _, r in sub.iterrows():
                lines.append(
                    f"- {r['Instrument']} | {r['Indicator']} | {r['Manual/reference region']}: "
                    f"evidence={int(r['Heterogeneity evidence count'])}; flags={r['Heterogeneity evidence flags present']}"
                )
        lines.append("")

    add_section("Mixed saturation + heterogeneity regions:", mixed)
    add_section("Structurally heterogeneous regions:", heter)
    add_section("Possible structural heterogeneity regions:", possible)

    HETEROGENEITY_SUMMARY_TXT.write_text("\n".join(lines), encoding="utf-8")

# ============================================================
# Display versions
# ============================================================
def build_display_versions(table1, table2, table3, table4):
    # Table I is already dynamic; all generated Band k columns are retained.
    table1_display = table1.copy()

    table2_display = table2[
        [
            "Instrument", "Indicator",
            "N non-missing", "% missing",
            "Mean", "Median", "SD", "IQR",
            "Observed minimum", "Observed maximum",
            "Primary cutoff", "Secondary cutoff",
            "% >= primary cutoff", "% in secondary manual region",
            "Shapiro-Wilk p", "D'Agostino-Pearson p", "Normality flag",
            "Tie proportion (%)",
            "KDE peak", "Compression Index (CI)",
            "Threshold Density Ratio (TDR) primary",
            "Threshold Density Ratio (TDR) secondary",
            "Variance homogeneity check",
        ]
    ].copy()

    table3_display = table3[
        [
            "Instrument",
            "Indicator",
            "Manual region",
            "Region n",
            "Region %",
            "Observed occupied start",
            "Observed occupied end",
            "Observed occupied width",
            "Regional unique values",
            "Regional tie proportion (%)",
            "Regional IQR",
            "Regional KDE peak",
            "Region pattern",
            "Interpretation note",
        ]
    ].copy()


    table4_display = table4[
        [
            "Instrument",
            "Indicator",
            "Manual region",
            "D_region",
            "Delta_min_mode",
            "Delta_min_empirical",
            "Delta_min_final",
            "K_prime_preservation_for_later_Omega",
            "G_manual",
            "G_occupied",
            "U_observed",
            "S_manual = D/G_manual",
            "S_occupied = D/G_occupied",
            "S_realized = D/U_observed",
            "R_occupied = 1/S_occupied",
            "Compression evidence count",
            "Compression evidence flags",
            "RSST structural saturation status",
            "Omega pathway status",
            "Refinement eligibility status",
        ]
    ].copy()

    return table1_display, table2_display, table3_display, table4_display


# ============================================================
# Figures
# ============================================================
def _extract_pct_from_string(txt):
    if isinstance(txt, str) and "(" in txt and "%)" in txt:
        try:
            return float(txt.split("(")[1].replace("%)", ""))
        except Exception:
            return 0.0
    return 0.0


def plot_figure1(table1):
    instruments = table1["Instrument"].dropna().astype(str).drop_duplicates().tolist()
    if len(instruments) == 0:
        return

    fig_height = max(4, 3.4 * len(instruments))
    fig, axes = plt.subplots(nrows=len(instruments), ncols=1, figsize=(12, fig_height), dpi=300)
    if len(instruments) == 1:
        axes = [axes]

    band_label_cols = [c for c in table1.columns if re.match(r"Band \d+ label", str(c))]

    for ax, inst in zip(axes, instruments):
        sub = table1[table1["Instrument"] == inst].reset_index(drop=True)
        y = np.arange(len(sub))
        labels = sub["Indicator"].tolist()
        left = np.zeros(len(sub))

        band_order = []
        for _, r in sub.iterrows():
            for lab_col in band_label_cols:
                lab = r.get(lab_col, "")
                if isinstance(lab, str) and lab and lab not in band_order:
                    band_order.append(lab)

        for band in band_order:
            vals = []
            for _, r in sub.iterrows():
                val = 0.0
                for lab_col in band_label_cols:
                    if r.get(lab_col, "") == band:
                        n_col = lab_col.replace(" label", " n (%)")
                        val = _extract_pct_from_string(r.get(n_col, ""))
                        break
                vals.append(val)

            ax.barh(
                y,
                vals,
                left=left,
                label=band,
                color=get_band_color(band),
                edgecolor="black",
                linewidth=0.5,
            )
            left += np.array(vals)

        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlim(0, 100)
        ax.set_xlabel("Percent of observations")
        ax.set_title(inst)
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(axis="x", alpha=0.2)

    fig.suptitle("Figure 1. Manual Threshold Occupancy by Instrument", y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(FIG1_PNG, bbox_inches="tight")
    fig.savefig(FIG1_PDF, bbox_inches="tight")
    plt.close(fig)


def plot_figure2(table2):
    metrics = pd.DataFrame({
        "Missing %": table2["% missing"].to_numpy(),
        "Non-normal flag": table2["Normality flag"].str.contains("Non-normal").astype(float).to_numpy() * 100.0,
        "Tie %": table2["Tie proportion (%)"].to_numpy(),
        "Outlier %": table2["IQR mild outliers %"].to_numpy(),
        ">= primary cutoff %": table2["% >= primary cutoff"].to_numpy(),
        "Secondary manual %": table2["% in secondary manual region"].fillna(0.0).to_numpy(),
        "CI x100": table2["Compression Index (CI)"].to_numpy() * 100.0,
        "TDR primary x100": table2["Threshold Density Ratio (TDR) primary"].to_numpy() * 100.0,
        "TDR secondary x100": table2["Threshold Density Ratio (TDR) secondary"].fillna(0.0).to_numpy() * 100.0,
    })

    metrics.index = [f"{r['Instrument']} | {r['Indicator']}" for _, r in table2.iterrows()]

    fig, ax = plt.subplots(figsize=(11, 8), dpi=300)
    im = ax.imshow(metrics.to_numpy(dtype=float), aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(metrics.shape[1]))
    ax.set_xticklabels(metrics.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(metrics.shape[0]))
    ax.set_yticklabels(metrics.index)
    ax.set_title("Figure 2. Data Health and Manual-Threshold Audit Heatmap")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Scaled descriptive magnitude")
    plt.tight_layout()
    fig.savefig(FIG2_PNG, bbox_inches="tight")
    fig.savefig(FIG2_PDF, bbox_inches="tight")
    plt.close(fig)


def plot_figure3(table3):
    fig, ax = plt.subplots(figsize=(11, 8), dpi=300)

    for _, r in table3.iterrows():
        x = 0.0 if pd.isna(r["Region %"]) else r["Region %"]
        y = 0.0 if pd.isna(r["Regional unique values"]) else r["Regional unique values"]
        size = 50 if pd.isna(r["Regional IQR"]) else max(50, r["Regional IQR"] * 60)
        color = get_instrument_color(r["Instrument"])

        ax.scatter(x, y, s=size, color=color, alpha=0.75, edgecolor="black", linewidth=0.8)
        ax.annotate(
            f"{short_label_for(r['Indicator'])} | {r['Manual region']}",
            (x, y),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
        )

    ax.set_xlabel("Manual-region occupancy (%)")
    ax.set_ylabel("Regional unique values")
    ax.set_title("Figure 3. Manual Region Occupancy and Structure")
    ax.grid(alpha=0.2)
    plt.tight_layout()
    fig.savefig(FIG3_PNG, bbox_inches="tight")
    fig.savefig(FIG3_PDF, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Text outputs
# ============================================================
def recommendation_lines(table1, table2, table3, table4):
    lines = ["Indicator recommendations", "========================", ""]

    for _, r in table2.iterrows():
        inst = r["Instrument"]
        ind = r["Indicator"]

        lines.append(f"{inst} | {ind}")
        status = table1.loc[
            (table1["Instrument"] == inst) & (table1["Indicator"] == ind),
            "Cutoff/Threshold status"
        ].iloc[0]
        lines.append(f"- Manual status: {status}")

        if "Non-normal" in r["Normality flag"]:
            lines.append("- Distribution is descriptively non-normal; clinically meaningful observations should still be retained.")
        else:
            lines.append("- Normality was not significantly rejected at the available test level.")

        if r["Tie proportion (%)"] >= 25:
            lines.append("- Tie burden is substantial and may indicate weakening visible score resolution.")

        if pd.notna(r["Compression Index (CI)"]) and r["Compression Index (CI)"] >= 0.08:
            lines.append("- Compression Index is relatively elevated; the arithmetic center and the densest region of the distribution are noticeably separated.")

        if pd.notna(r["Threshold Density Ratio (TDR) primary"]) and r["Threshold Density Ratio (TDR) primary"] >= 0.50:
            lines.append("- A substantial share of estimated density lies in the official threshold-defined region.")

        rel = table3[(table3["Instrument"] == inst) & (table3["Indicator"] == ind)]
        if len(rel) > 0:
            if rel["Interpretation note"].str.contains("fewer than 5").any():
                lines.append("- At least one elevated manual region contains fewer than 5 observations; descriptive interpretation is limited in that region.")
            else:
                lines.append("- Elevated manual regions can be descriptively reported; manual thresholds remain primary.")

            if "Possible saturation signal" in " | ".join(rel["Interpretation note"].astype(str).tolist()):
                lines.append("- The elevated region shows a pattern potentially compatible with emerging saturation or apparent homogeneity.")

        rel4 = table4[(table4["Instrument"] == inst) & (table4["Indicator"] == ind)]
        if len(rel4) > 0:
            status_counts = rel4["RSST structural saturation status"].value_counts().to_dict()
            status_summary = "; ".join([f"{k}: {v}" for k, v in status_counts.items()])
            lines.append(f"- RSST structural saturation status by manual region: {status_summary}.")
            if (rel4["RSST structural saturation status"] == "Structurally saturated").any():
                lines.append("- At least one manual region is structurally saturated under RSST Stage 1 criteria; auxiliary-structure assessment is required before any Ω-based refinement.")
            elif rel4["RSST structural saturation status"].str.contains("Indeterminate", regex=False).any():
                lines.append("- At least one manual region is too sparse for stable RSST structural classification.")

        lines.append("- Manual thresholds remain primary; CI, TDR, KDE, tie burden, and regional structure summaries are descriptive only.")
        lines.append("- This first-stage output does not solve saturation directly; it identifies where a later saturation-focused step may be needed.")
        lines.append("")

    return "\n".join(lines)


def write_user_guide():
    guide_text = """
RSST_Stage1_Diagnostics.py
User Guide and Variable Directory
=================================

This program is a manual-threshold-first, Stage 1 structural diagnostic report.

Its role is to:
- read Cutoff or Threshold/Reference values from the workbook's Cutoff sheet,
- preserve the input-defined Cutoff or Threshold/Reference interpretation,
- describe the structure of the observed data,
- identify clustering, repetition, compression, threshold concentration, and possible Metric Silence risk,
- and help detect where visible score uniformity may be masking retained heterogeneity.

This Stage 1 module is designed to take an initial structural "X-ray" of the data. It evaluates whether bounded observed-score regions show evidence of local compression, saturation, tie burden, threshold concentration, or loss of local ordinal distinguishability. It does not authorize Ω-based refinement, does not apply Ω, and does not derive or report RSST_Score.

This version is dynamic and scale-agnostic. Any dataset with the same workbook
structure can be analyzed if:
- the Data sheet contains the score columns,
- the Cutoff sheet contains matching score-column headers,
- each score column has at least one Cutoff or Threshold/Reference value,
- and each score column has a Min/Max scale definition.

Supported Cutoff or Threshold/Reference formats:
- one value: 65
- two values: 60/70
- more than two values: 10/20/30 or separate Cutoff 1, Cutoff 2, Cutoff 3 rows

For k supplied values, the program generates k+1 manual bands:
<c1, c1 to <c2, ..., >=ck.

Cutoff values remain primary when supplied; Threshold/Reference values are used when Cutoff is empty.
Observed data remain intact.
All additional metrics are descriptive only.

Key metrics:
- Regional structural heterogeneity profile: median/peak position, regional normality, lower/upper median-split structure, and mixed saturation-plus-heterogeneity classification
- Compression Index (CI): descriptive separation between the mean and the KDE peak
- Threshold Density Ratio (TDR): descriptive density concentration in official threshold-defined regions
- KDE Peak: most densely occupied score location
- Tie Proportion: extent of repeated score values
- Region pattern: structural description of elevated manual regions
- D_region: empirical density inside the indicator-region
- G_manual / G_occupied: available regional granularity under the observed metric step
- S = D/G: RSST saturation index
- R = 1/S: RSST resolvable distinction
- K′ preservation setting: reported only as a recommended preservation setting for possible later Ω-based refinement; not used in Stage 1 structural saturation classification
- RSST structural saturation status: indicator-region classification based on density-granularity imbalance plus observable compression evidence

Omega (Ω) is not applied in this Stage 1 report. RSST_Score is not produced here.
Structural saturation or compression does not by itself justify Ω-based refinement.
Ω-based refinement requires a separate auxiliary-structure and refinement-eligibility assessment using item-level, component-level, or eligible proxy-indicator information.

K′ is not used in Stage 1 structural saturation classification. It is reported only as a recommended preservation setting for a possible later Ω-based refinement module. The recommended K′ should be interpreted according to the native observed-score resolution, not merely according to whether a cutoff/reference boundary is expressed as an integer or decimal value.

A visible score region can appear homogeneous without being truly homogeneous.
This program helps identify where such apparent homogeneity may be occurring.
"""
    USER_GUIDE_TXT.write_text(textwrap.dedent(guide_text).strip(), encoding="utf-8")


def write_notes_and_summary(table1, table2, table3, table4):
    PIPELINE_NOTE_TXT.write_text(
        textwrap.dedent(
            """
            Pipeline note.
            This program is a manual-threshold-first, Stage 1 structural diagnostic report. It does not create new cutoffs and does not replace official manual interpretation. Manual bands are generated only from the cutoff values supplied in the workbook's Cutoff sheet.

            The purpose of this program is twofold. First, it determines where the observed data fall relative to the supplied Cutoff or Threshold/Reference boundaries. Second, it identifies whether the observed score structure shows clustering, repetition, compression, or apparent homogeneity that may be compatible with emerging regional score saturation.

            This dynamic version supports one, two, three, four, or more supplied cutoffs for each indicator. For k supplied values, the program creates k+1 manual bands: <c1, c1 to <c2, ..., >=ck. Elevated manual regions are all regions at or above the first supplied cutoff.

            This program includes an RSST structural saturation classification layer. For each indicator-region, the program computes empirical density (D), available regional granularity (G), saturation indices (S = D/G), resolvable distinction (R = 1/S), and observable compression markers. Structural saturation is classified only when density-granularity imbalance is accompanied by observable compression evidence. Table IV is the primary RSST structural saturation table and remains intact.

            This program does not apply Omega (Ω), does not derive or report RSST_Score, and does not determine refinement eligibility. Refinement eligibility requires a separate auxiliary-structure and refinement-eligibility assessment using item-level or eligible proxy-indicator information.

            K′ is reported only as the preservation parameter for the later Omega / RSST_Score refinement module. It is not used to classify structural saturation in Stage 1. Under Δmin AUTO mode, integer-valued clinical score metrics use K′ = 10000, whereas decimal/fixed-precision metrics use K′ = 1000000 to maintain scale-faithful preservation under finer empirical step sizes.

            All additional quantities, including KDE peak, Compression Index (CI), Threshold Density Ratio (TDR), tie burden, outlier flags, regional structure summaries, RSST structural saturation outputs, and regional heterogeneity outputs, are descriptive/diagnostic only and do not replace manual interpretation. Table V complements Table IV by evaluating whether the same cutoff/reference regions also contain regional structural heterogeneity.
            """
        ).strip(),
        encoding="utf-8",
    )

    TABLE1_NOTE_TXT.write_text(
        textwrap.dedent(
            """
            Table I note.
            Threshold-based status is defined exclusively by the Cutoff or Threshold/Reference values supplied in the Cutoff sheet. No percentile anchors, derived sub-bands, or additional cutoffs are created at this stage. This table describes only where the observed cohort falls relative to the supplied cutoff/reference boundaries.
            """
        ).strip(),
        encoding="utf-8",
    )

    TABLE2_NOTE_TXT.write_text(
        textwrap.dedent(
            """
            Table II note.
            Data-health and assumption diagnostics are descriptive only. Clinically meaningful observations are retained. Outliers are flagged but not removed, winsorized, or normalized. Normality was assessed with Shapiro-Wilk and, when sample size permitted, the D'Agostino-Pearson omnibus normality test. Homogeneity of variance is not assessed at this indicator-only stage because it requires a grouping variable and belongs to the later inferential layer.

            This table should also be read as a structural warning layer. High tie burden, strong threshold concentration, elevated Compression Index (CI), or elevated Threshold Density Ratio (TDR) may indicate that visible score differentiation is weakening within bounded score regions. These quantities do not establish saturation by themselves, but they help identify where a later saturation-focused analytic step may be needed.
            """
        ).strip(),
        encoding="utf-8",
    )

    TABLE3_NOTE_TXT.write_text(
        textwrap.dedent(
            """
            Table III note.
            Manual-region structure is reported only inside supplied manual elevated regions. "Observed occupied start" and "Observed occupied end" describe where data actually occur inside the manual region. Region n, regional unique values, regional tie proportion, and regional IQR are reported directly. If a manual region contains fewer than 5 observations, the region is still reported, but descriptive interpretation is explicitly limited.

            The purpose of this table is to determine whether elevated official regions appear sparse, broadly occupied, densely repeated, or potentially compatible with emerging regional score saturation. Visible uniformity within an elevated region should not automatically be interpreted as true homogeneity.
            """
        ).strip(),
        encoding="utf-8",
    )


    TABLE4_NOTE_TXT.write_text(
        textwrap.dedent(
            """
            Table IV note.
            RSST structural saturation classification is assigned at the indicator-region level. D_region denotes empirical density inside the manual region. G_manual denotes available metric granularity across the full official manual region, whereas G_occupied denotes available metric granularity across the actually occupied observed span. S = D/G indexes density-granularity imbalance, and R = 1/S indexes resolvable distinction.

            S > 1 alone is not treated as sufficient evidence of structural saturation. Structural saturation is classified as present only when density-granularity imbalance is accompanied by observable compression evidence, including high tie burden, reduced unique-score representation, restricted occupied width, low regional IQR, or a dense/restricted region pattern.

            K′ is reported only as the preservation parameter for the later Ω / RSST_Score refinement module. It is not used to classify structural saturation in Stage 1. Under Δmin AUTO mode, integer-valued clinical score metrics use K′ = 10000, whereas decimal/fixed-precision metrics use K′ = 1000000 to maintain scale-faithful preservation under finer empirical step sizes.

            Table IV does not apply Omega (Ω), does not derive or report RSST_Score, and does not determine refinement eligibility. It identifies whether a structurally saturated indicator-region should proceed to a separate auxiliary-structure and refinement-eligibility assessment. Table V complements Table IV by evaluating regional structural heterogeneity in the same cutoff/reference regions; Table V does not replace, modify, or reduce the role of Table IV.
            """
        ).strip(),
        encoding="utf-8",
    )

    summary_lines = []
    summary_lines.append("Executive summary")
    summary_lines.append("=================")
    summary_lines.append("")
    summary_lines.append(f"Indicators analyzed: {len(table1)}")
    summary_lines.append("")
    summary_lines.append("Cutoff/Threshold status summary:")

    for _, r in table1.iterrows():
        summary_lines.append(f"- {r['Instrument']} | {r['Indicator']}: {r['Cutoff/Threshold status']}")

    summary_lines.append("")
    summary_lines.append("Top Compression Index indicators:")
    top_ci = table2.sort_values("Compression Index (CI)", ascending=False).head(5)
    for _, r in top_ci.iterrows():
        summary_lines.append(f"- {r['Instrument']} | {r['Indicator']}: CI={r['Compression Index (CI)']:.3f}")

    summary_lines.append("")
    summary_lines.append("RSST structural saturation status summary:")
    status_counts = table4["RSST structural saturation status"].value_counts()
    for status, count in status_counts.items():
        summary_lines.append(f"- {status}: {int(count)} indicator-region(s)")

    saturated = table4[table4["RSST structural saturation status"] == "Structurally saturated"]
    summary_lines.append("")
    summary_lines.append("Structurally saturated indicator-regions:")
    if len(saturated) == 0:
        summary_lines.append("- None")
    else:
        for _, r in saturated.iterrows():
            summary_lines.append(
                f"- {r['Instrument']} | {r['Indicator']} | {r['Manual region']}: "
                f"S_occupied={r['S_occupied = D/G_occupied']:.3f}, flags={r['Compression evidence flags']}"
            )

    summary_lines.append("")
    summary_lines.append("Manual regions with n < 5:")
    small = table3[table3["Region n"] < 5]
    if len(small) == 0:
        summary_lines.append("- None")
    else:
        for _, r in small.iterrows():
            summary_lines.append(f"- {r['Instrument']} | {r['Indicator']} | {r['Manual region']}: n={int(r['Region n'])}")

    EXEC_SUMMARY_TXT.write_text("\n".join(summary_lines), encoding="utf-8")
    RECOMMENDATIONS_TXT.write_text(recommendation_lines(table1, table2, table3, table4), encoding="utf-8")


# ============================================================
# Save tables
# ============================================================
def save_tables(table1, table1_display, table2, table2_display, table3, table3_display, table4, table4_display, table5=None, table5_display=None):
    table1.to_excel(TABLE1_FULL_XLSX, index=False)
    table1.to_csv(TABLE1_FULL_CSV, index=False)
    table1_display.to_excel(TABLE1_DISPLAY_XLSX, index=False)
    table1_display.to_csv(TABLE1_DISPLAY_CSV, index=False)

    table2.to_excel(TABLE2_FULL_XLSX, index=False)
    table2.to_csv(TABLE2_FULL_CSV, index=False)
    table2_display.to_excel(TABLE2_DISPLAY_XLSX, index=False)
    table2_display.to_csv(TABLE2_DISPLAY_CSV, index=False)

    table3.to_excel(TABLE3_FULL_XLSX, index=False)
    table3.to_csv(TABLE3_FULL_CSV, index=False)
    table3_display.to_excel(TABLE3_DISPLAY_XLSX, index=False)
    table3_display.to_csv(TABLE3_DISPLAY_CSV, index=False)

    table4.to_excel(TABLE4_FULL_XLSX, index=False)
    table4.to_csv(TABLE4_FULL_CSV, index=False)
    table4_display.to_excel(TABLE4_DISPLAY_XLSX, index=False)
    table4_display.to_csv(TABLE4_DISPLAY_CSV, index=False)

    if table5 is not None and table5_display is not None:
        table5.to_excel(TABLE5_FULL_XLSX, index=False)
        table5.to_csv(TABLE5_FULL_CSV, index=False)
        table5_display.to_excel(TABLE5_DISPLAY_XLSX, index=False)
        table5_display.to_csv(TABLE5_DISPLAY_CSV, index=False)

# ============================================================
# Main
# ============================================================
def main():
    print("Running with:", sys.executable)

    parser = argparse.ArgumentParser(
        description="Dynamic manual-threshold-first RSST Stage 1 report."
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default=str(INPUT_FILE),
        help="Excel workbook with Data and Cutoff sheets. Default: RSST-Data-1_Cutoff_Filled.xlsx",
    )
    parser.add_argument(
        "--data-sheet",
        default=DATA_SHEET_NAME,
        help="Name of the data sheet. Default: Data",
    )
    parser.add_argument(
        "--cutoff-sheet",
        default=CUTOFF_SHEET_NAME,
        help="Name of the cutoff sheet. Default: Cutoff",
    )
    args = parser.parse_args()

    input_file = Path(args.input_file)
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file.resolve()}")

    global indicator_specs
    df, indicator_specs = load_indicator_specs_from_workbook(
        input_file,
        data_sheet=args.data_sheet,
        cutoff_sheet=args.cutoff_sheet,
    )
    initialize_dynamic_plot_labels(indicator_specs)

    print(f"Loaded indicators from Cutoff sheet: {len(indicator_specs)}")
    for spec in indicator_specs:
        cutoffs = "; ".join(_format_threshold_value(v) for v in spec.get("cutoff_values", []))
        print(f"- {spec['instrument']} | {spec['indicator']} | {spec['column']} | cutoff/threshold values: {cutoffs}")

    table1 = build_table1(df)
    table2 = build_table2(df)
    table3 = build_table3(df)
    table4 = build_table4(df, table3)
    table5 = build_table5(df, table2, table3, table4)

    table1_display, table2_display, table3_display, table4_display = build_display_versions(table1, table2, table3, table4)
    table5_display = build_table5_display(table5)

    save_tables(table1, table1_display, table2, table2_display, table3, table3_display, table4, table4_display, table5, table5_display)
    plot_figure1(table1)
    plot_figure2(table2)
    plot_figure3(table3)
    plot_figure4(table5)
    write_notes_and_summary(table1, table2, table3, table4)
    write_table5_note_and_summary(table5)
    write_user_guide()

    print("\nSaved outputs to:", OUTPUT_DIR.resolve())
    print("Table I full:", TABLE1_FULL_XLSX)
    print("Table II full:", TABLE2_FULL_XLSX)
    print("Table III full:", TABLE3_FULL_XLSX)
    print("Table IV full:", TABLE4_FULL_XLSX)
    print("Table V full:", TABLE5_FULL_XLSX)
    print("Figure 1:", FIG1_PNG)
    print("Figure 2:", FIG2_PNG)
    print("Figure 3:", FIG3_PNG)
    print("Figure 4:", FIG4_PNG)
    print("Executive summary:", EXEC_SUMMARY_TXT)
    print("Recommendations:", RECOMMENDATIONS_TXT)
    print("Heterogeneity summary:", HETEROGENEITY_SUMMARY_TXT)
    print("User guide:", USER_GUIDE_TXT)

    print("\nTable I preview:")
    print(table1_display.to_string(index=False))
    print("*** End of program ***")

if __name__ == "__main__":
    main()
