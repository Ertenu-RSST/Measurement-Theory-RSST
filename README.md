# Regional Score Saturation Theory (RSST) / SFORR Supplementary Code Package

This repository contains the public synthetic N = 347 code-execution package for the Regional Score Saturation Theory (RSST) / Scale-Faithful Ordinal Resolution Refinement (SFORR) supplementary computational suite.

The package provides reproducible public materials for executing and inspecting the three supplementary computational layers associated with the RSST framework:

- Supplementary File S1: Stage-1 structural saturation diagnostics.
- Supplementary File S2: ITEM_LEVEL_GRM-based GRM_RSST_Score derivation and reporting.
- Supplementary File S3: GRM refit-bootstrap ordering-stability diagnostics.

The synthetic N = 347 files are provided only for public code execution, output-structure verification, and workflow inspection. They are not the empirical N = 103 caregiver dataset used for the manuscript’s empirical analyses and should not be interpreted as participant-level empirical evidence, clinical evidence, or manuscript results.

## Repository purpose

This repository is intended to allow reviewers, readers, and other users to inspect and execute the public RSST/SFORR supplementary programs without access to the original empirical caregiver dataset.

The repository supports:

- inspection of the RSST/SFORR computational architecture;
- execution of the S1, S2, and S3 scripts on synthetic public data;
- review of expected input and output workbook structures;
- verification of ITEM_LEVEL_GRM input governance;
- inspection of SILENT-mode and Omega-withholding logic where applicable;
- review of K′ preservation-bound implementation;
- review of residual auxiliary-unresolved tie audit outputs;
- inspection of GRM refit-bootstrap ordering-stability diagnostics.

The synthetic data and synthetic outputs are not used to generate the manuscript’s empirical tables or empirical findings.

## Supplementary computational layers

### Supplementary File S1: Stage-1 structural diagnostic pipeline

`RSST_Stage1_Diagnostics_v1.py`

Supplementary File S1 implements the Stage-1 structural diagnostic layer of RSST/SFORR. It evaluates bounded observed-score indicators for regional saturation, compression, tied-score burden, manual threshold-region occupancy, density–granularity imbalance, Metric Silence screening, and regional heterogeneity.

S1 is diagnostic only. It does not estimate a graded response model, does not call R `mirt`, does not construct a GRM-derived auxiliary ordering vector `A_aux`, does not apply Omega (Ω), and does not derive or report `GRM_RSST_Score`.

### Supplementary File S2: ITEM_LEVEL_GRM GRM_RSST_Score derivation and reporting pipeline

`RSST-Scale-Faithful_Ordinal_Resolution_Refinement_SFORR_v1.py`

Supplementary File S2 implements the production ITEM_LEVEL_GRM derivation-and-reporting pathway. When valid target-specific ordinal item responses are available and all eligibility conditions are satisfied, S2 fits a unidimensional graded response model using R `mirt`, obtains the GRM-derived auxiliary ordering vector `A_aux`, and applies the Omega-governed, K′-bounded, boundary-preserving refinement pathway to derive and report `GRM_RSST_Score`.

If required eligibility conditions are not satisfied, S2 withholds Omega, activates SILENT mode, preserves the native observed score, and does not introduce unsupported within-cell distinctions.

S2 does not break residual auxiliary-unresolved ties arbitrarily. If observations share the same native observed target score and the same GRM-derived `A_aux` value, they remain tied in `GRM_RSST_Score` and are documented in the audit outputs.

### Supplementary File S3: GRM refit-bootstrap ordering-stability diagnostic

`RSST_GRM_refit_bootstrap_ordering_stability_diagnostic_program_v1.py`

Supplementary File S3 implements a supplementary GRM refit-bootstrap ordering-stability diagnostic for the ITEM_LEVEL_GRM pathway. It evaluates the sampling/refit sensitivity of the GRM-derived auxiliary ordering vector `A_aux` inside native observed-score tied cells.

S3 does not derive, report, or modify `GRM_RSST_Score`. It does not apply Omega (Ω), does not select or modify K′, does not alter native observed scores, does not modify thresholds or severity bands, and does not override SILENT governance from S2.

For each bootstrap/refit iteration, S3 refits the GRM through R `mirt`, scores the original case set under the refit model, and compares within-cell pairwise ordering directions against the reference `A_aux` ordering. Failed, nonconverged, or non-finite refit iterations are excluded from the stability summary and documented in the output workbook.

## Synthetic N = 347 public data

The synthetic item-level datasets are prepared solely for public code execution, output-structure verification, and workflow inspection.

The synthetic datasets preserve item-level structure and recompute target scores from synthetic item responses under the declared scoring rules. They do not contain original participant-level records or source-row identifiers.

Included indicators:

- BDI-II
- MMCGI_HSL
- MMCGI_PSB
- MMCGI_WFI
- PANAS_NA
- PANAS_PA_SeverityAligned
- ZBI

Additional details are provided in `README_SYNTHETIC_N347.md` and `Synthetic_N347_validation_summary.xlsx`.

## Expected repository structure

The repository is organized around the three supplementary computational layers:

```text
RSST_SYN_N347/
│
├─ README.md
├─ README_SYNTHETIC_N347.md
├─ Synthetic_N347_validation_summary.xlsx
├─ requirements.txt
├─ setup_R_mirt_environment.R
├─ LICENSE
├─ CITATION.cff
├─ .zenodo.json
├─ CHECKSUMS_SHA256.txt
│
├─ S1/
│  ├─ RSST_Stage1_Diagnostics_v1.py
│  ├─ RSST-Data_Cutoff_Filled_SYN.xlsx
│  └─ RSST_Stage1_Diagnostics-v1/
│
├─ S2/
│  ├─ RSST-Scale-Faithful_Ordinal_Resolution_Refinement_SFORR_v1.py
│  ├─ GRM-Target_with_Ordinal_Items_<INDICATOR>_Input-SYN.xlsx
│  └─ <INDICATOR>_Score/
│
├─ S3-SYN/
│  ├─ RSST_GRM_refit_bootstrap_ordering_stability_diagnostic_program_v1.py
│  ├─ GRM-Target_with_Ordinal_Items_<INDICATOR>_Input.xlsx
│  └─ <INDICATOR>/
│
└─ synthetic_N347/
   └─ <INDICATOR>/
```

## Software environment

The package was prepared and tested using:

- Python 3.10.7
- R 4.5.3
- R package: `mirt`
- Python packages listed in `requirements.txt`

The S2 and S3 pipelines use Python for workflow governance and R `mirt` as the graded response model backend. Users may need to update the local `Rscript.exe` path in the S2 and S3 scripts before execution.

## Python setup

Install the required Python packages:

```bash
pip install -r requirements.txt
```

## R setup

The file `setup_R_mirt_environment.R` provides a helper script for installing or verifying the R `mirt` package.

Run in R or through Rscript:

```bash
Rscript setup_R_mirt_environment.R
```

This script is a setup and verification helper. It does not create a locked R package snapshot.

## Running Supplementary File S1

From the `S1/` folder, run:

```bash
python RSST_Stage1_Diagnostics_v1.py
```

The expected public input workbook is:

```text
RSST-Data_Cutoff_Filled_SYN.xlsx
```

S1 generates Stage-1 diagnostic tables, CSV outputs, figures, methodology notes, table notes, executive summaries, recommendations, and user-guide materials. No `GRM_RSST_Score` is produced by S1.

## Running Supplementary File S2

From the `S2/` folder, run:

```bash
python RSST-Scale-Faithful_Ordinal_Resolution_Refinement_SFORR_v1.py
```

The public synthetic S2 input workbooks use the naming convention:

```text
GRM-Target_with_Ordinal_Items_<INDICATOR>_Input-SYN.xlsx
```

S2 validates target-specific item-level evidence, fits the GRM through R `mirt` when eligible, applies the Omega-governed boundary-preserving refinement pathway, and reports `GRM_RSST_Score` only when all eligibility and preservation-governance conditions are satisfied.

The output workbook includes score outputs, item-target consistency diagnostics, GRM/SILENT governance summaries, K′ preservation-bound diagnostics, Omega application summaries, auxiliary-ordering diagnostics, residual auxiliary-unresolved tie audits, and GRM-reference comparison workspaces where applicable.

## Running Supplementary File S3

From the `S3-SYN/` folder, set the desired values in the user-configuration section:

```python
INPUT_INDICATOR = "<INDICATOR>"
BOOTSTRAP_ITERATIONS = 20
```

Then run:

```bash
python RSST_GRM_refit_bootstrap_ordering_stability_diagnostic_program_v1.py
```

The public synthetic S3 input workbooks use the naming convention:

```text
GRM-Target_with_Ordinal_Items_<INDICATOR>_Input.xlsx
```

The diagnostic-intensity outputs are organized using:

```text
B20
B50
B100
B500
B1000
```

S3 output workbooks are written into indicator-specific folders and are paired with run logs such as:

```text
S3_Run_B20.txt
S3_Run_B50.txt
S3_Run_B100.txt
S3_Run_B500.txt
S3_Run_B1000.txt
```

## Output interpretation notes

`GRM_RSST_Score` is a secondary, scale-faithful ordinal refinement of the native observed score. It is not a replacement latent trait estimate, not a new primary outcome scale, and not an interval-scale measurement.

Residual auxiliary-unresolved ties may remain after S2. This is expected when the selected auxiliary ordering vector does not distinguish observations that share the same native observed target score. Such ties are retained and documented rather than broken by input row order.

Excel display precision is used for audit readability. `GRM_RSST_Score` values and GRM-reference comparison columns may be displayed to 14 decimal places in audit workbooks to preserve visibility of small preservation-scale offsets. Analytic comparisons, ranking, and validation use the stored full-precision numeric values; `GRM_RSST_Score` values are not rounded before ranking or validation.

## Checksum verification

A checksum file is provided:

```text
CHECKSUMS_SHA256.txt
```

Users may verify file integrity by computing SHA-256 hashes locally and comparing them with the values in this file.

## License

This repository is released under the MIT License. See `LICENSE`.

## Citation

Citation metadata are provided in `CITATION.cff`.

After Zenodo archival, the DOI will be added here:

```text
DOI: [to be added after Zenodo release]
```

## Important data-use note

The synthetic N = 347 data in this repository are provided only for public code execution, output-structure verification, and workflow inspection. They are not the empirical N = 103 caregiver dataset used in the manuscript’s empirical analyses and must not be interpreted as participant-level empirical evidence, clinical evidence, or manuscript results.
