# Synthetic N = 347 Item-Level Public Code-Execution Inputs

The synthetic data files included in this folder are item-level datasets prepared solely for public code execution, output-structure verification, and workflow inspection for the RSST/SFORR supplementary computational suite. These files are not the empirical N = 103 caregiver dataset used for the manuscript’s empirical analyses and should not be interpreted as participant-level empirical evidence. The synthetic datasets allow Supplementary Files S1-S3 to be executed and inspected without sharing the original study data.

Target scores are recomputed from synthetic item responses under the declared scoring rules. The files are intended to support reproducible code execution, input/output structure review, and verification of the public workflow. They are not intended to provide empirical findings, clinical evidence, or manuscript results.

## Governance constraints applied

* Synthetic sample size: N = 347 per indicator.
* Item-level structure is preserved for each indicator.
* Target scores are recomputed from synthetic item responses under the declared scoring rules.
* Each perturbed item differs from its source pattern by no more than 1 response category.
* Each synthetic record differs from the source-pattern item-sum total by no more than 2 points.
* The public files do not include source-row identifiers or participant-level records.
* The synthetic files are provided only for public code execution, output-structure verification, and workflow inspection.

## Included indicators

* BDI-II
* MMCGI_HSL
* MMCGI_PSB
* MMCGI_WFI
* PANAS_NA
* PANAS_PA_SeverityAligned
* ZBI

## Validation summary

See `Synthetic_N347_validation_summary.xlsx` for an audit of item ranges, target-score ranges, item-target consistency, and perturbation constraints.

Generated with deterministic synthetic-data seed: 20260602.

