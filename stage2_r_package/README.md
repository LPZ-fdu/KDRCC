# KDR-CC R package and independent validation

Kernel-Based Doubly Robust Case-Crossover Regression (KDR-CC), implemented in kdrcc 0.1.0.

Maintainer: Peizheng Li, pzli25@m.fudan.edu.cn. Package license: GPL-3. Requirements: 64-bit R >= 3.6.3 and CompQuadForm >= 1.4.4.

## Contents

| Directory | Contents |
|---|---|
| `package/kdrcc` | Installable source repository, R functions, help, package tests and CI specification |
| `package/installers` | Source archives and Windows R 4.4 binaries for kdrcc and its dependency |
| `package/kdrcc_user_guide_zh.docx` | Complete Chinese user guide |
| `package/kdrcc_user_guide_en.docx` | Complete English user guide |
| `package/examples` | Matched case data, complete analysis, configuration selection, custom graphics and worked outputs |
| `validation` | Independent numerical checks, frozen reference inputs, matched formal G1 reproduction and result reporting |
| `validation/results` | Completed 100-repetition Linux reproduction and its final figures, tables and numerical agreement |
| `validation/records` | Recorded installation workflow |

## Install and use

From this distribution directory, run in R:

```r
source("package/INSTALL_AND_TEST.R", encoding="UTF-8")
source("package/examples/RUN_EXAMPLE.R", encoding="UTF-8")
```

Read either Word user guide for data preparation, every function and argument, auxiliary specifications, prediction, curves, OR contrasts, variable and interaction tests, contribution shares, grouped tuning, diagnostics, saving models and custom ggplot2 output. The user guides describe ordinary package use independently of the simulation comparison.

Kernel learning records the candidate center count, retained basis dimension and feature screening decisions. In grouped tuning, eligibility requires finite held-out loss and an available coefficient covariance in every training fold. Final estimation computes covariance on the independent final sample. Package tests cover feature screening, candidate eligibility and prediction with a reduced basis.

## Reproduce the independent package assessment

```r
setwd("validation")
source("RUN_VALIDATION.R", encoding="UTF-8")
```

`RUN_VALIDATION.R` executes the numerical suite and matched formal G1 reproduction. `RUN_END_TO_END.R` runs the matched cohort directly. `PREPARE_RESULTS.R` redraws the final G1 figures from saved summaries. `END_TO_END_SETTINGS.R` sets repetition IDs, worker processes and execution options. The archived complete study uses the same 100 independent repetitions and fixed configuration as the first-stage G1 study: N=20,000, L=20,000, M=144, kappa=0.0001 and lambda=1. It records 2,900/2,900 numerical comparisons passing under Linux R 3.6.3.

Temporary per-repetition fits and logs support resumption during execution. Complete reporting retains aggregated results and clears those checkpoints. The original frozen reference records remain input data for reproducibility.

## Files and integrity

All runnable examples use paths relative to their own script. Source-package documentation lives in `package/kdrcc/man` and `package/kdrcc/inst/doc`. `CHECKSUMS.csv` records the distributed payload; `validation/VALIDATION_CHECKSUMS.csv` records the independent assessment materials. Source, example and validation code are supplied together with their reference data and completed formal outputs.
