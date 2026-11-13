# Variable identification in matched case-crossover simulations

This Stage 1 supplement compares identification of active and null exposures by KDR-CC, spline conditional logistic regression (Spline-CLR), and linear conditional logistic regression (Linear-CLR). The comparison uses the frozen simulation mechanism, repetition-level exposure truth, and individual matched sets.

| Setting | Comparison |
|---|---|
| B | Eight active exposures under four nuisance-model scenarios, with 100 independent repetitions |
| C | Eight active exposures and 8, 16, 40, or 80 null exposures, with 100 independent repetitions per setting |
| D | Eight active exposures under low, moderate, and high overall correlation, with 100 independent repetitions per setting and the low-correlation reference drawn from B |

KDR-CC uses the coordinate-derivative test from the formal study. Linear-CLR tests each exposure coefficient with a two-sided Wald test. Spline-CLR jointly tests the coefficients of each exposure's natural-spline basis. The CLR Wald tests use the matched-set sandwich covariance transformed to the additive basis. Each fit uses Holm correction across all candidate exposures at a significance level of 0.05. The active set includes relative humidity, which participates through an interaction with temperature in the generating response. The CLR fits use additive exposure terms. Correctly adjusted CLR fits supply the G1 and G2 comparison, and misspecified adjustment supplies the G3 and G4 comparison. C uses G1 and correctly adjusted CLR fits.

The analysis evaluates the proportion of active exposures identified, exposure-specific identification, and the proportion of null exposures selected. Repetition-level rates provide the Monte Carlo intervals for overall performance. Exposure-specific intervals use the Wilson method. Paired differences use common usable repetitions and are reported in percentage points. Fit and test availability are recorded for each planned repetition.

## Publication report

[results/index.html](results/index.html) links the figures, tables, and captions. The two figures report overall active- and null-exposure identification and individual active-exposure identification across the four dimensions in C. The tables provide the B, C, and D comparisons, exposure-specific rates, and paired differences. `results/results_manifest.csv` provides sizes and SHA-256 digests for the publication report.

## Analysis source and records

`RUN_SELECTION_COMPARISON.R` is the submission entry point. `R/selection.R` contains the matched simulation plan, comparator fits, and statistical tests. `R/report.R` contains the aggregation and reporting procedures. `tools/package_results.py` assembles the publication report and release assets. The analysis uses the Stage 1 run archive, R 3.6.3 or later with CompQuadForm, and Python 3.6 or later for archiving.

The supplemental release assets contain the task plan, input-run references, execution log, checkpoint records, fit availability, repetition-level tests, and the 800 planned comparator fits for C. `release_assets.csv` and `runtime_manifest.csv` record sizes and SHA-256 digests. The original fitted objects and evaluation datasets are indexed by [the Stage 1 run archive](../GitHub_repository/docs/ARCHIVE_LAYOUT.md).
