# Matched reproduction of formal KDR-CC G1 results

## Scope and configuration

The study reproduces KDR-CC G1 from formal Module B `module_B_6fab5de7950a`. The default cohort is original repetition IDs 1 through 100. Each repetition has 20,000 independent learning matched sets and 20,000 estimation matched sets, 144 Gaussian kernel centers, kernel-learning penalty kappa = 0.0001 and interaction-complexity penalty lambda = 1. The full specification is stored in `reference/stage_one_G1/configuration.rds`; the selected task rows are in `plan.rds`.

## Exact matching of the study inputs

1. **Matched data.** The archived generating functions regenerate each original 2,000-set block with its named `data_pool` learning or estimation stream. Matched-set IDs are offset when blocks are combined. Original cross-fitting fold labels are preserved.
2. **Auxiliary models.** The adjustment is `~ W + I(W^2)`. The G1 exposure dictionary consists of the six prespecified pollutant anchor functions, cold and hot temperature squares, a humidity square, standardized exposures and all pairwise products. The package estimates the g coefficients by cross-fitting. The f model estimates saturated ordered W-cell means from the training folds.
3. **Independent kernel learning.** The adapter calls the installed package's `.kdr_validate()` and `.kdr_learn()` with the original named center seed. It then calls the installed package's `kdr_kernel_penalty()` with the separate original named penalty seed. This explicitly supplies the two sampling streams used in the formal implementation. Kernel geometry, learned distances, feature whitening and moment transformation are learned from the regenerated sample. The archived fitted coefficients supply comparison values after fitting.
4. **Final fitting.** Public `kdrcc()` receives the independently learned kernel and original matched-set fold IDs. It estimates auxiliaries, regularized coefficients and the full coefficient covariance. Each parallel worker fits one repetition with numerical-library thread counts set to one.
5. **Variable tests.** Public `kdrcc_variables()` computes derivative quadratic statistics and weighted chi-square tails. For Monte Carlo fallback tails, the adapter supplies the original named variable-test stream using the installed package tail routine, then applies the same Holm adjustment and numerical-stability rule.
6. **Response evaluation.** Public prediction, single-exposure, interaction and contribution interfaces evaluate the archived reference vector, concentration grids, observed combinations, four interaction pairs, quantitative OR contrasts and active-variable contribution distribution.
7. **Comparison.** Each completed refit is compared with its corresponding archived result. Checks cover kernel components, coefficients, covariance, predictions, standard errors, contributions, variable decisions, interaction decisions and directions. `END_TO_END_SETTINGS.R` specifies numerical comparison tolerances. Summary tables report the number of checks passing and maximum discrepancies across the cohort.

The adapter is versioned for installed package 0.1.0. `source_manifest.csv` identifies the original fitted records and kernels by filename and MD5. `VALIDATION_CHECKSUMS.csv` identifies the validation inputs and scripts by SHA-256.

## Execution

```r
source("RUN_END_TO_END.R", encoding = "UTF-8")
```

Default execution uses 32 PSOCK R processes. Adjust `workers` in `END_TO_END_SETTINGS.R` to match the computer. Set `replicates` to original repetition IDs when selecting a smaller cohort; both archived and package summaries then use that exact same cohort. Keep `replicates = 1:100` to reproduce the full formal mean curves.

The console reports completion count, elapsed time, estimated remaining time and current stages. Interrupted runs resume saved checkpoints. `prepare_results_after_fit = TRUE` produces final reports after fitting. With this setting set to `FALSE`, run the independent result entry when ready:

```r
source("PREPARE_RESULTS.R", encoding = "UTF-8")
```

`require_complete_cohort = TRUE` requires valid estimation, intervals, response outputs, variable tests and interaction tests for every selected repetition before generating final comparison figures. Incomplete runs retain their computation checkpoints for resumption. Successful complete reports remove temporary per-repetition checkpoints and retain final aggregated results. Subsequent calls to either entry reuse the completed summary. Use a new output directory or `resume = FALSE` for a new complete run.

## Final figures and tables

The figures show generating truth, archived KDR-CC G1 and installed-package KDR-CC G1. They include the joint curve, eight single-exposure curves, predictions at observed combinations, interaction contrasts, variable identification, overall interaction detection and correct significant directions, stacked contribution shares, quantitative ORs and absolute relative OR errors. Numerical diagnostics have a dedicated directory.

Mean confidence bounds average the fitted lower and upper pointwise limits separately. Repetition bounds are empirical pointwise 2.5th and 97.5th percentiles. Contributions average fitted percentage shares over repetitions. Tables report RMS Bias, RMSE, empirical SD, model SE, OR error, interval performance, detection and numerical agreement using the original study's aggregation definitions.

`analysis/summary.rds` supports independent plot regeneration. `analysis/index.html` provides the result gallery. No individual-repetition figures or individual fitted models form part of the final output.
