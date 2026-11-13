# Active-exposure contribution evaluation

The contribution target contains the eight prespecified active exposures. All additional coordinates are fixed at the marginal-median reference throughout the Shapley coalitions. Mean absolute contributions are normalized within the eight active exposures. Fitted functions retain all candidate coordinates used during estimation.

## Saved-model execution

Set `code_root`, `work_root` and `workers` in `05_submit.R`. The default is 32 PSOCK processes. The `contributions` action reads existing fitted coefficients, kernel objects, evaluation combinations and references. It writes the contribution fields into each corresponding `fits/<key>.rds` record. Other fitted and inference fields retain their values. Each completed record is committed before the next task; matching completed contributions are reused on a repeated submission.

From the source directory:

```bash
nohup Rscript --vanilla 05_submit.R contributions C,E > CE_contributions.txt 2>&1 &
```

After this command completes:

```bash
nohup Rscript --vanilla 05_submit.R report C,E > CE_report.txt 2>&1 &
```

`06_prepare_contributions.R` is a convenience entry for the same calculation and defaults to C,E. `report_run_directories` in `05_submit.R` can specify complete paths to selected runs. Otherwise, the latest requested module runs are resolved from the work directory.

Progress records include completed counts, elapsed time and estimated remaining time. `contribution_processing.csv` contains per-record status; `contribution_processing.rds` identifies the evaluation scope and code. The source snapshot contains the contribution and reporting implementation. The report action writes the standard analysis tables, figures and captions in place. Completed source and result directories each contain one contribution definition and one set of reported contribution results.

## Evaluation quantities

For active set A and reference b, the allocated contrast is `h(x_A,b_other)-h(b)`. The eight Shapley values sum to this contrast at every evaluation point. The share for exposure j is `100*mean(abs(phi_j))/sum_{k in A}mean(abs(phi_k))`. Contribution error is computed over these eight shares in every repetition, then summarized across repetitions. Variable identification uses the full candidate set; joint and observed-combination effect assessments use their specified prediction inputs.

Four-node Gaussian-kernel quadrature evaluates the eight-variable game exactly up to floating-point error after absorbing the fixed-background factors into the center coefficients. Algebraic tests compare this calculation with subset enumeration. Saved-model tests verify normalization, numerical closure, field preservation, parallel execution and resumability.
