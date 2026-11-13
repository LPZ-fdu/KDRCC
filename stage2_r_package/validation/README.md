# Independent validation of kdrcc 0.1.0

This project validates the installed R package against the formal first-stage KDR-CC G1 study. The principal reproduction uses the same 100 independent repetition IDs, generating law, learning and estimation samples, kernel-learning streams, auxiliary folds, evaluation points and statistical summaries as the archived Module B results.

## Run the matched formal reproduction

Install the package with `../package/INSTALL_AND_TEST.R`, then run:

```r
source("RUN_END_TO_END.R", encoding = "UTF-8")
```

`END_TO_END_SETTINGS.R` sets `replicates = 1:100`, `workers = 32L`, checkpoint resumption and report appearance. All statistical settings come from the frozen first-stage configuration: N = 20,000 estimation matched sets, L = 20,000 independent learning matched sets, M = 144 centers, kappa = 0.0001 and lambda = 1. The original configuration identifier is A035.

The installed package independently learns each kernel and fits each final model. Archived coefficients supply the comparison target after estimation. The adapter synchronizes the center and penalty-row seeds separately through the installed package's learning and penalty functions, then uses the public `kdrcc()` fitting interface. The G1 auxiliary dictionary and original matched-set fold IDs are carried explicitly. Monte Carlo fallback variable tests use the original variable-specific streams. Details appear in `END_TO_END_README.md`.

## Final outputs

`results/LATEST_G1_REPRODUCTION.txt` identifies the output directory. Its `analysis` directory contains:

- `index.html`: a gallery of the final G1 comparison figures and tables;
- `figures`: joint and single-exposure curves, observed-combination predictions, interactions, contribution shares, quantitative ORs, OR errors and active-variable identification;
- `tables`: aggregated agreement and performance results in CSV and titled three-line HTML tables;
- `model_diagnostics`: summarized numerical diagnostics and their comparison graphic;
- `summary.rds`: the aggregated data needed to regenerate the figures;
- `REPORT.md` and `Figure_and_table_captions.md`: execution summary and output definitions.

The archived and package curves average the same selected repetition IDs. Light blue bands average fitted pointwise 95% interval endpoints; darker blue bands show the central 95% between-repetition range. Both implementations use the same generating truth and reference combinations.

Per-repetition fits, learned kernels and worker logs are temporary computation checkpoints. They are removed after the complete cohort is successfully summarized and rendered. The final aggregated data and reports remain. To redraw final figures independently:

```r
source("PREPARE_RESULTS.R", encoding = "UTF-8")
```

## Numerical and installation records

The 26-check numerical suite validates moments, exact curvature, auxiliary-adjusted covariance, derivatives, serialization, contribution allocation, testing, tuning and saved formal outputs. Execute it with `source("scripts/01_validate.R", encoding="UTF-8")` from this directory. `RUN_VALIDATION.R` runs this numerical suite followed by the full matched reproduction. The recorded installation workflow is retained in `records/installation_test`. The completed formal reproduction is retained in `results/G1_reproduction_0a040df3a874`. Its 100 matched repetitions record 2,900/2,900 numerical checks passing under Linux R 3.6.3.

`reference/stage_one_G1` contains the input comparison archive: 100 frozen formal records and kernels, the original evaluation distribution, configuration, task plan, filtered G1 summaries and source checksums. These are reference inputs to the reproduction. `reference/formal_simulation` and `reference/formal_report` preserve the generating and summary definitions; `reference/frozen_engine` supports numerical checks. RDS inputs and outputs use serialization version 2. The package supports R >= 3.6.3; platform-specific executions record their R version and session information.
