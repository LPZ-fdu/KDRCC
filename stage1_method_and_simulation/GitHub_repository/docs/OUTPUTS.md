# Output reference

A module run contains `config.rds`, `plan.rds`, `manifest.rds`, task/configuration CSV files, session information, source snapshots, kernels, nuisance estimates, fit records, evaluation coordinates and logs. `worker_processes.csv` records active PSOCK PIDs. Result directories are identified by the scientific configuration and implementation hashes.

## Shared tables

Every module reports availability, configuration, joint/observed-combination effect accuracy, absolute OR relative errors, interval calibration, variable identification and computational resources. HTML tables have titles, named columns and explanatory notes; CSV versions support downstream analysis.

Interpretation tables include:

- `Table_<module>_contribution_shares`: true and estimated shares, bias, RMSE, SD and MCSE.
- `Table_<module>_quantitative_OR`: selected contrasts, true OR, mean/median estimated OR, mean interval limits, errors, coverage and direction detection.
- `Table_<module>_interaction_detection`: overall raw and pair-family Holm detection, usable tests and Wilson Monte Carlo intervals.
- `Table_<module>_interaction_direction`: nominal and full-grid-Holm correct/wrong direction rates and Monte Carlo error.
- `Table_<module>_interpretation_metrics`: configuration-level inputs to the extended assessment score.

Repetition-level tables retain the underlying contributions, contrasts and decisions. Single-exposure point tables contain both repetition quantiles and mean CI endpoints with interval counts.

## Figures

A provides one-dimensional parameter/sample-size plots and two-dimensional parameter grids for effect accuracy and interpretation accuracy. Computational costs accompany the effect plots. `Table_A_configuration_ranking` and `recommended_configuration.rds` support reviewed parameter freezing.

B provides joint-effect recovery, observed-combination prediction, OR-error distributions, single-exposure curves G1–G4, interaction-contrast recovery, contribution-share stacked bars, quantitative OR forests and interaction significance/direction recovery. C provides the same kinds of effect and interpretation outputs for G1 KDR-CC across null counts, plus identification and false-positive trends. C single/interaction recovery figures are indexed by null count.

Confidence-limit shading and repetition-range shading have different meanings, specified in the captions. Quantitative OR forest segments show the average interval endpoints; their full error and coverage summaries are in the accompanying tables.

## Model diagnostics

`analysis/model_diagnostics/` contains:

- Fit completion, stopping gradients, equation residuals, iteration counts, objective values, Jacobian condition numbers, covariance eigenvalues and errors.
- Iteration histories and convergence figures by method and nuisance scenario.
- Kernel-learning projected gradients, optimizer statuses and distance-parameter boundary rates.
- Formatted tables, a local caption registry and an index page.

## Figure and table captions

`analysis/Figure_and_table_captions.csv` and `.txt` consolidate every figure/table title and note, including the diagnostics subdirectory. `analysis/index.html` provides a browsable report. The reporting source and effective font/settings are copied into `analysis/report_source/` and `report_metadata.rds` for reproducibility.

## Module D correlation outputs

`Figure_D_correlations.png` displays the original B, moderate and high measured overall Pearson matrices. B fits are reused; only the moderate and high settings generate new fits. `Table_D_correlations` records target means, pairwise floors, achieved coefficients, sources and pairwise ranges; `Table_D_pairwise_correlations` records all pairs. `Table_D_performance` links every performance row to its correlation setting and achieved overall correlation.

Trend figures compare joint and observed-combination RMSE, RMS bias, OR error, SE calibration, coverage and significant effect directions across correlation levels for three methods under G1–G4. Other figures summarize contribution-share error, interaction detection and G1 exposure identification. `analysis/correlation_scenarios/<ID>/` contains the full within-setting joint, actual-combination, single-exposure, contribution, quantitative-OR and interaction displays and their tables. The main caption registry and HTML index link these nested reports. Numerical diagnostics remain in `analysis/model_diagnostics/`.

## Module E selection and paired evaluation

`choices/` stores selected configurations and ranked validation scores for every independent repetition. `cross_validation/` stores candidate-by-fold scores and failures. `Table_E_performance` summarizes the full adaptive procedure alongside the fixed configuration; `Table_E_paired_changes` reports within-repetition effect-metric differences and Monte Carlo errors. Other tables contain interval calibration, quantitative ORs, contribution shares, interaction tests and directions, variable identification, selection frequencies, fallback availability and time costs.

Figures compare accuracy, interval calibration, overall and individual exposure identification, interpretation metrics, configuration frequencies and selection time. Scenario-specific figures display joint curves and observed-combination predictions. `Table_E_joint_support` records whether joint evaluation points lie within a prespecified empirical exposure ellipsoid. This support summary uses the eight active coordinates. Numerical diagnostics are stored in `analysis/model_diagnostics/`. The common caption registry and HTML index list all generated outputs.


## Saved-model contribution evaluation

Contribution shares evaluate the eight prespecified active exposures at the reference background, normalized within those eight variables. Use `05_submit.R contributions C,E` followed by `05_submit.R report C,E` to calculate contributions from saved fitted functions and write the standard reports. Execution settings remain in `05_submit.R`. See [Contribution evaluation](CONTRIBUTIONS.md).
