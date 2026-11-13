# KDR-CC simulation study

Reproducible simulation code for Kernel-Based Doubly Robust Case-Crossover Regression (KDR-CC). The study evaluates effect estimation, uncertainty, exposure identification, contribution allocation and interaction detection in individual time-stratified case-crossover data.

## Frozen stage-one materials

This repository contains the KDR-CC method specification, Chinese appendix, simulation implementation and A–E assessment records. It supplies the development-stage evidence and saved-model reference for the companion R package and package-validation materials. The stage identifier is `stage1-20261107`.

The complete publication reports are in `records/results/`. Read each module's `analysis/index.html` for its figures, tables, captions and model diagnostics. Actual submission settings and the reviewed configuration are in `records/submissions/` and `records/config/`.

Detailed fitted records, learned kernels, evaluation objects, cross-validation records and execution logs are distributed as independent ZIP assets attached to the same GitHub Release. Download all `stage1_runtime_*.zip` files and follow [ARCHIVE_LAYOUT.md](docs/ARCHIVE_LAYOUT.md). The file manifests in `archive/` link every retained runtime file to its asset and SHA-256 digest.

```bash
python tools/verify_stage1.py
python tools/verify_stage1.py --assets ../release_assets
python tools/restore_records.py --assets ../release_assets --destination ..
```

The restored `KDRCC_formal_20261027` directory sits beside the source repository. Submission and reporting use `05_submit.R`. [FINAL_REVIEW.md](docs/FINAL_REVIEW.md) states the statistical scope and review results. [PUBLICATION_METADATA.md](docs/PUBLICATION_METADATA.md) lists the publication and licensing fields maintained by the project owners.

## Requirements and execution

Use 64-bit R 3.6.3 or newer with `CompQuadForm`, `parallel` and `splines`. The CRAN source archive of CompQuadForm is provided in `vendor/` for offline installation. Linux source installation requires a C/C++ toolchain. PNG reporting uses Cairo and Times New Roman when available, with the effective font recorded in report metadata.

**Edit `05_submit.R` for every submission setting.** Set `code_root` to this source directory and `work_root` to a separate writable directory. With both set to NULL, the script locates its source directory and uses the sibling `KDRCC_formal_20261027` directory. Source and work directories are separate and non-nested. The default is 32 independent PSOCK R processes with one BLAS/OpenMP thread each.

From the directory containing `05_submit.R`:

```bash
Rscript --vanilla 05_submit.R install
Rscript --vanilla 05_submit.R check
Rscript --vanilla 05_submit.R plan A
nohup Rscript --vanilla 05_submit.R run A > A_run.txt 2>&1 &
```

After fitting finishes, prepare A's figures and tables:

```bash
nohup Rscript --vanilla 05_submit.R report A > A_report.txt 2>&1 &
```

Review the results. Set `submission$freeze$config_id` to the chosen configuration ID, or keep NULL to use the eligible recommendation. Keep `freeze$source_run=NULL` to use the latest A run in the work directory. Then:

```bash
Rscript --vanilla 05_submit.R freeze
Rscript --vanilla 05_submit.R plan B,C
nohup Rscript --vanilla 05_submit.R run B,C > BC_run.txt 2>&1 &
```

Reporting is a separate action. B can be reported while C is running:

```bash
nohup Rscript --vanilla 05_submit.R report B > B_report.txt 2>&1 &
nohup Rscript --vanilla 05_submit.R report C > C_report.txt 2>&1 &
```

The frozen file is `<work_root>/config/fixed_parameters.R`. `common_fixed` overrides its parameters for B, C and D; `B_overrides`, `C_overrides` and `D_overrides` take precedence over the common settings. All five parameters can also be provided explicitly. The generating-mechanism identity is checked when reading a frozen file. Settings, overrides, code hashes, logs and session information are recorded for every submission.

## Study organization

| Module | Purpose | Default independent repetitions and fits |
|---|---|---|
| A | Parameter configuration, sample size, contribution accuracy, interaction detection and computational cost | 50 unique configurations × 100 repetitions: 5,000 fits |
| B | Fixed configuration; KDR-CC, Spline-CLR and Linear-CLR under G1–G4 | 8 distinct fits × 100 repetitions: 800 fits |
| C | G1 KDR-CC with 8, 16, 40 or 80 null exposures | 4 dimensions × 100 repetitions: 400 fits |
| D | Fixed configuration at moderate/high overall correlation, with original B reference reused | 2 new settings × 8 distinct fits × 100 repetitions: 1,600 new fits |
| E | Data-driven selection versus the frozen configuration across seven settings | 700 paired repetitions: 37,800 validation fits and 1,400 final fits |

The full A–D default plan contains 7,800 new fits; modules are submitted independently. CLR uses two adjustment specifications per method; their outputs are displayed under the corresponding four nuisance scenarios. A uses an independent simulation batch from B/C. Within a repetition, candidate dimensions share the active-exposure data and existing null columns. Every repetition generates new learning and estimation data.

A evaluates joint and observed-combination errors, model SE versus empirical variability, interval coverage, contribution-share errors, overall interaction detection, direction detection and resource use. A's explicit weighted recommendation is reviewed before freezing. B and C calculate the same effect and interpretation targets. C additionally measures true-positive and false-positive variable-identification rates.

## Results

Each run is stored under `<work_root>/results/module_<module>_<identifier>/`. The `analysis/index.html` page links publication figures, formatted tables and their captions. Model fitting diagnostics have their own `analysis/model_diagnostics/` directory. Tables include CSV data and formatted HTML with titles and notes. Figures are PNG files. `Figure_and_table_captions.csv` and `.txt` provide the shared caption registry.

Single-exposure figures show mean fitted curves, light shading from separately averaged pointwise 95% interval limits, and darker shading from the 2.5th–97.5th percentiles of replicate point estimates. Contribution shares are displayed as stacked bars. Quantitative OR contrasts are displayed as forest plots and tables. Interaction reports distinguish overall pair detection from correct significant directions at individual grid points.

The source tree contains implementation, configuration, documentation, numerical fixtures and executable checks. Generated fits, runtime records and analyses reside in the work directory. To regenerate figures locally, copy a complete module run and set `report_run_directories` to that path in `05_submit.R`.

- [Simulation protocol](docs/SIMULATION_PROTOCOL.md)
- [Module D: design, submission and outputs](docs/MODULE_D.md)
- [Module E: data-driven configuration selection](docs/MODULE_E.md)
- [Parameter configuration and recommendation](docs/PARAMETERS.md)
- [Output definitions](docs/OUTPUTS.md)
- [Server execution](docs/SERVER_RUN.md)
- [Reproducibility](docs/REPRODUCIBILITY.md)
- [Implementation checks](docs/IMPLEMENTATION_CHECKS.md)
- [R compatibility](docs/R_COMPATIBILITY.md)

## Overall-correlation assessment

Module D reads the reviewed configuration from the existing work directory. `D_correlation_levels` and `repetitions$D` in `05_submit.R` control its two new correlation settings and independent repetitions. The original-correlation reference reuses completed B fits. Moderate and high mean absolute Pearson targets are 0.60 and 0.85, with pairwise floors of 0.45 and 0.75; the mean tolerance is 0.03. Fixed empirical marginal maps and Gaussian scores calibrate the final lag-defined exposure dependence. Temporal and exposure–confounder dependence also change under this recoupling. `D_overrides` permits explicit parameter changes. D has its own random-stream namespace and scenario-specific caches. The full correlation matrices and achieved mean absolute Pearson coefficients accompany the performance comparisons.

```bash
nohup Rscript --vanilla 05_submit.R plan D > D_plan.txt 2>&1 &
nohup Rscript --vanilla 05_submit.R run D > D_run.txt 2>&1 &
```

Review the plan before submitting the run. After completion:

```bash
nohup Rscript --vanilla 05_submit.R report D > D_report.txt 2>&1 &
```

Use the same `work_root` that contains the reviewed `config/fixed_parameters.R`. The D submission checks the completed B reference and copies it into the D run. The report combines all three correlation settings, with full three-method G1–G4 assessments in `analysis/correlation_scenarios/`.

## Application-oriented configuration selection

Module E compares the frozen configuration with three-fold selection by held-out conditional negative log likelihood. Every repetition has independent learning, configuration-selection and final-estimation samples. The two strategies share the final data. E uses training-data spline and kernel nuisance features. Seven default settings cover eight active exposures, 8/16/40/80 additional null exposures, and moderate/high overall correlation. Each setting has 100 new independent repetitions. Candidate configurations, all validation scores, selected settings and paired final outputs are retained.

`05_submit.R` starts with module A and 32 workers. Select E explicitly on the command line after freezing the reviewed configuration. `E_settings` controls the candidate set, folds and scenarios; `E_overrides` changes inherited fixed parameters. Review `plan E`, submit `run E`, and then run `report E` separately. Details are in [MODULE_E.md](docs/MODULE_E.md). The default 18-candidate plan requires 37,800 validation fits and 1,400 final fits; validation and final-fitting progress are shown independently.

## Method specification and release

The statistical definitions and implementation map are in [METHOD_SPECIFICATION.md](docs/METHOD_SPECIFICATION.md). The complete Chinese theory and simulation appendix is in [manuscript/Appendix_KDRCC_zh.docx](manuscript/Appendix_KDRCC_zh.docx), with its structured source and reproducible document builder beside it. [RELEASE_VALIDATION.md](docs/RELEASE_VALIDATION.md) records the review findings and source-to-results inventory. `SOURCE_CHECKSUMS.csv` freezes the distributed files.

`config/fixed_parameters.R` contains the reviewed A035 settings. For a complete new A–E study, use the review-and-freeze workflow above. To reproduce B–E directly with that reviewed configuration, copy this file to `<work_root>/config/fixed_parameters.R`; the run-specific overrides remain available. D additionally uses the completed B reference. Submission records preserve the effective settings in both cases.

The base moment has a double-robust zero under a correct finite response representation and either correct nuisance model. The implemented fixed-penalty estimator targets the corresponding regularized moment criterion. The specification separately states these targets, the conditional covariance, and the simulation comparisons against the generating function.


## Saved-model contribution evaluation

Contribution shares evaluate the eight prespecified active exposures at the reference background, normalized within those eight variables. Use `05_submit.R contributions C,E` followed by `05_submit.R report C,E` to calculate contributions from saved fitted functions and write the standard reports. Execution settings remain in `05_submit.R`. See [Contribution evaluation](docs/CONTRIBUTIONS.md).
