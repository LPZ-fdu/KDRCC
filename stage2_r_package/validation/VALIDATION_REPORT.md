# Independent validation records

Package: kdrcc 0.1.0. Maintainer: Peizheng Li, pzli25@m.fudan.edu.cn. GPL-3.

## Completed formal reproduction

`results/G1_reproduction_0a040df3a874/analysis/REPORT.md` records 100/100 original G1 repetitions completed and 2,900/2,900 numerical comparisons passing. The executed environment is Linux x86_64, Ubuntu 18.04.5 LTS, R 3.6.3, kdrcc 0.1.0 and CompQuadForm 1.4.4. The fixed configuration is N=20,000, L=20,000, M=144, kappa=0.0001 and lambda=1.

The installed package independently regenerates the formal data, learns the kernel and estimates the model using the original named random streams. Comparisons cover kernel components, coefficients, full covariance, joint and observed-combination effects, single-exposure curves, interaction contrasts, variable statistics and decisions, contribution shares, OR contrasts, interaction significance and direction decisions.

The analysis directory contains the execution report, final gallery, PNG figures, aggregated agreement and performance tables, captions, model diagnostics and summary.rds. Input comparison records are retained in reference/stage_one_G1. The source manifest and checksum files identify the inputs. `records/installation_test` retains the installation workflow record.

## Executable checks

`scripts/01_validate.R` implements 26 numerical and saved-output checks covering moments, Jacobians, coefficients, exact curvature, full auxiliary covariance, derivatives, serialization, Shapley allocation, Holm adjustment, interaction anchors and serial/PSOCK tuning. `scripts/03_check_source.R` runs standard package build and check commands. Package-level test scripts are supplied in `../package/kdrcc/tests`.

`RUN_VALIDATION.R` runs the numerical suite and matched reproduction. `RUN_END_TO_END.R` runs the matched cohort directly. `PREPARE_RESULTS.R` renders the final plots from the aggregated summary. Details of inputs, tolerances and execution appear in END_TO_END_README.md and END_TO_END_SETTINGS.R.
