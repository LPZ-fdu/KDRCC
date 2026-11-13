# kdrcc

**Kernel-Based Doubly Robust Case-Crossover Regression (KDR-CC)**

Author and maintainer: Peizheng Li (李培政), pzli25@m.fudan.edu.cn. Version 0.1.0. License: GPL-3. R >= 3.6.3.

The package fits nonlinear multiexposure associations in individually matched case-crossover data. It combines independent Gaussian-kernel learning, matched-set auxiliary cross-fitting and a regularized criterion built from a doubly robust moment. It provides relative effect curves, quantitative odds ratios (ORs), pointwise intervals, derivative-based participation tests, four-point interaction tests, exact baseline Shapley allocations and grouped configuration selection.

## User guides

The distribution includes complete Chinese and English Word guides in the parent package directory, with runnable case data and example scripts in `../examples`.

## Installation

Copy the complete release directory to the target computer. In R, set the working directory to the distribution's `package` directory and run:

```r
source("INSTALL_AND_TEST.R", encoding = "UTF-8")
```

The script installs dependencies, installs the package, runs an independent matched-data example and writes an installation report, a fitted model and PNG figures to `installation_test`. `settings` at the beginning controls the library, archive type and dependency installation. Windows R 4.4 can use the included binary archives. Other supported R versions install the pure-R kdrcc source archive. CompQuadForm source installation uses the platform's C/C++ build tools; its source archive is included in `installers`.

An installed CompQuadForm >= 1.4.4 is reused. The script first uses a compatible bundled Windows binary, then the configured CRAN repository, and finally the included dependency source archive. A network-connected computer can install CompQuadForm through its standard R repository. On Linux the source archive builds with a C/C++ compiler and R development tools. On Windows, source dependency compilation uses the Rtools release matching the installed R version.

## Data and a complete fit

Each row represents an event or reference date. `strata` identifies the individual's matched set; `case` contains one event indicator and zero-valued reference indicators. Exposure columns are numeric. A one-sided adjustment formula specifies time-varying covariates. Complete matched sets define all data partitions and folds.

```r
library(kdrcc)
d <- kdrcc_example(800, seed = 702)
parts <- kdrcc_split(d, "set_id", learning = 250, seed = 703)
control <- kdrcc_control(centers = 10, kappa = 0.001, lambda = 1,
                        kernel_tolerance = 0.0001)
fit <- kdrcc(parts$estimation,
             exposures = c("X1", "X2", "X3", "X4"),
             adjust = ~ W + I(W^2), learning_data = parts$learning,
             nuisance = kdrcc_nuisance(g_df = 2, kernel_rank = 4, f_df = 2),
             control = control)
summary(fit)
predict(fit, parts$estimation[1:5, ], type = "or")
joint <- kdrcc_joint(fit)
single <- kdrcc_single(fit)
plot(joint)
plot(single)
```

The example contains harmful X1, protective X2, an X1–X3 interaction, a null X4 and time-varying W. Increasing and decreasing associations use the same unconstrained response family. Concentration curves and specified OR contrasts describe effect direction and magnitude.

With `learning_data` omitted, `kdrcc()` reserves 25% of matched sets for independent learning. For a fixed learned kernel, pass `kernel = previous_fit$kernel`. Its exposure and auxiliary designs are then reused. `reference` sets a named complete exposure combination; prediction also accepts row-specific reference combinations.

## Auxiliary models

`g` is estimated by conditional logistic regression using the adjustment design and a finite exposure dictionary. Its adjustment coefficients define the ordered-pair covariate function. `g_exposure` accepts `"spline_kernel"`, `"spline"`, `"kernel"`, a one-sided formula, or a finite-matrix callback.

`f` estimates the conditional kernel-feature difference among event-second ordered pairs. `f_method = "basis"` uses finite weighted least squares on ordered covariate spline terms. `f_formula` supports a prespecified finite design with names such as `plus.W` and `minus.W`. `f_method = "cells"` uses saturated indicators for discrete ordered covariate combinations. Each training fold requires information for the chosen finite design. Spline knots, column selection and whitening are learned independently. The covariance propagates the full g and f estimation influence through the final criterion.

```r
nu <- kdrcc_nuisance(g_exposure = ~ X1 + I(X1^2) + X2 + X3,
                    f_method = "basis",
                    f_formula = ~ plus.W + minus.W + I(plus.W^2) + I(minus.W^2))
```

## Configuration selection

Use three independently allocated samples: kernel learning, configuration selection and final estimation. The selection sample uses outer matched-set validation, with auxiliary cross-fitting inside each training fold. The score is held-out conditional negative log likelihood per matched set.

```r
d <- kdrcc_example(1000, seed = 85)
parts <- kdrcc_split(d, "set_id", learning = 250, tuning = 400, seed = 86)
nu <- kdrcc_nuisance(g_df = 2, kernel_rank = 4, f_df = 2)
grid <- expand.grid(M = c(8, 12), kappa = c(0.001, 0.003), lambda = c(0.3, 1))
cc <- kdrcc_control(centers = 8, kappa = 0.001, lambda = 1, workers = 2)
selection <- kdrcc_tune(parts$tuning, parts$learning,
  exposures = c("X1", "X2", "X3", "X4"), adjust = ~ W + I(W^2),
  nuisance = nu, grid = grid, control = cc)
fit <- kdrcc(parts$estimation, exposures = c("X1", "X2", "X3", "X4"),
  adjust = ~ W + I(W^2), kernel = selection$kernel, control = selection$control)
```

`workers` controls independent kernel-profile tasks in selection. Each final fit uses one process. A candidate must produce finite held-out scores and a usable sandwich covariance from every training fold. The tuning object stores fold losses, covariance diagnostics, eligibility and the selected configuration. If all candidates fail and `fallback = TRUE`, the package warns and returns the supplied control configuration. Final estimation checks covariance again on the independent estimation sample.

## Interpretation

```r
participation <- kdrcc_variables(fit)
interaction <- kdrcc_interaction(fit, pairs = matrix(c("X1", "X3"), 1, 2))
contribution <- kdrcc_contribution(fit, variables = c("X1", "X2", "X3"))
diagnostics <- kdrcc_diagnostics(fit)
saveRDS(fit, "fitted_kdrcc.rds", version = 2)
```

Variable participation uses the evaluation-average squared derivative of the regularized working response. The working null is a zero derivative on that grid; weighted chi-square tail probabilities and Holm-adjusted p values are returned with numerical stability information. The finite positive-distance Gaussian representation and the grid rank determine the null's geometry. Generating-label recognition is assessed by simulation in the development archive.

Interaction uses the anchored four-point log-OR contrast. Output includes overall pair tests and pointwise effect estimates, intervals, adjusted p values and directions. Pair tests and pointwise tests form separate multiplicity families.

Baseline Shapley allocation decomposes a specified exposure contrast relative to a complete reference. `variables` identifies the contribution set. The remaining coordinates stay at reference. Signed allocations sum to the corresponding log OR; overall shares normalize mean absolute allocations within that set.

## Statistical and computational target

The base moment has the double-robust population identity under the logistic pair model and either correct auxiliary function. The finite Gaussian representation supplies the response coordinates. Final estimation minimizes one half of the squared transformed moment plus a mixed-derivative quadratic penalty. Positive `lambda` defines a regularized population target. The reported covariance uses exact local criterion curvature and finite-dimensional auxiliary influence, conditional on independently learned and selected objects. Wald intervals and quadratic-form references use the resulting local asymptotic approximation. See `inst/doc/method.html` and the help pages for definitions and assumptions.

`centers` gives the candidate kernel center count, `kappa` stabilizes independent kernel learning and `lambda` controls mixed derivative complexity in final estimation. Learning retains the kernel features with usable within-set variation and a numerically independent column span. The retained dimension and per-feature status appear in `kdrcc_diagnostics()`. Package defaults provide starting settings for grouped validation. Dense feature and auxiliary arrays scale with matched dates, centers and auxiliary dimension. `keep_data = FALSE` retains fitted components and evaluation rows.

## Reproducibility

Public objects retain their specification, seeds, matched-set folds, learned geometry, auxiliary fits, convergence history and software version. Randomized package operations restore the caller's random-number state. `../../validation` contains executable numerical comparisons and actual execution records. Standard package tests run through `R CMD check --no-manual`. The source repository includes Windows and Linux CI checks.

## Citation and license

```r
citation("kdrcc")
```

Peizheng Li (2026). *kdrcc: Kernel-Based Doubly Robust Case-Crossover Regression*. R package version 0.1.0. Distributed under GNU General Public License version 3.
