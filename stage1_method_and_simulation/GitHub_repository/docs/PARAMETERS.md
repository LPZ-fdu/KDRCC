# Parameters and module A assessment

All run-specific choices belong in `05_submit.R`. Defaults are specified in `00_parameters.R`; the complete 50-configuration design is in `default_module_A_configurations.csv`, with view memberships in `default_module_A_views.csv`.

N is the number of estimation matched sets; L is the number of independent learning matched sets; M is the number of Gaussian kernel centers. Kappa penalizes coefficients during outcome-supervised kernel learning. Lambda multiplies the normalized mixed-derivative penalty during final coefficient estimation. The module A grid reference is N=20,000, L=20,000, M=144, kappa=0.0001 and lambda=10. The reviewed configuration A035 used by B–D is N=20,000, L=20,000, M=144, kappa=0.0001 and lambda=1; E uses it as its fixed comparator.

A varies N from 1,000 to 50,000, L from 2,000 to 50,000, M from 24 to 192, kappa from 0.000001 to 0.001, and lambda from 0 to 20. Two-dimensional grids compare kappa with lambda and L with M. Duplicate configurations are fitted once and shared between the appropriate figures. `A_overrides` can change reference values, grid values or the entire `manual_design` data frame.

## Primary recommendation score

Each term is averaged or aggregated across independent repetitions as specified below; smaller scores rank first. The weights are explicitly editable through `A_overrides$recommendation$weights`.

| Term | Definition | Weight |
|---|---|---:|
| joint_nrmse | Joint RMSE / RMS true joint contrast | 0.20 |
| individual_nrmse | Observed-combination RMSE / RMS true contrast | 0.20 |
| joint_nrms_bias | Joint RMS bias / RMS true contrast | 0.10 |
| individual_nrms_bias | Observed-combination RMS bias / RMS true contrast | 0.10 |
| joint_or_error | Mean joint absolute OR relative error as a fraction | 0.05 |
| individual_or_error | Mean observed-combination absolute OR relative error as a fraction | 0.05 |
| contribution_share_error | Mean across repetitions of RMS share error over the eight prespecified active exposures, in fractional units | 0.10 |
| interaction_miss | One minus the mean Holm-corrected pair detection rate | 0.10 |
| direction_miss | One minus the mean correct significant direction rate, with full-grid Holm correction | 0.05 |
| direction_wrong | Mean wrong significant direction rate, with full-grid Holm correction | 0.05 |

The four interactions receive equal weight when summarizing pair detection or direction rates. A configuration must have at least 80 usable repetitions and at least 95% usable results for point estimation and interpretation to enter the recommendation. The complete score components and availability are retained. Coverage, SE calibration, variable-identification rates and computational cost accompany the recommendation for scientific review.

Review A's tables and figures, choose `freeze$config_id`, then run the freeze action. B, C and D inherit the selected N,L,M,kappa,lambda. Their `common_fixed`, `B_overrides`, `C_overrides` and `D_overrides` entries allow explicit alternatives, which are saved with each run.

## Interpretation and reporting settings

`advanced_overrides$interpretation` controls contribution evaluation points (default 1,000), chunk size (128), single-exposure grid points (41), interaction grid points per axis (9), interaction alpha (0.05), numerical tail accuracy (1e-8), integration term limits (100,000/1,000,000), and joint OR quantiles (0.25,0.75,0.90). Exact reference concentrations are added to the grids. These settings define evaluation targets and are fixed before fitting.

`report_overrides` controls font, resolution, figure width, point-table output, scatter display density and C figure pagination. Numerical model outputs remain available for independent reporting. Parameters that define scientific evaluation grids are recorded with the fits.

Kernel centers use k-means++ coverage selection with 25% coverage anchors; distance parameters are learned from the independent outcome sample within [0.01,100] times the geometry-based inverse-distance scale. The initial distance is 0.0625 times that same scale. Each matching group stays together in three-fold nuisance cross-fitting. The source configuration lists optimization tolerances and iteration limits.

## Module D settings

`repetitions$D` defaults to 100. `D_correlation_levels` specifies `correlation_id`, `correlation_label`, `target_mean`, `min_abs` and `mean_tolerance`. Moderate/high target mean absolute Pearson correlations are 0.60/0.85, the minimum absolute pairwise values are 0.45/0.75, and both mean tolerances are 0.03. Internal score-mixing parameters are calibrated independently and stored. IDs must be unique safe names; D0 is reserved for the completed B reference. `D_baseline_run=NULL` discovers B in the work directory. Model settings are checked against that B run. See MODULE_D.md for the marginal maps, calibration and evaluation definitions.

The reviewed work-directory configuration supplies D's five model parameters. `D_overrides` is empty by default. All entries in that list are recorded explicitly as run-specific choices. Parameter settings remain constant across D's correlation levels.

## Module E settings

`repetitions$E=100` and `E_settings` define seven default exposure settings, a 30,000-set configuration pool, three validation folds, six M–kappa pairs, and lambda values 0.3, 1, 3. The frozen configuration is always a candidate. `E_overrides` can change the fixed baseline and the common N/L budgets. Default reviewed budgets are N=20,000 and L=20,000. Each fold trains on 20,000 configuration-pool sets and scores 10,000 held-out sets. The independent learning pool is shared across folds and strategies. An explicit `manual_candidates` data frame with M, kappa and lambda replaces the Cartesian candidate list. `kernel_candidates` replaces the complete M–kappa table. All choices and effective parameters are recorded. See MODULE_E.md for the parameter table and data roles.
