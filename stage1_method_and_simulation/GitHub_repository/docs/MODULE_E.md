# Module E: data-driven configuration selection

## Objectives

Module E evaluates a configuration-selection procedure that uses observed matched-set outcomes. It quantifies whether this procedure improves effect estimation, variable identification and interpretation relative to the frozen configuration when exposure dimension or dependence increases. All scientific evaluations use new independent repetitions. The fixed and adaptive strategies are compared on identical final-estimation data and evaluation combinations within a repetition.

## Settings and computation

| Setting | Active exposures | Null exposures | Dependence |
|---|---:|---:|---|
| E0 | 8 | 0 | Original |
| EQ8 | 8 | 8 | Original |
| EQ16 | 8 | 16 | Original |
| EQ40 | 8 | 40 | Original |
| EQ80 | 8 | 80 | Original |
| EDM | 8 | 0 | Moderate |
| EDH | 8 | 0 | High |

Each setting has 100 independent repetitions. The moderate/high copula specification and calibration use the D correlation definitions. The conditional response surface and coefficients remain those in `config/mechanism.R`. Overall dependence calibration is saved with E and is reproducible independently of stored D fits.

The fixed comparator reads `<work_root>/config/fixed_parameters.R`, with optional `common_fixed` and `E_overrides`. The reviewed configuration is N=20,000, L=20,000, M=144, kappa=0.0001, lambda=1. The effective values are printed by the plan action. The candidate kernel configurations are:

| M | kappa |
|---:|---:|
| 96 | 0.0001 |
| 144 | 0.00003 |
| 144 | 0.0001 |
| 144 | 0.0003 |
| 192 | 0.0001 |
| 288 | 0.0001 |

Each kernel configuration is crossed with lambda=0.3, 1 and 3, producing 18 candidates. The fixed configuration is always included; an explicit override outside this set adds a candidate. Kernel centers, distance learning, nuisance fits and features are shared where their inputs and definitions agree. Each candidate is evaluated in three folds. The default study contains 700 paired repetitions, 4,200 learned kernels, 37,800 validation fits and 1,400 final fits. Final fits additionally evaluate uncertainty, contributions, variable tests and interactions. Workloads and elapsed times therefore differ between validation and final fits.

## Independent data roles

Each repetition generates three mutually independent data pools:

1. **Kernel-learning pool, L complete matched sets.** Learn exposure scaling, centers, distance parameters, whitening, the penalty operator and the nuisance dictionary. This pool is part of the training information shared by every validation fold. It is independent of every held-out outcome.
2. **Configuration pool, 30,000 complete matched sets.** Assign sets to three balanced folds. In each round, fit the full KDR-CC estimator on two folds (20,000 sets by default), including its internal nuisance cross-fitting; score the remaining 10,000 sets. The kernel-learning pool is unchanged across folds. Reusing its fitted objects is exact computational reuse for this partitioned learning scheme.
3. **Final-estimation pool, N complete matched sets.** After selection, estimate the selected model and the fixed comparator on the same fresh data. Both use their corresponding independently learned kernels. N defaults to the frozen value 20,000. The complete default repetition uses 70,000 sets across the three fitting roles, plus the shared independent evaluation distribution.

Matched sets are never split across folds. Scenario, repetition, role and block determine semantic random seeds. Scheduling and worker count do not enter these seeds. Generated true responses and outcome probabilities are removed from the data supplied to fitting and validation.

## Data-derived nuisance specification

E uses the full W and W-squared adjustment and the nine ordered-pair confounder cells for f. The exposure component used to estimate the confounder coefficients comprises natural cubic splines of every candidate exposure and a low-rank kernel basis. Interior spline knots are the learning-sample 25th, 50th and 75th percentiles; boundary knots are its 0.5th and 99.5th percentiles. Up to 24 retained kernel basis columns supplement the additive dictionary. All transformations and rank reduction are fixed from the independent learning pool. The same prescription is used for both strategies, including null exposures.

The internal G1 label denotes full adjustment and the full f cell specification. The exposure nuisance regression is a data-derived approximation. It is fitted from observed data; the simulation truth is used to define the evaluation targets and score their recovery. A–D retain their established nuisance specifications. The paired fixed arm in E provides the direct comparator for the application-oriented procedure.

## Validation criterion and selection rule

For candidate c and held-out matched set i, define eta_ij = h_hat_c(A_ij) + gamma_hat_c1 W_ij + gamma_hat_c2 W_ij^2. The probability of the observed case day is its exp(eta) divided by the sum over its matched set. Computation uses a log-sum-exp transformation.

The score is minus the sum of log case-day probabilities across all held-out sets, divided by their total count. It is a held-out conditional negative log likelihood. The candidate's final regularized KDR coefficients enter this score; the kernel-learning objective alone does not determine selection. Validation loss contains no penalty term. The three folds are weighted by their actual matched-set counts. A candidate is eligible only when all folds return finite scores. The minimum score determines the selected configuration; exact ties use the recorded candidate order. If every candidate fails, the predeclared rule attempts the frozen configuration and records `fallback=TRUE` and `selection_ok=FALSE`.

The score targets prediction on the scenario's observed exposure distribution. Joint-path accuracy, interaction detection, contribution allocation and variable identification are evaluated separately after selection. This distinguishes an improved validation score from improvement in each scientific target.

## Final evaluation and aggregation

Every repetition yields two final records, labelled Fixed and Adaptive. The adaptive record stores its actual M, kappa, lambda and kernel key. Scientific summaries aggregate all selected configurations within the Adaptive strategy, so they evaluate the entire selection procedure.

For joint and observed-combination effects, report RMS bias, RMSE, absolute relative OR errors, empirical SD, model SE, pointwise coverage, interval width and correct/wrong significant directions. Fixed evaluation points within a scenario permit pointwise bias and empirical variability across repetitions. The joint coordinates of the eight active exposures are common across dependence settings; actual-combination points follow their own scenario distributions. Null variables are included in model inputs and joint coordinates. Contribution evaluation varies the eight active coordinates while fixing all additional coordinates at their reference values.

Variable identification reports each active-exposure rate and the mean null inclusion rate. Contribution shares normalize the mean absolute baseline-Shapley allocation within the eight prespecified active exposures at the reference background. Interaction summaries preserve the four-pair Holm test and the gridwise corrected direction criterion used in B. Quantitative OR contrasts and numerical convergence diagnostics are retained. Covariance calculations condition on the independent learned objects and selected configuration; the simulation summaries include their repetition-to-repetition variation.

Paired changes are Adaptive minus Fixed for the same scenario, repetition and target. Their Monte Carlo uncertainty uses complete pairs. Marginal performance tables retain each strategy's availability counts. Missing fits and tests have explicit failure records. Configuration selection frequencies, fallback rates, all fold scores, candidate rankings and selected values are saved for audit and reproduction.

The support table describes squared Mahalanobis distances of the joint path in the eight active coordinates relative to an independent scenario sample. Its empirical 95th percentile supplies a geometric reference. Errors at sparse combinations remain visible in the full joint-path evaluation.

## Submission and output

Edit only `05_submit.R`. `E_settings` controls the scenarios, configuration-pool size, validation folds and candidates; `repetitions$E` controls repetition count; `E_overrides` controls the fixed comparator and N/L budget. `manual_candidates` accepts a data frame with M, kappa and lambda columns. `workers` controls independent PSOCK processes, default 32. E settings also accept the nuisance dictionary controls documented above.

From the source directory:

```bash
nohup Rscript --vanilla 05_submit.R check E > E_check.txt 2>&1 &
nohup Rscript --vanilla 05_submit.R plan E > E_plan.txt 2>&1 &
```

Review each log after its process completes. Then run:

```bash
nohup Rscript --vanilla 05_submit.R run E > E_run.txt 2>&1 &
tail -f E_run.txt
```

After fitting, report separately:

```bash
nohup Rscript --vanilla 05_submit.R report E > E_report.txt 2>&1 &
```

The work directory remains the sibling `KDRCC_formal_20261027` by default. E results are saved under `results/module_E_<identifier>/`. `analysis/index.html` links the English PNG figures, formatted tables and captions. `analysis/model_diagnostics/` contains numerical diagnostics. `choices/` stores selected configurations; `cross_validation/` stores fold assignments, per-fold scores and candidate rankings. Checkpoints are resumed with the same settings. Set `clear_stale_lock=TRUE` only after confirming all previous processes for that run have stopped. Reporting can be repeated from a complete copied E run.

## Methodological reference

Cawley GC, Talbot NLC. On over-fitting in model selection and subsequent selection bias in performance evaluation. Journal of Machine Learning Research. 2010;11:2079–2107. https://www.jmlr.org/papers/v11/cawley10a.html
