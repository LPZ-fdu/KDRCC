# Implementation checks

`Rscript --vanilla 05_submit.R check` runs interpreter/dependency checks and numerical unit checks in the writable work directory.

The numerical checks cover calendar matching, one event per set, pair orientation and normalization, deterministic streams, dimension-preserving null augmentation, fixed-response fixtures, Gaussian derivatives, penalized objective gradients and numerical quadratic-form tails. Analytic true Shapley allocations and Gaussian quadrature allocations are compared with independent subset enumeration. Contrast SEs are checked against the full covariance expression.

`tests/smoke_workflow.R` defines `run_smoke_workflow(source_root, output_root)`, an explicit small parallel A→freeze→B/C workflow. It checks result availability, reporting and resume behavior. The smoke configuration uses two repetitions and small samples to exercise the implementation. Formal settings are obtained from `00_parameters.R` through the submission script.

`tests/response_reporting_checks.R` checks analytic interaction contrasts, additive-comparator identities, saved predictions, independent-repetition summaries, model reuse and figure outputs for a completed B run. `tests/server_workflow_checks.R` checks source/work separation and submission orchestration. Validation records are written outside the source project when a separate output root is supplied.

## Module D checks

The submission `check` action also evaluates D's 1,600-new-fit/200-kernel plan, scenario keys, Gaussian-score transforms, marginal variances, signed dependence, identity behavior, matching, one event per set and independent seed namespace. `run_module_D_smoke()` is an explicit small two-worker workflow test. It creates a small B reference, imports the saved fits, runs moderate/high D scenarios, checks common targets and resumption, verifies that original B files remain unchanged, and validates the combined report. Its output directory is separate from formal results.

## Module E checks

`tests/module_E_checks.R` validates candidate construction, complete-fold eligibility, fallback behavior, grouped subsetting, truth-independent nuisance features, conditional likelihood scoring and independent random namespaces. `tests/module_E_smoke.R` runs two settings and two independent repetitions with two workers, three candidates and three validation folds: 36 validation fits and eight paired final fits. It checks covariance and interpretation outputs, aggregation over selected configurations, unchanged fit checkpoints on resume, and the full report. All validation output is written to the supplied separate work directory. The submission check action includes E's inexpensive unit checks.

## Statistical contract checks

`tests/method_contract_checks.R` verifies the exact double-robust product identity, both nuisance-correct branches, the penalized exact Hessian, all-fold g/f gradient sensitivities, the effective penalty scaling, the local penalty-target displacement, and the positive-semidefinite trace-normalized mixed-derivative matrix. The submission `check` action includes these deterministic checks. Numerical fixtures test the actual objective and covariance derivatives at nonzero moment residuals.
