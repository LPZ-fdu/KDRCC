# Final methods and implementation review

Review date: 2026-11-07

The review covered the appendix definitions and proofs, their implemented special cases, the formal A–E plans and fitted records, and the source snapshots accompanying each run. The review identified a coherent mathematical specification and implementation for the declared KDR-CC procedure. The method is a regularized finite-kernel estimator built from a doubly robust base moment.

## Statistical conclusions

The symmetric ordered-pair construction preserves conditional direction odds and original matched-set clustering. The base moment has the stated product-bias identity under a correct finite response representation. The fixed-penalty criterion has its own population target; the local target-displacement result follows by the implicit function theorem. Conditional covariance uses the exact criterion Hessian and both auxiliary-estimation paths. Local consistency and normality are stated for fixed dimensions, regular auxiliary probability limits, positive local curvature and vanishing numerical error at the required rate.

Derivative and interaction tests use centered Gaussian quadratic-form references under their stated working-response nulls. The finite positive-distance Gaussian space gives the coordinate-zero null its stated structural meaning. Holm control uses valid working-null tail probabilities. Simulation recognition rates evaluate generating active/null labels. Contribution summaries allocate the specified active-coordinate contrast and normalize mean absolute allocations within that set.

The appendix therefore supplies explicit identification, target, estimation, computation and conditional-inference statements for the implemented analysis. Theoretical scope and empirical performance are specified separately for effect estimation, variable participation and contribution allocation. The code retains empirical bias, coverage and recognition outcomes across the tested dimension and correlation scenarios.

## Implementation evidence

All 59 R source files parsed successfully. Algebraic and numerical checks passed: 30 numerical units, 20 method-contract checks, 22 D checks, E design/grouped-scoring checks, 10 active-contribution checks and 60 source/archive checks. The archive checks covered all five modules. The nine core numerical, kernel, estimation, covariance, derivative and interaction files were byte-identical to the corresponding formal run snapshots. Conditional generator and nuisance-design fixtures agreed across snapshots. Complete file-level comparison is recorded in `archive/snapshot_comparison.csv`.

The review used deterministic algebraic fixtures and saved records. All formal fitted records were read successfully. Their counts are:

| Module | Planned | Saved | Fit flag successful | Eight-variable contribution records |
|---|---:|---:|---:|---:|
| A | 5000 | 5000 | 4998 | 4998 |
| B | 800 | 800 | 800 | 800 |
| C | 400 | 400 | 399 | 399 |
| D | 1600 | 1600 | 1600 | 1600 |
| E | 1400 | 1400 | 1396 | 1396 |

The seven unsuccessful fit records remain in the runtime archive with their errors and availability flags. The definition-to-code matrix is `archive/theory_code_map.csv`. The latest project appendix is preserved byte for byte in `manuscript/Appendix_KDRCC_zh.docx`; its mathematical content agrees with the structured document source. Formula numbering and equation objects were checked structurally.

## Frozen scientific settings

The generating mechanism identifier is `9d3750065660d0322c16b0fd90ff134c`. The reviewed fixed configuration is A035: N=20,000; L=20,000; M=144; kappa=0.0001; lambda=1. Module E additionally retains every validation candidate score, selected setting and paired final fit. Configuration and reference-function definitions accompany the results.

The present review supports freezing the method and A–E development-stage evaluation archive under these definitions. Package-level behavior and reproducibility can be assessed against the saved configurations, evaluation inputs, fitted functions and covariance objects.
