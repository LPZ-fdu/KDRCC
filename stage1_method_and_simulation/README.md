# Stage 1 Method development and simulation

This stage contains the KDR-CC method specification, complete Chinese theory appendix, simulation implementation, reviewed parameter configuration and completed A–E assessment records. Its frozen study identifier is stage1-20261107.

## Repository contents

`GitHub_repository` contains the executable study source, configuration, documentation, appendix, numerical checks, archive tools and publication reports. Start with its [README](GitHub_repository/README.md). `records/results/module_*/analysis/index.html` contains the figures and tables for each completed module. Submission and fixed-configuration records are in `records/submissions` and `records/config`.

The modules evaluate parameter configuration and sample size (A), auxiliary-model robustness and three-method comparisons (B), additional null exposures (C), stronger exposure correlation (D), and independent data-driven configuration selection (E). Model interpretation includes joint and single-exposure contrasts, observed-combination predictions, interval performance, variable participation, contribution shares and interaction detection.

The [variable identification supplement](variable_identification_supplement/README.md) compares KDR-CC with spline and linear conditional logistic regression across B, C, and D. Its [publication report](variable_identification_supplement/results/index.html) includes overall and exposure-specific identification rates, null-exposure inclusion rates, and paired method differences. Detailed supplementary records have a separate release asset manifest.

## Complete runtime assets

Seven independent `stage1_runtime_*.zip` assets accompany this repository through GitHub Releases. Download all seven and place them in `release_assets` beside `GitHub_repository`, then follow [ARCHIVE_LAYOUT.md](GitHub_repository/docs/ARCHIVE_LAYOUT.md). Source and runtime manifests preserve sizes and SHA-256 hashes.

```bash
cd GitHub_repository
python tools/verify_stage1.py
python tools/verify_stage1.py --assets ../release_assets
python tools/restore_records.py --assets ../release_assets --destination ..
```

Restoration supplies complete fitted records, learned kernels, evaluation objects, selection records and execution logs for reproducibility. These runtime assets supplement the repository's core source and publication reports.

## Relationship to Stage 2

Stage 1 provides the method implementation and statistical simulation evidence. Stage 2 distributes the public R package, user documentation and an independent matched reproduction of the formal KDR-CC G1 study. Both stages are parts of the same project and retain separate execution records.
