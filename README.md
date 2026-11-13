# KDR-CC Method development and R package

Kernel-Based Doubly Robust Case-Crossover Regression (KDR-CC) models nonlinear multiexposure associations in individually matched case-crossover data. The project combines the statistical method and simulation study with a public R package and independently executed package assessment.

Maintainer: Peizheng Li, pzli25@m.fudan.edu.cn. Package version: 0.1.0. Package license: GPL-3. Supported R versions: 3.6.3 and later.

## Project stages

| Stage | Directory | Contents |
|---|---|---|
| Method development and simulation | [stage1_method_and_simulation](stage1_method_and_simulation/README.md) | Theory appendix, method specification, study code, A–E designs, fixed configuration, reports, actual execution records and runtime-archive manifests |
| R package and independent assessment | [stage2_r_package](stage2_r_package/README.md) | Public package source and installers, Chinese and English user guides, case data and examples, independent checks, formal matched G1 reproduction and final outputs |

The method-development record and package-development record are complementary project components. Stage 1 evaluates statistical performance; Stage 2 documents public software use and verifies numerical agreement with the formal study.

## Use the package

```r
source("stage2_r_package/package/INSTALL_AND_TEST.R", encoding="UTF-8")
source("stage2_r_package/package/examples/RUN_EXAMPLE.R", encoding="UTF-8")
```

Read the [English guide](stage2_r_package/package/kdrcc_user_guide_en.docx) or [Chinese guide](stage2_r_package/package/kdrcc_user_guide_zh.docx) for complete use. The source package is [stage2_r_package/package/kdrcc](stage2_r_package/package/kdrcc). The guides cover every exported function, configuration, returned object and the connection to user-defined tables and ggplot2 graphics.

## Review and reproduce results

- Stage 1: open `stage1_method_and_simulation/GitHub_repository/records/results/module_*/analysis/index.html` and read the study source README for submission and restoration.
- Stage 1 supplement: see [variable identification in matched simulations](stage1_method_and_simulation/variable_identification_supplement/README.md) for the three-method comparison and its publication report and runtime records.
- Stage 2: open `stage2_r_package/validation/results/G1_reproduction_0a040df3a874/analysis/index.html`. The completed Linux R 3.6.3 run includes 100 matched formal repetitions and 2,900 numerical comparisons, all passing.
- To run the independent package assessment, execute `RUN_VALIDATION.R` from `stage2_r_package/validation`. Plot regeneration is a separate entry, `PREPARE_RESULTS.R`.

The seven complete Stage 1 runtime ZIP assets are distributed through the accompanying GitHub Release. Archive layout and SHA-256 manifests are in Stage 1. Downloading and restoring those assets supplies the full saved models and execution material used by the statistical study.

## Integrity and citation

`SUBMISSION_CHECKSUMS.csv` lists repository files and their SHA-256 values. Each stage also has source or validation manifests. The package GPL-3 license and CITATION are supplied in its source. Use `citation("kdrcc")` for the software citation. Study publication metadata can be maintained in this README together with the manuscript citation.
