# kdrcc package distribution

Author and maintainer: Peizheng Li, pzli25@m.fudan.edu.cn. Version 0.1.0. GPL-3. R >= 3.6.3.

## User documentation

- [Chinese user guide](kdrcc_user_guide_zh.docx)
- [English user guide](kdrcc_user_guide_en.docx)
- [Package source and function help](kdrcc/README.md)

The guides cover installation, matched data, sample allocation, fitting, every public function and parameter, diagnostics, OR prediction, joint and single-exposure curves, interaction and variable participation tests, contribution allocation, tuning and custom reporting.

## Installation

```r
source("INSTALL_AND_TEST.R", encoding="UTF-8")
```

The installer selects the appropriate archive, installs or reuses CompQuadForm, and writes installation workflow results to `installation_test`. Configure the library, archive, repository and output directory in the settings at the top of the script. The distribution includes source archives for supported platforms and Windows R 4.4 binary archives.

## Examples

`examples/example_matched_data.csv` contains 800 matched sets with one event and three reference dates each. `examples/RUN_EXAMPLE.R` fits a complete analysis and saves CSV tables, PNG figures and RDS objects under `examples/output`. `examples/RUN_TUNING.R` demonstrates independent learning, configuration selection and final estimation. `examples/CUSTOM_PLOTS.R` generates ggplot2 figures from the example CSV outputs. `examples/README.md` lists files and output definitions.

## Source repository

`kdrcc` contains DESCRIPTION, NAMESPACE, R code, help pages, method documentation, GPL-3 license, citation and package tests. `installers` contains its distributable archives and the bundled CompQuadForm dependency. Source R files define the statistical implementation; the guides and examples use these public interfaces.
