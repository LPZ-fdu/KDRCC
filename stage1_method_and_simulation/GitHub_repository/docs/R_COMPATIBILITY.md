# R compatibility

The implementation targets R 3.6.3 and newer. Code uses base R, splines, parallel PSOCK workers and CompQuadForm. On R 3.6 the loader sets stringsAsFactors=FALSE to match character-column handling across supported versions. Saved objects use serialization version 3. Text files are UTF-8; source files are loaded with explicit encoding.

The install action uses the bundled unmodified CompQuadForm CRAN source archive for offline Linux installation when available. The vendor directory includes its provenance. Standard recommended R packages supply parallel and splines.

The check action parses code on the actual interpreter, verifies dependencies, validates the weighted-chi-square calculation, tests serialization and launches two independent workers to check package versions and seed parity. Unit tests verify matching, analytic derivatives, objective gradients, response fixtures and contribution identities. Run this check with the same Rscript executable used for the study.

Cairo supplies PNG output without an interactive display. The requested font and its effective Linux fallback are recorded in each report. Memory measurements depend on available platform samplers; their availability is explicit in resource tables.

Module D uses base R empirical quantile maps, Gaussian transforms and scalar root solving and the same PSOCK, estimation and reporting interfaces. Its implementation targets R 3.6.3 and newer. The `check D` action runs the numerical checks on the actual server interpreter before the plan is executed.

Module E uses base R and `splines::ns` for its training-data nuisance basis, with the existing CompQuadForm inference and PSOCK infrastructure. Its syntax targets R 3.6.3 and newer. `check E` verifies the installed server interpreter and dependencies before fitting.
