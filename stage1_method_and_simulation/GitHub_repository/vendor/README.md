# Offline dependency: CompQuadForm 1.4.4

This directory contains the **unmodified CRAN source archive**, not a Windows or macOS binary. It is compiled for the R interpreter and OS on which the install action is run.

- Package: CompQuadForm, version 1.4.4, published 2025-07-13.
- Authors: Pierre Lafaye De Micheaux, Pierre Duchesne, Robert B. Davies, R. W. Farebrother, J. Sheil and I. O Muircheartaigh, as listed in the original package DESCRIPTION.
- License: GPL (>= 2), as specified in the original archive.
- Project page: https://cran.r-project.org/package=CompQuadForm
- Download source: https://cran.r-project.org/src/contrib/CompQuadForm_1.4.4.tar.gz
- Archive size: 19,797 bytes.
- MD5: `44e119935fa15aa034de317d41e801d8`
- SHA256: `84af2087bf7cd6dfb7aa9468fd88939267559f46c0dc7bea9d010314f42fa143`

The original DESCRIPTION contains no Depends, Imports or LinkingTo declarations; it declares `NeedsCompilation: yes`. It does not declare an R 4.x minimum. Actual installation and the project's compatibility check must still succeed on the server's R 3.6.3 build.

Run `Rscript --vanilla 05_submit.R install`. The installer verifies the archive checksum, installs into the R user library with `repos=NULL`, and requires no repository connection. An already available package is retained and its version is reported.

The server needs make, C/C++ compilation tools and the development headers for its R installation. No administrator access is needed to install into an already writable user library; missing system compilation tools require the administrator's assistance. Do not copy the existing Windows DLL package to Linux.
