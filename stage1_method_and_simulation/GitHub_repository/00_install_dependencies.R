# Only one non-recommended CRAN package is required. R >= 3.6.3.
if(getRversion()<"3.6.3")stop("R 3.6.3 or newer is required.")
if(!requireNamespace("CompQuadForm",quietly=TRUE)) {
  user_library <- path.expand(Sys.getenv("R_LIBS_USER"))
  if(!nzchar(user_library))stop("Set R_LIBS_USER to a writable package library.")
  dir.create(user_library,recursive=TRUE,showWarnings=FALSE)
  .libPaths(c(user_library,.libPaths()))
  install.packages("CompQuadForm",lib=user_library,repos="https://cloud.r-project.org")
}
if(!requireNamespace("CompQuadForm",quietly=TRUE))stop("CompQuadForm installation did not complete.")
cat("Dependency available: CompQuadForm",as.character(packageVersion("CompQuadForm")),"\n")
cat("R",as.character(getRversion()),"|",R.version$platform,"\n")
cat("Use the same R and package versions across machines when reproducing results.\n")
