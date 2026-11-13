# Evaluate contributions from saved models using settings in 05_submit.R.
# Usage: Rscript --vanilla 06_prepare_contributions.R [C,E]
.entry_root <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  path <- if(length(paths)) tail(paths,1)[[1]] else if(length(arg)) sub("^--file=","",arg[1]) else "06_prepare_contributions.R"
  normalizePath(dirname(path), winslash="/", mustWork=TRUE)
})
.entry_env <- new.env(parent=.GlobalEnv)
.entry_env$.submission_action <- "contributions"
.entry_args <- commandArgs(trailingOnly=TRUE)
if(length(.entry_args)>1L)stop("Usage: Rscript 06_prepare_contributions.R [C,E]")
.entry_env$.submission_modules <- if(length(.entry_args))strsplit(.entry_args[1],",",fixed=TRUE)[[1]] else c("C","E")
source(file.path(.entry_root,"05_submit.R"),local=.entry_env,encoding="UTF-8")
