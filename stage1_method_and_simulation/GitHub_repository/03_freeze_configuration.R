# Convenience entry: settings are read from 05_submit.R without editing this file.
.entry_root <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  path <- if(length(paths)) tail(paths,1)[[1]] else if(length(arg)) sub("^--file=","",arg[1]) else "03_freeze_configuration.R"
  normalizePath(dirname(path), winslash="/", mustWork=TRUE)
})
.entry_env <- new.env(parent=.GlobalEnv)
.entry_env$.submission_action <- "freeze"
source(file.path(.entry_root,"05_submit.R"),local=.entry_env,encoding="UTF-8")
