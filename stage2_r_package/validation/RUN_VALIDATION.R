# Numerical validation and matched formal G1 reproduction.
run_kdrcc_validation <- function(include_end_to_end=TRUE) {
  script<-NULL
  for(fr in sys.frames())if(!is.null(fr$ofile))script<-fr$ofile
  args<-grep("^--file=",commandArgs(FALSE),value=TRUE)
  if(is.null(script)&&length(args))script<-sub("^--file=","",args[1])
  root<-if(is.null(script))getwd() else dirname(normalizePath(script,winslash="/",mustWork=TRUE))
  original<-getwd();on.exit(setwd(original),add=TRUE);setwd(root)
  source(file.path(root,"scripts/01_validate.R"),local=new.env(parent=globalenv()),encoding="UTF-8")
  source(file.path(root,"scripts/02_figures.R"),local=new.env(parent=globalenv()),encoding="UTF-8")
  if(isTRUE(include_end_to_end))source(file.path(root,"RUN_END_TO_END.R"),local=new.env(parent=globalenv()),encoding="UTF-8")
  cat("Validation results:",file.path(root,"results"),"\n")
}
run_kdrcc_validation()
