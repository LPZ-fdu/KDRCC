prepare_stage_one_G1_results <- function() {
  script<-NULL
  for(fr in sys.frames())if(!is.null(fr$ofile))script<-fr$ofile
  args<-grep("^--file=",commandArgs(FALSE),value=TRUE)
  if(is.null(script)&&length(args))script<-sub("^--file=","",args[1])
  root<-if(is.null(script))getwd() else dirname(normalizePath(script,winslash="/",mustWork=TRUE))
  env<-new.env(parent=globalenv());sys.source(file.path(root,"END_TO_END_SETTINGS.R"),env)
  sys.source(file.path(root,"scripts/04_reproduce_G1.R"),env)
  sys.source(file.path(root,"scripts/05_prepare_G1_results.R"),env)
  run<-env$rp_find_run(root,env$reproduction_settings)
  env$rp_report(root,run,env$reproduction_settings)
}
prepare_stage_one_G1_results()
