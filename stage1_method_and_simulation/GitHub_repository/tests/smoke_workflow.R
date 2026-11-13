run_smoke_workflow <- function(root, output_root=tempdir()) {
  # This creates its own small project and never writes into formal results/.
  dest<-file.path(output_root,"validation","workflow_check")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(d in c("R","config","monitor")) {
    dir.create(file.path(dest,d),showWarnings=FALSE)
    file.copy(list.files(file.path(root,d),full.names=TRUE),file.path(dest,d),overwrite=TRUE,recursive=TRUE)
  }
  file.copy(file.path(root,c("load.R","00_parameters.R")),dest,overwrite=TRUE)
  source(file.path(root,"00_parameters.R"),encoding="UTF-8")
  cfg<-simulation_config(dest);cfg$schema<-"implementation-check-v1";cfg$seed<-947121L
  cfg$execution$workers<-2L;cfg$execution$block_groups<-200L;cfg$execution$progress_seconds<-15
  cfg$evaluation$groups<-250L;cfg$evaluation$function_points<-60L;cfg$evaluation$coverage_points<-60L
  cfg$evaluation$quantiles<-c(.1,.25,.5,.75,.9);cfg$selection$evaluation_points<-60L
  cfg$interpretation$contribution_points<-60L;cfg$interpretation$single_points<-7L;cfg$interpretation$interaction_points<-4L
  cfg$kernel$penalty_rows<-50L;cfg$report$dpi<-100L
  cfg$A$manual_design<-data.frame(N=c(400L,600L),L=600L,M=12L,kappa=.001,lambda=3)
  cfg$A$repetitions<-2L;cfg$A$recommendation$min_repetitions<-2L
  cfg$B$repetitions<-2L;cfg$C$repetitions<-2L;cfg$C$null_counts<-c(4L,8L)
  a<-run_simulation(cfg,dest,"A")[[1]];prepare_results(a);freeze_configuration(dest,a,"A001")
  b<-run_simulation(cfg,dest,c("B","C"));lapply(b,prepare_results)
  all<-lapply(c(a,unlist(b)),sim_collect)
  check<-kdr_bind(lapply(all,function(x)x$status))
  kdr_csv(check,file.path(output_root,"validation","workflow_task_status.csv"))
  stopifnot(all(check$point_ok),all(check$interval_ok),
    all(check$selection_ok[check$method=="KDR"&check$scenario_gf=="G1"]),
    sum(check$reused)==0L)
  for(run in c(a,unlist(b))) {
    rr<-lapply(list.files(file.path(run,"fits"),full.names=TRUE),readRDS)
    stopifnot(all(vapply(rr,function(x)isTRUE(x$interpretation_ok),logical(1))),
      all(vapply(rr,function(x)isTRUE(x$interpretation$detection$ok),logical(1))))
  }
  # Resume must preserve fit files and must not trigger new computation.
  files<-list.files(a,"[.]rds$",full.names=TRUE,recursive=TRUE)
  fits<-files[grepl("/fits/",gsub("\\\\","/",files))];before<-tools::md5sum(fits)
  run_simulation(cfg,dest,"A")
  stopifnot(identical(before,tools::md5sum(fits)))
  kdr_log("Parallel end-to-end, interpretation and resume checks passed.")
  invisible(dest)
}
