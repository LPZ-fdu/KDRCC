run_module_E_smoke <- function(root,output_root=tempdir()) {
  source(file.path(root,"00_parameters.R"),encoding="UTF-8",local=TRUE)
  cfg<-simulation_config(root);cfg$schema<-"module-E-workflow-check-v1";cfg$seed<-772091L
  cfg$modules<-"E";cfg$E$repetitions<-2L;cfg$E$null_counts<-2L;cfg$E$dependence<-"Original"
  cfg$configuration$overrides$E<-list(N=600L,L=600L,M=8L,kappa=.001,lambda=1)
  cfg$E$tuning_sets<-900L;cfg$E$validation_folds<-3L
  cfg$E$manual_candidates<-data.frame(M=c(8L,8L,10L),kappa=.001,lambda=c(1,3,1))
  cfg$E$nuisance_kernel_rank<-4L
  cfg$execution$workers<-2L;cfg$execution$block_groups<-200L
  cfg$execution$progress_seconds<-15;cfg$resources$enabled<-FALSE
  cfg$evaluation$groups<-300L;cfg$evaluation$function_points<-20L;cfg$evaluation$coverage_points<-20L
  cfg$evaluation$quantiles<-c(.1,.25,.5,.75,.9);cfg$selection$evaluation_points<-20L
  cfg$interpretation$contribution_points<-20L;cfg$interpretation$single_points<-5L
  cfg$interpretation$interaction_points<-3L;cfg$kernel$penalty_rows<-25L
  cfg$report$dpi<-72L;cfg$report$font<-"serif"
  work<-file.path(output_root,"module_E_smoke");dir.create(work,recursive=TRUE,showWarnings=FALSE)
  run<-run_simulation(cfg,work,"E")[[1]]
  z<-sim_collect(run);kdr_csv(z$status,file.path(work,"workflow_status.csv"))
  stopifnot(nrow(z$status)==8L,all(z$status$point_ok),all(z$status$interval_ok),all(z$status$interpretation_ok))
  choices<-lapply(list.files(file.path(run,"choices"),pattern="[.]rds$",full.names=TRUE),readRDS)
  stopifnot(length(choices)==4L,all(vapply(choices,function(a)a$selection_ok,logical(1))))
  for(a in choices)stopifnot(sum(a$ranking$eligible)==3L)
  # The reporting groups span selected configurations; they remain strategy-level groups.
  stopifnot(nrow(z$effects)==8L)
  fits<-list.files(file.path(run,"fits"),full.names=TRUE);before<-tools::md5sum(fits)
  resumed<-run_simulation(cfg,work,"E")[[1]]
  stopifnot(identical(run,resumed),identical(before,tools::md5sum(fits)))
  prepare_results(run)
  captions<-read.csv(file.path(run,"analysis","Figure_and_table_captions.csv"),stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(captions$id),all(file.exists(file.path(run,"analysis",captions$file))))
  kdr_log("Module E paired fits, parallel execution, full reporting and resume checks passed.")
  invisible(run)
}
