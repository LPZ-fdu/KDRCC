# Explicit opt-in implementation check. No formal defaults are changed.
run_server_workflow_checks <- function(source_root,work_root) {
  source_root<-normalizePath(source_root,winslash="/",mustWork=TRUE)
  source(file.path(source_root,"load.R"),encoding="UTF-8");load_kdrcc_simulation(source_root)
  e<-new.env(parent=.GlobalEnv)
  expr<-parse(file.path(source_root,"05_submit.R"),encoding="UTF-8")
  # The first three expressions contain only editable paths and settings.
  for(i in 1:3)eval(expr[[i]],e)
  s<-e$submission;s$action<-"check";s$workers<-2L
  s$seed<-984719L;s$repetitions<-list(A=2L,B=2L,C=2L);s$null_counts<-c(0L,4L,8L)
  s$A_overrides<-list(manual_design=data.frame(N=400L,L=600L,M=12L,kappa=.001,lambda=3),
    recommendation=list(min_repetitions=2L))
  s$advanced_overrides<-list(schema="server-implementation-check-v1",
    execution=list(block_groups=200L,progress_seconds=15),resources=list(enabled=FALSE),
    evaluation=list(groups=250L,function_points=60L,coverage_points=60L,quantiles=c(.1,.25,.5,.75,.9)),
    interpretation=list(contribution_points=60L,single_points=7L,interaction_points=4L),
    selection=list(evaluation_points=60L),kernel=list(penalty_rows=50L))
  s$report_overrides<-list(dpi=100L,width_inches=12,base_cex=1.12,
    scatter_points=60L,save_pointwise_tables=TRUE,C_panels_per_page=2L)
  files<-list.files(source_root,recursive=TRUE,full.names=TRUE,all.files=TRUE)
  files<-files[!file.info(files)$isdir];before<-tools::md5sum(files)
  execute_submission(s,source_root,work_root)
  s$action<-"run";s$modules<-"A";a<-execute_submission(s,source_root,work_root)[[1]]
  s$action<-"report";execute_submission(s,source_root,work_root)
  s$action<-"freeze";s$freeze$config_id<-"A001";execute_submission(s,source_root,work_root)
  stopifnot(file.exists(file.path(work_root,"config","fixed_parameters.R")))
  s$action<-"run";s$modules<-c("B","C");runs<-execute_submission(s,source_root,work_root)
  s$action<-"report";execute_submission(s,source_root,work_root)
  result<-sim_collect(runs$C)
  stopifnot(nrow(result$effects)==6L,all(result$effects$available_effect_repetitions==2L),
    identical(sort(unique(result$effects$null_count)),c(0L,4L,8L)),sum(result$status$reused)==0L,
    all(result$status$point_ok),all(result$status$interval_ok),all(result$status$selection_ok))
  for(q in s$null_counts) {
    jobs<-result$plan$tasks[result$plan$tasks$null_count==q,,drop=FALSE]
    rec<-lapply(jobs$key,function(k)readRDS(file.path(runs$C,"fits",paste0(k,".rds"))))
    direct<-sim_effect_summary(rec,2L)$metrics
    for(target in c("joint","individual")) {
      got<-result$effects[result$effects$null_count==q&result$effects$target==target,]
      want<-direct[direct$target==target,]
      stopifnot(abs(got$RMSE-want$RMSE)<1e-14,abs(got$RMS_Bias-want$RMS_Bias)<1e-14)
    }
  }
  stopifnot(length(result$or_errors)==6L)
  table<-read.csv(file.path(runs$C,"analysis","tables","Table_C_effect_accuracy.csv"),stringsAsFactors=FALSE)
  stopifnot(all(c("null_count","p")%in%names(table)),nrow(table)==6L)
  figures<-list.files(file.path(runs$C,"analysis","figures"))
  stopifnot(sum(grepl("Figure_C_joint_page",figures))==2L,sum(grepl("Figure_C_individual_page",figures))==2L,
    all(paste0(c("Figure_C_effect_accuracy","Figure_C_uncertainty","Figure_C_OR_errors"),".png")%in%figures))
  # Different B/C overrides are resolved independently and recorded explicitly.
  alt<-s;alt$common_fixed<-list(N=500L);alt$C_overrides<-list(N=700L)
  cfg<-sim_submission_config(alt,source_root,work_root)
  stopifnot(sim_resolve_fixed(cfg,work_root,"B")$parameters$N==500L,
    sim_resolve_fixed(cfg,work_root,"C")$parameters$N==700L)
  fitpaths<-list.files(file.path(runs$C,"fits"),full.names=TRUE);hash<-tools::md5sum(fitpaths)
  s$action<-"run";s$modules<-"C";execute_submission(s,source_root,work_root)
  stopifnot(identical(hash,tools::md5sum(fitpaths)),identical(before,tools::md5sum(files)),
    identical(files,list.files(source_root,recursive=TRUE,full.names=TRUE,all.files=TRUE)[!file.info(list.files(source_root,recursive=TRUE,full.names=TRUE,all.files=TRUE))$isdir]))
  # Changing plot selection affects figures only; the tables retain all q values.
  s$action<-"report";s$report_overrides$C_null_counts<-4L
  execute_submission(s,source_root,work_root)
  table2<-read.csv(file.path(runs$C,"analysis","tables","all_effect_metrics.csv"))
  stopifnot(nrow(table2)==6L)
  s$report_overrides$C_null_counts<-NULL;execute_submission(s,source_root,work_root)
  kdr_log("Server workflow checks passed: external workspace, A/freeze/B/C, C effects, overrides, reuse, resume and unchanged source.")
  invisible(list(A=a,B=runs$B,C=runs$C))
}
