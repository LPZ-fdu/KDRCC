prepare_results <- function(run,report_overrides=list()) {
  run<-normalizePath(run,winslash="/",mustWork=TRUE)
  lock<-file.path(run,"report_lock")
  if(!dir.create(lock,showWarnings=FALSE))stop("Report lock exists; another report may be running.")
  on.exit(sim_unlock(lock,run),add=TRUE)
  analysis<-file.path(run,"analysis");tables<-file.path(analysis,"tables");figures<-file.path(analysis,"figures")
  dir.create(tables,recursive=TRUE,showWarnings=FALSE);dir.create(figures,recursive=TRUE,showWarnings=FALSE)
  old_captions<-sim_read_resource_csv(file.path(analysis,"Figure_and_table_captions.csv"))
  kdr_log("Reading saved results:",basename(run))
  result<-sim_collect(run);result$config$report<-sim_merge(result$config$report,report_overrides)
  if(!isTRUE(capabilities("cairo")))stop("This R build lacks Cairo PNG support. Run reporting on a Cairo-enabled computer; fits need not be rerun.")
  result$config$report<-sim_report_font(result$config$report)
  module<-readRDS(file.path(run,"manifest.rds"))$module
  captions<-new.env(parent=emptyenv())
  emit<-function(d,id,title,note)sim_export_table(d,id,title,note,tables,captions)
  availability<-aggregate(result$status[c("completed","point_ok","interval_ok","selection_ok","interpretation_ok","interaction_ok")],
    by=result$status[c("config_id","method","scenario_gf","null_count","p")],FUN=sum)
  planned<-aggregate(rep(1L,nrow(result$status)),by=result$status[c("config_id","method","scenario_gf","null_count","p")],FUN=sum)
  names(planned)[ncol(planned)]<-"planned_repetitions"
  availability<-merge(planned,availability,by=c("config_id","method","scenario_gf","null_count","p"),sort=FALSE)
  emit(availability,paste0("Table_",module,"_availability"),"Availability of estimates and tests",
    "Counts refer to complete independent repetitions. Completed records include failures. Point, covariance and variable-test availability are reported separately. Selection is requested for KDR-CC G1 only. Failed and unavailable tests retain explicit status and availability counts.")
  kdr_csv(result$status,file.path(tables,"task_status.csv"))
  emit(result$plan$design$configurations,paste0("Table_",module,"_configurations"),"Model configurations",
    "N and L count complete matched sets. M is the number of Gaussian-kernel centers. Kappa is the coefficient penalty during outcome-supervised kernel learning. Lambda multiplies the normalized mixed-derivative penalty in final estimation.")
  d<-result$effects
  if(nrow(d)) {
    idcols<-c("config_id","method","scenario_gf","null_count","p","target","available_effect_repetitions")
    emit(d[c(idcols,"RMS_Bias","MCSE_RMS_Bias","RMSE","MCSE_RMSE","correct_significant_percent","MCSE_correct_significant_percent","wrong_significant_percent")],
      paste0("Table_",module,"_effect_accuracy"),"Accuracy of estimated exposure-effect contrasts",
      "RMS bias is the square root of the mean squared pointwise bias across common evaluation points. RMSE averages squared error over repetitions and points before taking the square root. Both are on the log OR scale. MCSE for these aggregate measures uses deletion of one entire repetition. Correct-direction detection requires the nominal 95% interval to exclude zero on the true-effect side; exact zero-reference contrasts are excluded. Wrong-direction significant detection is also reported. Independent repetitions determine Monte Carlo uncertainty.")
    emit(d[c(idcols,"mean_OR_absolute_relative_error","median_OR_absolute_relative_error","OR_error_Q25","OR_error_Q75","OR_error_P95","MCSE_mean_OR_error","nonfinite_OR_errors")],
      paste0("Table_",module,"_OR_errors"),"Distribution of absolute relative odds-ratio errors",
      "Error (%) = 100 * abs(exp(estimated log OR - true log OR) - 1). Distributional summaries pool evaluation points and repetitions. MCSE of the mean is computed across repetition-level means. No acceptance threshold or clipping of large errors is used.")
    emit(d[c(idcols,"available_interval_repetitions","RMS_empirical_SD","RMS_model_SE","SE_ratio","coverage_percent","MCSE_coverage_percent","mean_interval_width")],
      paste0("Table_",module,"_uncertainty"),"Uncertainty estimation and interval calibration",
      "Empirical SD is computed across independent repetitions at each fixed point and summarized by its root mean square. Model SE is summarized on the same squared scale. SE ratio is model SE divided by empirical SD; its reference value is one. Coverage is the average pointwise coverage of nominal 95% intervals and is a secondary assessment. Intervals condition on the independently learned kernel; approximation and regularization bias are assessed separately through their empirical errors and coverage. Empirical variability across repetitions can also include kernel learning. Width and coverage must be interpreted together with bias and estimation error.")
    kdr_csv(d,file.path(tables,"all_effect_metrics.csv"))
    kdr_csv(result$repetitions,file.path(tables,"effect_metrics_by_repetition.csv"))
    if(result$config$report$save_pointwise_tables)kdr_csv(result$points,file.path(tables,"pointwise_effect_summary.csv"))
  }
  if(nrow(result$cost)) {
    cost<-result$cost;groupcols<-c("config_id","method","scenario_gf",sim_parameter_names,"p")
    summary<-kdr_bind(lapply(sim_groups(cost,groupcols),function(ix){
      z<-cost[ix,];cbind(z[1,groupcols,drop=FALSE],data.frame(n=nrow(z),
        mean_learning_seconds=sim_mean(z$learning_seconds),mean_auxiliary_seconds=sim_mean(z$auxiliary_seconds),
        mean_estimation_seconds=sim_mean(z$estimation_seconds),mean_component_total_seconds=sim_mean(z$component_total_seconds),
        median_component_total_seconds=median(z$component_total_seconds,na.rm=TRUE),
        mean_sampled_peak_memory_GB=sim_mean(z$sampled_peak_memory_GB),
        memory_available_repetitions=sum(is.finite(z$sampled_peak_memory_GB)),mean_cache_GB=sim_mean(z$cache_GB)))
    }))
    emit(summary,paste0("Table_",module,"_resources"),"Computational time and storage requirements",
      "Times are recorded elapsed seconds for individual components on the running hardware. Shared learning and nuisance components are reused. Row-level component sums describe a fit in isolation; actual total elapsed time is recorded in the run log. G-family reuse can reduce observed nuisance timings. Memory is the sampled process working-set maximum across the relevant stages. R may retain allocator memory between tasks. Cache is the measured working-cache size before removal. Missing measurements are reported as NA, never zero.")
    kdr_csv(cost,file.path(tables,"resources_by_repetition.csv"))
  }
  selection<-sim_selection_summary(result$selection,result$status)
  emit(selection$overall,paste0("Table_",module,"_selection_rates"),"Mean identification and false-positive rates",
    "Rates are computed within each usable repetition and then averaged. TPR uses eight active exposures; FPR uses the candidate null exposures and is undefined when none are present. Monte Carlo intervals use repetition-level rates. Usable and planned repetitions are shown. These are empirical operating characteristics of the regularized test conditional on the learned kernel.")
  emit(selection$variables,paste0("Table_",module,"_active_exposure_rates"),"Identification rates for each active exposure",
    "Identification requires a Holm-adjusted numerical quadratic-form tail probability below 0.05. Intervals are 95% Wilson Monte Carlo intervals across usable repetitions. Humidity acts through a temperature interaction, without a main effect.")
  if(nrow(result$selection))kdr_csv(result$selection,file.path(tables,"variable_tests_by_repetition.csv"))
  sim_report_interpretation(run,result,tables,figures,captions,module)
  sim_report_model_diagnostics(run,result,captions,module)
  if(module=="A") {
    recommendation<-sim_recommend(result,run)
    if(!is.null(recommendation))emit(recommendation$ranking,"Table_A_configuration_ranking","Configuration recommendation scores",
      "Eligibility and weights were specified in 00_parameters.R. The score combines joint and observed-combination normalized RMSE, normalized RMS bias, mean OR error, contribution-share error, overall interaction miss rate, correct-direction miss rate and wrong-direction rate. Weights are recorded in the configuration. Coverage is a secondary measure. B and C read the explicitly frozen configuration and record any manual overrides.")
    sim_plot_A(result,figures,captions)
  } else if(module=="B") {
    sim_plot_B(result,figures,captions)
    sim_report_B_responses(run,result,tables,figures,captions)
  } else {
    sim_plot_C(result,selection,figures,captions)
    sim_plot_C_effects(result,figures,captions)
  }
  cap<-kdr_bind(lapply(sort(ls(captions)),function(id){a<-captions[[id]];data.frame(id=id,type=a$type,file=a$file,title=a$title,note=a$note)}))
  # Remove obsolete generated pages after a successful render, never raw results.
  if(nrow(old_captions)&&"file"%in%names(old_captions)) {
    stale<-setdiff(old_captions$file,cap$file)
    stale<-stale[grepl("^figures/Figure_[A-Za-z0-9_]+[.]png$",stale)]
    for(f in stale) {
      path<-file.path(analysis,f)
      if(file.exists(path)&&identical(dirname(normalizePath(path,winslash="/",mustWork=TRUE)),
        normalizePath(figures,winslash="/",mustWork=TRUE)))unlink(path)
    }
  }
  kdr_csv(cap,file.path(analysis,"Figure_and_table_captions.csv"))
  lines<-unlist(lapply(seq_len(nrow(cap)),function(i)c(paste(cap$id[i],cap$title[i]),cap$note[i],"")))
  writeLines(lines,file.path(analysis,"Figure_and_table_captions.txt"),useBytes=TRUE)
  html<-c("<!doctype html><html lang='en'><head><meta charset='utf-8'><title>KDR-CC simulation results</title><style>body{font-family:'Times New Roman',serif;margin:32px;max-width:1200px}img{max-width:100%}p{line-height:1.5}</style></head><body>",
    paste0("<h1>Module ",module," results</h1>"))
  for(i in seq_len(nrow(cap))) {
    html<-c(html,paste0("<h2>",sim_escape(cap$title[i]),"</h2>"),
      if(cap$type[i]=="Figure")paste0("<a href='",cap$file[i],"'><img src='",cap$file[i],"' alt='",sim_escape(cap$title[i]),"'></a>") else
        paste0("<p><a href='",cap$file[i],"'>Open formatted table</a> | <a href='",sub("html$","csv",cap$file[i]),"'>CSV</a></p>"),
      paste0("<p>",sim_escape(cap$note[i]),"</p>"))
  }
  writeLines(c(html,"</body></html>"),file.path(analysis,"index.html"),useBytes=TRUE)
  # Preserve the exact reporting code and options independently of fitting code.
  src<-getOption("kdr.root");if(!is.null(src)) {
    fs<-list.files(file.path(src,"R"),pattern="^(80|81|85|86|87|88|89|90).*R$",full.names=TRUE)
    dir.create(file.path(analysis,"report_source"),showWarnings=FALSE);file.copy(fs,file.path(analysis,"report_source"),overwrite=TRUE)
  }
  kdr_save(list(report_settings=result$config$report,generated_at=Sys.time()),file.path(analysis,"report_metadata.rds"))
  kdr_log("Report saved:",file.path(analysis,"index.html"));invisible(analysis)
}

freeze_configuration <- function(root,run=sim_pointer(root,"A"),config_id=NULL,manual_overrides=list()) {
  if(is.null(run))stop("No module A run found")
  pp<-readRDS(file.path(run,"plan.rds"));d<-pp$design$configurations
  if(is.null(config_id)) {
    rec<-sim_read(file.path(run,"analysis","recommended_configuration.rds"))
    if(is.null(rec))stop("No eligible recommendation. Review module A and specify a configuration ID explicitly.")
    config_id<-rec$config_id
  }
  if(!config_id%in%d$config_id)stop("Unknown configuration ID")
  pars<-as.list(d[d$config_id==config_id,sim_parameter_names,drop=FALSE]);pars<-sim_merge(pars,manual_overrides)
  if(length(setdiff(names(pars),sim_parameter_names)))stop("Unknown fixed parameter")
  value<-list(parameters=pars,source_run=basename(run),source_config_id=config_id,
    manual_overrides=manual_overrides,mechanism_id=object_id(readRDS(file.path(run,"config.rds"))$baseline),frozen_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"))
  target<-file.path(root,"config","fixed_parameters.R")
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  if(file.exists(target)) {
    archive<-file.path(root,"config","history");dir.create(archive,showWarnings=FALSE)
    file.copy(target,file.path(archive,paste0("fixed_",format(Sys.time(),"%Y%m%d_%H%M%S"),".R")))
  }
  writeLines(c("# Frozen configuration. Module B/C overrides are recorded separately.","fixed_parameters <- function() {",
    capture.output(dput(value,control="all")),"}"),target,useBytes=TRUE)
  kdr_log("Fixed configuration saved:",target);print(pars);invisible(value)
}
