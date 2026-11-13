# Fixed-configuration assessment across overall exposure dependence structures.
sim_D_levels <- function(cfg) {
  d<-cfg$D$correlation_levels
  required<-c("correlation_id","correlation_label","target_mean","min_abs","mean_tolerance")
  if(!is.data.frame(d)||!all(required%in%names(d))||!nrow(d))
    stop("D_correlation_levels requires correlation_id, correlation_label, target_mean, min_abs and mean_tolerance.")
  d<-d[,required,drop=FALSE]
  d$correlation_id<-as.character(d$correlation_id);d$correlation_label<-as.character(d$correlation_label)
  if(anyNA(d)||anyDuplicated(d$correlation_id)||any(d$correlation_id=="D0")||anyDuplicated(d$target_mean)||
     any(!grepl("^[A-Za-z][A-Za-z0-9_]*$",d$correlation_id))||
     any(!is.finite(d$target_mean)|d$target_mean<=0|d$target_mean>=1)||
     any(!is.finite(d$min_abs)|d$min_abs<0|d$min_abs>d$target_mean)||
     any(!is.finite(d$mean_tolerance)|d$mean_tolerance<=0))
    stop("Use unique non-D0 IDs and valid target means, lower bounds and tolerances.")
  r<-cfg$D$repetitions
  if(length(r)!=1L||!is.finite(r)||r<1||r!=floor(r))stop("D repetitions must be a positive integer.")
  c<-cfg$D$calibration
  for(n in c("reference_groups","calibration_groups","validation_groups")) {
    v<-c[[n]]
    if(length(v)!=1L||!is.finite(v)||v<100L||v!=floor(v))
      stop("D calibration sample sizes must be integers of at least 100 matched sets.")
  }
  if(length(c$maximum_mix)!=1L||!is.finite(c$maximum_mix)||c$maximum_mix<=0||c$maximum_mix>=1||
     length(c$root_tolerance)!=1L||!is.finite(c$root_tolerance)||c$root_tolerance<=0)
    stop("D calibration requires 0 < maximum_mix < 1 and a positive root_tolerance.")
  d<-d[order(d$target_mean),,drop=FALSE];rownames(d)<-NULL
  d$correlation_mix<-if(is.null(cfg$D$copula))NA_real_ else
    cfg$D$copula$levels$correlation_mix[match(d$correlation_id,cfg$D$copula$levels$correlation_id)]
  d
}

sim_D_level_config <- function(cfg,level) {
  if(is.null(cfg$D$copula))stop("Prepare D calibration before generating data.")
  i<-match(level$correlation_id,cfg$D$copula$levels$correlation_id)
  if(is.na(i))stop("D correlation setting is not present in the calibration.")
  cfg$baseline$dgp$copula_eta<-cfg$D$copula$levels$correlation_mix[i]
  cfg$batch<-"correlation"
  cfg
}

sim_plan_D <- function(cfg,resolved) {
  levels<-sim_D_levels(cfg)
  cf<-cfg;cf$B$repetitions<-cfg$D$repetitions
  base<-sim_plan(cf,"B",resolved);tasks<-configs<-list()
  for(i in seq_len(nrow(levels))) {
    level<-levels[i,,drop=FALSE];d<-base$tasks
    for(n in names(level))d[[n]]<-level[[n]]
    d$config_id<-paste0(level$correlation_id,"_",d$config_id)
    d$key<-paste0(level$correlation_id,"_",d$key)
    d$kernel_id<-paste0(level$correlation_id,"_",d$kernel_id)
    d$kernel_key<-paste0(level$correlation_id,"_",d$kernel_key)
    d$aux_key<-paste0(level$correlation_id,"_",d$aux_key)
    d$module<-"D";d$batch<-"correlation"
    tasks[[i]]<-d
    a<-base$design$configurations;a$config_id<-paste0(level$correlation_id,"_",a$config_id)
    configs[[i]]<-cbind(a,level)
  }
  tasks<-kdr_bind(tasks);tasks<-tasks[order(tasks$replicate,tasks$target_mean,tasks$method,tasks$scenario_gf),]
  rownames(tasks)<-NULL;stopifnot(!anyDuplicated(tasks$key))
  list(tasks=tasks,design=list(configurations=kdr_bind(configs),membership=data.frame(),correlation_levels=levels))
}

sim_data_key <- function(job) paste0(if(is.null(job$correlation_id))"" else paste0(job$correlation_id,"_"),
  "data_r",sprintf("%03d",job$replicate))

sim_D_cache <- function(cache,id) {
  if(length(id)!=1L||!grepl("^[A-Za-z][A-Za-z0-9_]*$",id))stop("Invalid correlation cache ID")
  path<-file.path(cache,id);dir.create(path,recursive=TRUE,showWarnings=FALSE)
  writeLines("owned",file.path(path,".owned_simulation_cache"));path
}

sim_D_data_jobs <- function(plan) lapply(sim_groups(plan,c("correlation_id","replicate")),function(ix) {
  a<-plan[ix,,drop=FALSE]
  list(replicate=a$replicate[1],L=max(a$L),N=max(a$N),correlation_id=a$correlation_id[1],
    target_mean=a$target_mean[1],correlation_mix=a$correlation_mix[1],correlation_label=a$correlation_label[1])
})

sim_evaluation_path <- function(run,p,correlation_id=NULL) file.path(run,"evaluation",
  paste0(if(is.null(correlation_id))"" else paste0(correlation_id,"_"),"p",p,".rds"))

sim_prepare_D_evaluations <- function(cfg,run) {
  levels<-sim_D_levels(cfg);cfg$batch<-"correlation"
  reference<-readRDS(file.path(run,"baseline_B","evaluation","p8.rds"))
  kdr_save(reference,file.path(run,"evaluation","common_reference.rds"))
  for(i in seq_len(nrow(levels))) {
    level<-levels[i,,drop=FALSE];path<-sim_evaluation_path(run,8L,level$correlation_id)
    if(!is.null(sim_read(path)))next
    kdr_log("Preparing evaluation distribution:",level$correlation_id,level$correlation_label)
    cf<-sim_D_level_config(cfg,level)
    d<-kdr_generate(cf$evaluation$groups,8L,cf,kdr_seed(cf,"evaluation"),events=FALSE)
    rows<-match(seq_len(cf$evaluation$function_points),d$group)
    out<-reference;out$individual<-d$A[rows,,drop=FALSE]
    out$true_individual<-kdr_truth(out$individual,cf)-kdr_truth(matrix(out$reference,1,8),cf)
    out$overall_correlation<-cor(d$A)
    out$correlation_sample_rows<-nrow(d$A);out$correlation_sample_sets<-d$n
    out$correlation_setting<-level
    if(!sim_D_check_target(sim_D_cor_summary(d$A),level))
      stop("Independent evaluation sample does not meet D correlation targets: ",level$correlation_id)
    out$marginal_summary<-data.frame(variable=colnames(d$A),mean=colMeans(d$A),SD=apply(d$A,2,sd))
    out$interpretation<-sim_interpretation_context(out,cf)
    kdr_save(out,path)
    kdr_log("Overall mean absolute Pearson correlation:",
      format(mean(abs(out$overall_correlation[upper.tri(out$overall_correlation)])),digits=4))
  }
  invisible(NULL)
}
