sim_parameter_names <- c("N","L","M","kappa","lambda")
sim_parameter_key <- function(d) do.call(paste,c(d[sim_parameter_names],sep="|"))
sim_design_A <- function(cfg) {
  a<-cfg$A$reference;g<-cfg$A$grid;views<-list()
  make<-function(N=a$N,L=a$L,M=a$M,kappa=a$kappa,lambda=a$lambda)
    expand.grid(N=N,L=L,M=M,kappa=kappa,lambda=lambda,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  if(is.null(cfg$A$manual_design)) {
    for(v in intersect(sim_parameter_names,cfg$A$views))views[[v]]<-do.call(make,setNames(list(g[[v]]),v))
    if("kappa_lambda"%in%cfg$A$views)views$kappa_lambda<-make(kappa=g$crossed_kappa,lambda=g$crossed_lambda)
    if("L_M"%in%cfg$A$views)views$L_M<-make(L=g$crossed_L,M=g$crossed_M)
    views$reference<-make()
  } else views$custom<-cfg$A$manual_design[,sim_parameter_names,drop=FALSE]
  d<-kdr_bind(views);d<-d[!duplicated(sim_parameter_key(d)),,drop=FALSE];rownames(d)<-NULL
  d$config_id<-sprintf("A%03d",seq_len(nrow(d)))
  membership<-kdr_bind(lapply(names(views),function(v)data.frame(
    view=v,config_id=d$config_id[match(sim_parameter_key(views[[v]]),sim_parameter_key(d))])))
  list(configurations=d,membership=unique(membership))
}
sim_resolve_fixed <- function(cfg,root,module) {
  path<-file.path(root,cfg$configuration$fixed_file)
  e<-new.env();if(file.exists(path))sys.source(path,e,keep.source=FALSE)
  base<-if(exists("fixed_parameters",e,inherits=FALSE))e$fixed_parameters() else NULL
  if(!is.null(base)&&!identical(base$mechanism_id,object_id(cfg$baseline)))
    stop("Frozen configuration belongs to a different generating mechanism. Freeze the current module A configuration.")
  if(!is.null(base)&&!is.null(base$parameters))base<-base$parameters
  override<-cfg$configuration$overrides[[module]]
  if(is.null(override))override<-list()
  if(length(setdiff(names(override),sim_parameter_names)))stop("Unknown manual override parameter")
  effective<-sim_merge(if(is.null(base))list() else base,override)
  if(!all(sim_parameter_names%in%names(effective)))stop(
    "Module ",module,": freeze a configuration with 03_freeze_configuration.R, or supply all five manual overrides.")
  effective<-effective[sim_parameter_names]
  list(parameters=effective,base=base,manual_overrides=override,
    source_file=cfg$configuration$fixed_file,configuration_id=substr(object_id(effective),1,12))
}
sim_plan <- function(cfg,module,resolved=NULL) {
  if(module=="E")return(sim_plan_E(cfg,resolved))
  if(module=="D")return(sim_plan_D(cfg,resolved))
  if(module=="A") {
    design<-sim_design_A(cfg);d<-design$configurations
  } else {
    d<-as.data.frame(resolved$parameters,stringsAsFactors=FALSE)
    d$config_id<-paste0("F_",resolved$configuration_id);design<-list(configurations=d,membership=data.frame())
  }
  rows<-list()
  for(i in seq_len(nrow(d)))for(q in if(module=="C")cfg$C$null_counts else 0L)
    for(r in seq_len(cfg[[module]]$repetitions)) {
      a<-d[i,,drop=FALSE];a$null_count<-q;a$p<-8L+q;a$replicate<-r
      settings<-if(module=="B")data.frame(method=c(rep("KDR",4),rep("SplineCLR",2),rep("LinearCLR",2)),
        scenario_gf=c("G1","G2","G3","G4","correct","wrong","correct","wrong")) else
        data.frame(method="KDR",scenario_gf="G1")
      for(j in seq_len(nrow(settings))) {
        b<-a;b$method<-settings$method[j];b$scenario_gf<-settings$scenario_gf[j]
        b$module<-module;b$batch<-if(module=="A")"configuration" else "validation"
        b$key<-paste0(b$config_id,"_q",q,"_r",sprintf("%03d",r),"_",b$method,"_",b$scenario_gf)
        b$kernel_id<-substr(object_id(as.list(b[1,c("L","M","kappa","p")])),1,12)
        b$kernel_key<-paste0("K",b$kernel_id,"_r",sprintf("%03d",r))
        b$aux_key<-paste0(b$kernel_key,"_N",b$N)
        rows[[length(rows)+1L]]<-b
      }
    }
  plan<-kdr_bind(rows);plan<-plan[order(plan$replicate,plan$null_count,plan$config_id,plan$method,plan$scenario_gf),]
  rownames(plan)<-NULL;stopifnot(!anyDuplicated(plan$key))
  list(tasks=plan,design=design)
}
sim_job_config <- function(cfg,job) {
  cfg$N<-as.integer(job$N);cfg$L<-as.integer(job$L);cfg$M<-as.integer(job$M)
  cfg$kernel$kappa<-job$kappa;cfg$penalty<-list(a=job$lambda,reference_N=job$N)
  cfg$batch<-job$batch
  if(!is.null(job$correlation_id))cfg<-sim_D_level_config(cfg,job)
  cfg
}
sim_validate <- function(cfg,plan) {
  stopifnot(cfg$interpretation$contribution_points>=1L,
    cfg$interpretation$contribution_points<=cfg$evaluation$function_points,
    cfg$interpretation$single_points>=3L,cfg$interpretation$interaction_points>=3L,
    cfg$interpretation$alpha>0,cfg$interpretation$alpha<1)
  integer_positive<-function(x)length(x)>0 && all(is.finite(x)&x>=1&x==floor(x))
  stopifnot(nrow(plan)>0,integer_positive(plan$N),integer_positive(plan$L),integer_positive(plan$M),
    all(plan$N>=100),all(plan$L>=100),all(plan$M>=4 & plan$M<plan$L),
    all(is.finite(plan$kappa)&plan$kappa>0),all(is.finite(plan$lambda)&plan$lambda>=0),
    all(is.finite(plan$null_count)&plan$null_count>=0&plan$null_count==floor(plan$null_count)),
    integer_positive(cfg$execution$workers),integer_positive(cfg$execution$block_groups),
    integer_positive(cfg$nuisance$folds),integer_positive(cfg$evaluation$groups),
    integer_positive(cfg$evaluation$function_points),integer_positive(cfg$evaluation$coverage_points),
    integer_positive(cfg$selection$evaluation_points),cfg$evaluation$confidence>0,cfg$evaluation$confidence<1,
    is.finite(cfg$execution$raw_memory_MB),cfg$execution$raw_memory_MB>=0,
    is.finite(cfg$execution$feature_memory_MB),cfg$execution$feature_memory_MB>=0,
    cfg$kernel$rho_bounds[1]>0,cfg$kernel$rho_bounds[2]>cfg$kernel$rho_bounds[1],
    cfg$kernel$coverage_fraction>=0,cfg$kernel$coverage_fraction<1,
    cfg$nuisance$folds>=2,cfg$evaluation$groups>=cfg$evaluation$function_points,
    cfg$selection$evaluation_points<=cfg$evaluation$function_points,
    cfg$evaluation$coverage_points<=cfg$evaluation$function_points,
    cfg$selection$alpha>0,cfg$selection$alpha<1,identical(cfg$selection$adjust,"holm"),
    cfg$selection$accuracy>0,cfg$selection$accuracy<cfg$selection$alpha/max(plan$p)/100,
    cfg$noise$group_fraction>=0,cfg$noise$group_fraction<1)
  if(!requireNamespace("CompQuadForm",quietly=TRUE))stop("Install CompQuadForm with 00_install_dependencies.R first.")
  invisible(TRUE)
}
