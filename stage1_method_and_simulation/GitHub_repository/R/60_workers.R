sim_event <- function(run,kind,key,ok,seconds,job=NULL,error="",interval_ok=FALSE,selection_ok=FALSE) {
  kdr_save(list(kind=kind,key=key,ok=ok,elapsed_seconds=seconds,job=job,error=error,
    interval_ok=interval_ok,selection_ok=selection_ok,completed_at=Sys.time()),
    file.path(run,"events",paste0(kind,"_",key,".rds")))
}
sim_worker_context <- function(run,key,cfg,slot) {
  e<-new.env(parent=emptyenv());e$key<-key;e$slot<-slot;e$started<-Sys.time();e$updated<-e$started
  e$path<-file.path(run,"status",paste0(key,".rds"));e$refresh<-cfg$execution$worker_status_seconds
  e$stage<-"Starting";e$number<-0L;e$detail<-"";e$i<-0L;e$n<-1L
  options(kdr.worker=e);e
}
sim_prepare_worker <- function(job,cfg,run,cache,slot) {
  key<-sim_data_key(job)
  if(!is.null(job$correlation_id)) {
    cfg<-sim_D_level_config(cfg,job);cache<-sim_D_cache(cache,job$correlation_id)
  }
  sim_worker_context(run,key,cfg,slot)
  start<-proc.time()[3];sim_resource_context(run,"data",key,job,cache,"Data generation")
  ans<-kdr_attempt({
    total<-ceiling(job$L/cfg$execution$block_groups)+ceiling(job$N/cfg$execution$block_groups);i<-0L
    sim_major(1,"Preparing independent data blocks")
    for(role in c("learning","estimation")) {
      size<-if(role=="learning")job$L else job$N
      for(b in seq_len(ceiling(size/cfg$execution$block_groups))) {
        path<-sim_pool_file(cache,role,job$replicate,b)
        if(is.null(sim_read(path)))kdr_save(kdr_generate(cfg$execution$block_groups,8L,cfg,
          kdr_seed(cfg,"data_pool",role,job$replicate,b)),path)
        i<-i+1L;sim_status(paste("Generated",role,"block",b),i,total,TRUE)
      }
    }
  })
  sim_event(run,"data",key,ans$ok,proc.time()[3]-start,job,ans$error)
  sim_resource_context(run,"idle",paste0("idle_",Sys.getpid()),phase="Idle")
  list(key=key,ok=ans$ok,error=ans$error)
}
sim_setup <- function(cfg,run,job,cache) {
  model<-file.path(run,"kernels",paste0(job$kernel_key,".rds"))
  record_path<-file.path(run,"kernel_records",paste0(job$kernel_key,".rds"))
  old<-sim_read(record_path)
  if(!is.null(old)&&isTRUE(old$ok)&&file.exists(model))return(readRDS(model))
  if(!is.null(old)&&!isTRUE(old$ok)&&!cfg$execution$retry_failed)stop(old$error)
  start<-proc.time()[3];sim_resource_context(run,"kernel",job$kernel_key,job,phase="Kernel learning")
  ans<-kdr_attempt({
    sim_major(1,"Independent learning sample")
    d<-sim_learning_data(cfg,cache,job$L,job$p,job$replicate)
    sim_major(2,"Kernel geometry and distance learning")
    s<-kdr_geometry(d,cfg,kdr_seed(cfg,"centers",job$replicate))
    s$learning<-kdr_learning_fit(d,s,cfg);s$rho<-s$learning$rho
    s$B<-kdr_whiten(kdr_kernel_raw(d$A,s),d$group);s$dimension<-ncol(s$B)
    s$penalty<-kdr_kernel_penalty(d,s,cfg,kdr_seed(cfg,"penalty",job$replicate))
    if(isTRUE(cfg$nuisance$data_derived)) {
      s$nuisance_raw_spec<-sim_E_nuisance_spec(d,s,cfg)
      cfg$nuisance$raw_spec<-s$nuisance_raw_spec
    }
    s$nuisance<-kdr_nuisance_spec(d,cfg)
    d$X<-kdr_basis(d$A,s);d$pairs$Z<-d$X[d$pairs$plus,,drop=FALSE]-d$X[d$pairs$minus,,drop=FALSE]
    cf<-cfg;cf$nuisance$folds<-1L;d$fold[]<-0L;each<-function(FUN,progress=FALSE)FUN(d)
    f<-kdr_f_fit(each,cf,d$n);d$fold[]<-1L
    aux<-list(g=list(list(gamma=tail(s$learning$beta,2))),f=f)
    st<-kdr_dr_state(numeric(s$dimension),each,aux,rep(1,9),d$n)
    s$T<-kdr_inverse(st$J);s$learning_N<-job$L;s$kappa<-job$kappa
    s$learning_seconds<-unname(proc.time()[3]-start)
    kdr_save(s,model);s
  })
  seconds<-unname(proc.time()[3]-start)
  kdr_save(list(job=job,ok=ans$ok,error=ans$error,elapsed_seconds=seconds),record_path)
  sim_event(run,"kernel",job$kernel_key,ans$ok,seconds,job,ans$error)
  if(!ans$ok)stop(ans$error)
  ans$value
}
sim_complete <- function(record,retry) {
  if(is.null(record))return(FALSE)
  if(!retry)return(TRUE)
  detection<-record$interpretation$detection
  pairs<-detection$value$pairs
  interaction_ok<-isTRUE(detection$ok) && !is.null(pairs) && nrow(pairs)>0L &&
    all(!is.na(pairs$test_ok) & pairs$test_ok)
  isTRUE(record$ok) && isTRUE(record$interval_ok) && isTRUE(record$interpretation_ok) &&
    interaction_ok && (!isTRUE(record$selection_required)||isTRUE(record$selection$valid))
}
sim_save_fit <- function(run,record) {
  kdr_save(record,file.path(run,"fits",paste0(record$job$key,".rds")))
  sim_event(run,"fit",record$job$key,isTRUE(record$ok),record$elapsed_seconds,record$job,
    record$error,isTRUE(record$interval_ok),isTRUE(record$selection$valid))
}
sim_failed_fit <- function(run,job,error,stage,seconds=0) sim_save_fit(run,
  list(job=job,ok=FALSE,interval_ok=FALSE,error=error,failure_stage=stage,
       elapsed_seconds=seconds,completed_at=Sys.time(),warnings=character()))

sim_auxiliary <- function(raw,each,s,cfg,job,work,run,gm,fw) {
  key<-paste0(job$aux_key,"_",gm,"_f",as.integer(fw));path<-file.path(run,"auxiliaries",paste0(key,".rds"))
  old<-sim_read(path);if(!is.null(old)&&isTRUE(old$ok))return(old)
  sim_resource_context(run,"auxiliary",key,job,work,"Nuisance estimation")
  start<-proc.time()[3];sim_major(3,"Cross-fitted nuisance estimation")
  ans<-kdr_attempt({
    g<-kdr_g_fit(raw,s,cfg,gm,file.path(work,"g"))
    fp<-file.path(work,paste0("f_",fw,".rds"));f<-sim_read(fp)
    if(is.null(f)){f<-kdr_f_fit(each,cfg,job$N,wrong=fw);kdr_save(f,fp)}
    list(g=g,f=f)
  })
  result<-list(ok=ans$ok,error=ans$error,aux=ans$value,seconds=unname(proc.time()[3]-start),job=job)
  kdr_save(result,path);sim_resource_checkpoint(run,"auxiliary",key,work)
  sim_event(run,"auxiliary",key,ans$ok,result$seconds,job,ans$error)
  if(!ans$ok)stop(ans$error)
  result
}

sim_one_fit <- function(job,cfg,s,each,aux,ev,run,work,setup_seconds=0) {
  start<-proc.time()[3];warnings<-character();gc(FALSE)
  sim_resource_context(run,"fit",job$key,job,work,"Point estimation")
  ans<-withCallingHandlers(kdr_attempt({
    sim_major(4,paste("Estimating",sim_label(job$method),job$scenario_gf))
    cp<-file.path(run,"points",paste0(job$key,".rds"));point<-sim_read(cp)
    if(is.null(point)) {
      t<-proc.time()[3]
      fit<-if(job$method=="KDR")kdr_fit(each,aux$aux,s,cfg,job$N) else
        kdr_cl_model(each,s,cfg,wrong=job$scenario_gf=="wrong")
      point<-list(fit=fit,seconds=unname(proc.time()[3]-t));kdr_save(point,cp)
    }
    fit<-point$fit
    sim_major(5,"Covariance and effect evaluation");sim_resource_phase("Covariance")
    t<-proc.time()[3]
    infer<-if(job$method=="KDR")kdr_attempt(kdr_dr_covariance(fit,each,aux$aux,s,cfg,job$N)) else
      list(ok=TRUE,value=list(V=fit$V,scope="Matched-set conditional likelihood sandwich covariance"),error="")
    infer_seconds<-unname(proc.time()[3]-t);V<-if(infer$ok)infer$value$V else NULL
    t<-proc.time()[3];out<-kdr_outputs(fit,s,ev,cfg,V);effect_seconds<-unname(proc.time()[3]-t)
    sim_resource_phase("Contributions and interaction inference");t<-proc.time()[3]
    ext<-kdr_attempt(sim_extended_fit(fit,s,ev,cfg,V))
    interpretation_seconds<-unname(proc.time()[3]-t)
    required<-job$method=="KDR" && job$scenario_gf=="G1"
    sim_major(6,if(required)"Variable tests and checkpoint" else "Saving completed fit")
    sim_resource_phase("Variable testing");t<-proc.time()[3]
    sel<-if(required)kdr_attempt(kdr_selection(fit,s,ev,cfg,V,job)) else NULL
    selection_seconds<-unname(proc.time()[3]-t)
    selection<-if(!required)NULL else if(sel$ok)sel$value else list(table=data.frame(),valid=FALSE,error=sel$error)
    list(fit=fit,outputs=out,interpretation=ext$value,interpretation_ok=ext$ok,
      interpretation_error=ext$error,interpretation_seconds=interpretation_seconds,selection=selection,selection_required=required,
      interval_ok=infer$ok,interval_error=infer$error,inference=infer$value,
      dimension=s$dimension,kernel_key=if(job$method=="KDR")job$kernel_key else NA_character_,
      point_seconds=point$seconds,inference_seconds=infer_seconds,effects_seconds=effect_seconds,
      selection_seconds=selection_seconds,learning_seconds=if(job$method=="KDR")s$learning_seconds else setup_seconds,
      auxiliary_seconds=if(job$method=="KDR")aux$seconds else 0,
      computation_seconds=point$seconds+infer_seconds+effect_seconds+selection_seconds+interpretation_seconds)
  }),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  record<-c(list(job=job,ok=ans$ok,error=ans$error,elapsed_seconds=unname(proc.time()[3]-start),
    completed_at=Sys.time(),warnings=unique(warnings),pid=Sys.getpid(),reused=FALSE),ans$value)
  sim_resource_checkpoint(run,"fit",job$key,work);sim_save_fit(run,record)
  sim_status("Result saved",1,1,TRUE);invisible(record)
}

sim_worker <- function(jobs,cfg,run,cache,slot) {
  first<-as.list(jobs[1,,drop=FALSE]);key<-first$kernel_key
  logfile<-file(file.path(run,"logs",paste0(key,".txt")),"at",encoding="UTF-8")
  sink(logfile);on.exit({sink();close(logfile)},add=TRUE)
  sim_worker_context(run,key,cfg,slot)
  on.exit({options(kdr.worker=NULL);sim_resource_context(run,"idle",paste0("idle_",Sys.getpid()),phase="Idle")},add=TRUE)
  if(!is.null(first$correlation_id))cache<-sim_D_cache(cache,first$correlation_id)
  cc<-sim_job_config(cfg,first);work<-file.path(cache,"work",key)
  need<-vapply(jobs$key,function(k)!sim_complete(sim_read(file.path(run,"fits",paste0(k,".rds"))),cfg$execution$retry_failed),logical(1))
  if(!any(need))return(list(key=key,ok=TRUE,error=""))
  jobs<-jobs[need,,drop=FALSE]
  setup<-if(any(jobs$method=="KDR"))kdr_attempt(sim_setup(cc,run,first,cache)) else list(ok=FALSE,error="Not required")
  ev<-readRDS(sim_evaluation_path(run,first$p,first$correlation_id))
  comparators<-list()
  for(N in sort(unique(jobs$N))) {
    sub<-jobs[jobs$N==N,,drop=FALSE];sub<-sub[order(sub$method,sub$scenario_gf,sub$lambda),]
    job<-as.list(sub[1,,drop=FALSE]);cc<-sim_job_config(cfg,job)
    nw<-file.path(work,paste0("N",N));raw<-sim_raw(cc,cache,N,job$p,job$replicate)
    for(i in seq_len(nrow(sub))) {
      jj<-as.list(sub[i,,drop=FALSE]);cf<-sim_job_config(cfg,jj)
      z<-kdr_attempt({
        if(jj$method=="KDR") {
          if(!setup$ok)stop(setup$error)
          s<-setup$value;each<-kdr_feature_iterator(raw,s,file.path(nw,"kernel_features"),cf)
          gm<-if(jj$scenario_gf%in%c("G1","G2"))"correct" else "wrong"
          fw<-jj$scenario_gf%in%c("G2","G4")
          aux<-sim_auxiliary(raw,each,s,cf,jj,nw,run,gm,fw);setup_time<-0
        } else {
          if(is.null(comparators[[jj$method]])) {
            t<-proc.time()[3];d<-sim_learning_data(cf,cache,jj$L,jj$p,jj$replicate)
            spec<-kdr_comparator_setup(d,cf,if(jj$method=="SplineCLR")"spline" else "linear")
            comparators[[jj$method]]<-list(spec=spec,seconds=unname(proc.time()[3]-t))
            kdr_save(spec,file.path(run,"comparators",paste0(jj$kernel_key,"_",jj$method,".rds")))
          }
          s<-comparators[[jj$method]]$spec;setup_time<-comparators[[jj$method]]$seconds
          each<-kdr_feature_iterator(raw,s,file.path(nw,paste0(jj$method,"_features")),cf);aux<-NULL
        }
        sim_one_fit(jj,cf,s,each,aux,ev,run,nw,setup_time)
      })
      if(!z$ok)sim_failed_fit(run,jj,z$error,"setup_or_nuisance")
    }
    if(!cfg$execution$keep_work_cache)sim_safe_remove(nw,cache)
  }
  gc(FALSE);list(key=key,ok=TRUE,error="")
}
