# Each parallel bundle is an independent scenario/repetition with paired final fits.
sim_E_cv <- function(jobs,cfg,run,cache,candidates) {
  job<-as.list(jobs[1,,drop=FALSE]);id<-job$kernel_key
  directory<-file.path(run,"cross_validation",id);dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  K<-cfg$E$validation_folds;T<-cfg$E$tuning_sets
  outer<-with_seed(kdr_seed(cfg,"validation_partition",job$replicate),sample(rep(seq_len(K),length.out=T)))
  kdr_save(outer,file.path(directory,"fold_assignment.rds"))
  scores<-list();total<-nrow(candidates)*K;done<-0L
  context<-getOption("kdr.worker")
  if(!is.null(context)){context$validation_done<-0L;context$validation_total<-total}
  for(kid in unique(candidates$kernel_id)) {
    cs<-candidates[candidates$kernel_id==kid,,drop=FALSE];kj<-job
    for(n in sim_parameter_names)kj[[n]]<-cs[[n]][1]
    kj$kernel_key<-paste0(id,"_K",kid);kj$aux_key<-kj$kernel_key
    cf<-sim_E_config(cfg,job,kj)
    paths<-unlist(lapply(seq_len(K),function(k)file.path(directory,paste0(cs$candidate_id,"_fold",k,".rds"))))
    cached<-lapply(paths,sim_read)
    if(all(vapply(cached,function(a)!is.null(a)&&(!cfg$execution$retry_failed||isTRUE(a$ok)),logical(1)))) {
      scores<-c(scores,cached);done<-done+length(cached)
      if(!is.null(context))context$validation_done<-done
      sim_status("Validation scores restored",done,total,TRUE);next
    }
    setup<-kdr_attempt(sim_setup(cf,run,kj,cache))
    for(fold in seq_len(K)) {
      trainN<-sum(outer!=fold);cw<-file.path(cache,"work",id,paste0("K",kid,"_fold",fold))
      ir<-file.path(directory,paste0("K",kid,"_fold",fold));dir.create(ir,recursive=TRUE,showWarnings=FALSE)
      subpaths<-file.path(directory,paste0(cs$candidate_id,"_fold",fold,".rds"))
      old<-lapply(subpaths,sim_read)
      need<-vapply(old,function(a)is.null(a)||(cfg$execution$retry_failed&&!isTRUE(a$ok)),logical(1))
      if(any(need)) {
        train<-sim_E_raw(cf,cache,T,job$p,job$replicate,fold=fold,outer=outer)
        valid<-sim_E_raw(cf,cache,T,job$p,job$replicate,fold=fold,outer=outer,validation=TRUE)
        prep<-kdr_attempt({
          if(!setup$ok)stop(setup$error)
          s<-setup$value;cf$nuisance$raw_spec<-s$nuisance_raw_spec
          each<-kdr_feature_iterator(train,s,file.path(cw,"features"),cf)
          aj<-kj;aj$N<-as.integer(trainN);aj$aux_key<-paste0(kj$kernel_key,"_fold",fold)
          aux<-sim_auxiliary(train,each,s,cf,aj,cw,ir,"correct",FALSE)
          list(s=s,each=each,aux=aux)
        })
      }
      for(i in seq_len(nrow(cs))) {
        if(need[i]) {
          start<-proc.time()[3];warnings<-character()
          z<-withCallingHandlers(kdr_attempt({
            if(!prep$ok)stop(prep$error)
            cc<-cf;cc$N<-as.integer(trainN);cc$penalty<-list(a=cs$lambda[i],reference_N=trainN)
            sim_major(4,paste("Validation",cs$candidate_id[i],"fold",fold,"of",K))
            fit<-kdr_fit(prep$value$each,prep$value$aux$aux,prep$value$s,cc,trainN)
            sc<-sim_E_score(fit,prep$value$s,valid)
            c(sc,list(iterations=fit$iterations,gradient=fit$gradient))
          }),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
          a<-data.frame(candidate_id=cs$candidate_id[i],fold=fold,ok=z$ok,
            loss_sum=if(z$ok)z$value$loss_sum else NA_real_,n=sum(outer==fold),training_sets=trainN,
            elapsed_seconds=unname(proc.time()[3]-start),error=z$error,
            warnings=paste(unique(warnings),collapse=" | "),stringsAsFactors=FALSE)
          kdr_save(a,subpaths[i])
        } else a<-old[[i]]
        scores[[length(scores)+1L]]<-a;done<-done+1L
        if(!is.null(context))context$validation_done<-done
        sim_status(paste("Validation fits completed",done,"of",total),done,total,TRUE)
      }
      if(!cfg$execution$keep_work_cache)sim_safe_remove(cw,cache)
      gc(FALSE)
    }
  }
  scores<-kdr_bind(scores);kdr_csv(scores,file.path(directory,"fold_scores.csv"))
  selected<-sim_E_choose(scores,candidates,K,cfg$E$fallback_to_fixed)
  selected$scenario_id<-job$scenario_id;selected$replicate<-job$replicate
  selected$rule<-"Minimum matched-set-weighted held-out conditional negative log likelihood; all folds required"
  selected$final_data_role<-"estimation";selected$tuning_sets<-T;selected$kernel_learning_sets<-job$L
  kdr_save(selected,file.path(run,"choices",paste0(id,".rds")))
  kdr_csv(selected$ranking,file.path(directory,"candidate_ranking.csv"));selected
}

sim_E_worker <- function(jobs,cfg,run,cache,slot) {
  first<-as.list(jobs[1,,drop=FALSE]);id<-first$kernel_key;start<-proc.time()[3]
  logfile<-file(file.path(run,"logs",paste0(id,".txt")),"at",encoding="UTF-8")
  sink(logfile);on.exit({sink();close(logfile)},add=TRUE)
  sim_worker_context(run,id,cfg,slot)
  on.exit({options(kdr.worker=NULL);sim_resource_context(run,"idle",paste0("idle_",Sys.getpid()),phase="Idle")},add=TRUE)
  cache<-sim_D_cache(cache,id);cf<-sim_E_config(cfg,first,first)
  sim_resource_context(run,"selection",id,first,cache,"Configuration selection")
  plan<-readRDS(file.path(run,"plan.rds"));candidates<-plan$design$configurations
  need<-vapply(jobs$key,function(k)!sim_complete(sim_read(file.path(run,"fits",paste0(k,".rds"))),cfg$execution$retry_failed),logical(1))
  if(!any(need))return(list(key=id,ok=TRUE,error=""))
  sim_E_prepare_data(cf,cache,first)
  choicepath<-file.path(run,"choices",paste0(id,".rds"));choice<-sim_read(choicepath)
  if(is.null(choice)||cfg$execution$retry_failed)choice<-sim_E_cv(jobs,cf,run,cache,candidates)
  selection_seconds<-unname(proc.time()[3]-start)
  sim_resource_checkpoint(run,"selection",id,cache)
  ev<-readRDS(file.path(run,"evaluation",paste0(first$scenario_id,".rds")))
  raw<-sim_E_raw(cf,cache,first$N,first$p,first$replicate,role="estimation")
  actual<-list();fit_start<-proc.time()[3]
  for(i in which(need)) {
    j<-as.list(jobs[i,,drop=FALSE]);c<-if(j$strategy=="Fixed")candidates[candidates$is_fixed,,drop=FALSE] else choice$candidate
    for(n in sim_parameter_names)j[[n]]<-c[[n]][1]
    j$kernel_key<-paste0(id,"_K",c$kernel_id[1]);j$aux_key<-paste0(j$kernel_key,"_final")
    cc<-sim_E_config(cfg,j,j);work<-file.path(cache,"work",id,paste0("final_",c$kernel_id[1]))
    ans<-kdr_attempt({
      s<-sim_setup(cc,run,j,cache);cc$nuisance$raw_spec<-s$nuisance_raw_spec
      each<-kdr_feature_iterator(raw,s,file.path(work,"features"),cc)
      aux<-sim_auxiliary(raw,each,s,cc,j,work,run,"correct",FALSE)
      record<-sim_one_fit(j,cc,s,each,aux,ev,run,work)
      record$configuration_selection<-list(candidate_id=c$candidate_id[1],strategy=j$strategy,
        fallback=if(j$strategy=="Adaptive")choice$fallback else FALSE,
        selection_seconds=if(j$strategy=="Adaptive")selection_seconds else 0)
      sim_save_fit(run,record);record
    })
    if(!ans$ok)sim_failed_fit(run,j,ans$error,"E_final_estimation")
    actual[[length(actual)+1L]]<-as.data.frame(j,stringsAsFactors=FALSE)
  }
  choice$total_elapsed_seconds<-unname(proc.time()[3]-start)
  choice$selection_elapsed_seconds<-selection_seconds;choice$final_fits_elapsed_seconds<-unname(proc.time()[3]-fit_start)
  kdr_save(choice,choicepath)
  kdr_csv(kdr_bind(actual),file.path(run,"choices",paste0(id,"_final_parameters.csv")))
  if(!cfg$execution$keep_work_cache)sim_safe_remove(file.path(cache,"work",id),cache)
  gc(FALSE);list(key=id,ok=TRUE,error="")
}

run_module_E <- function(cfg,root) {
  started<-Sys.time();source_root<-getOption("kdr.root",root)
  resolved<-sim_resolve_fixed(cfg,root,"E");pp<-sim_plan_E(cfg,resolved);plan<-pp$tasks
  code<-list.files(file.path(source_root,"R"),pattern="[.]R$",full.names=TRUE)
  code<-code[as.integer(substr(basename(code),1,2))<80];hashes<-tools::md5sum(code);names(hashes)<-basename(code)
  id<-substr(object_id(list(science=sim_science(cfg),plan=pp,code=hashes)),1,12)
  folder<-paste0("module_E_",id);if(!cfg$execution$resume)folder<-paste0(folder,"_",format(Sys.time(),"%Y%m%d_%H%M%S"))
  run<-file.path(root,"results",folder)
  for(d in c("fits","kernels","kernel_records","auxiliaries","points","events","status","logs","evaluation","choices","cross_validation"))
    dir.create(file.path(run,d),recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(run,"run_lock")
  if(dir.exists(lock)&&cfg$execution$clear_stale_lock)sim_unlock(lock,run)
  if(!dir.create(lock,showWarnings=FALSE))stop("Run lock exists. Confirm no process is active before enabling clear_stale_lock.")
  on.exit(sim_unlock(lock,run),add=TRUE);writeLines(as.character(Sys.getpid()),file.path(lock,"owner_pid.txt"))
  logfile<-file(file.path(run,"logs",paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),".txt")),"at",encoding="UTF-8")
  sink(logfile,split=TRUE);on.exit({sink();close(logfile)},add=TRUE)
  if(any(pp$design$scenarios$dependence!="Original"))cfg<-sim_D_prepare_generator(cfg,run)
  cfg$batch<-"application_selection"
  manifest<-list(module="E",run_id=id,created_at=Sys.time(),effective_configuration=resolved,
    code_hash=hashes,planned_fits=nrow(plan),planned_validation_fits=nrow(plan)/2*pp$cv_fits_per_repetition)
  kdr_save(cfg,file.path(run,"config.rds"));kdr_save(pp,file.path(run,"plan.rds"));kdr_save(manifest,file.path(run,"manifest.rds"))
  kdr_csv(plan,file.path(run,"task_manifest.csv"));kdr_csv(pp$design$configurations,file.path(run,"parameter_configurations.csv"))
  sim_snapshot(source_root,run);capture.output(sessionInfo(),file=file.path(run,"sessionInfo.txt"))
  writeLines(folder,file.path(root,"results","LATEST_E.txt"))
  kdr_log("Module E | independent repetitions",nrow(plan)/2,"| final fits",nrow(plan),
    "| validation fits",manifest$planned_validation_fits,"| candidate configurations",nrow(pp$design$configurations))
  sim_E_prepare_evaluations(cfg,run,pp$design$scenarios)
  cache<-if(is.null(cfg$execution$cache_parent))file.path(run,"cache") else file.path(cfg$execution$cache_parent,folder)
  dir.create(cache,recursive=TRUE,showWarnings=FALSE);writeLines("owned",file.path(cache,".owned_simulation_cache"))
  needed<-vapply(plan$key,function(k)!sim_complete(sim_read(file.path(run,"fits",paste0(k,".rds"))),cfg$execution$retry_failed),logical(1))
  if(!any(needed)){kdr_log("All E final records are available; run report E separately.");return(invisible(run))}
  ids<-unique(plan$kernel_key[needed]);bundles<-lapply(ids,function(k)plan[plan$kernel_key==k,,drop=FALSE])
  for(k in plan$key[needed]) {path<-file.path(run,"events",paste0("fit_",k,".rds"));if(file.exists(path))unlink(path)}
  vars<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
  prior<-Sys.getenv(vars,unset=NA_character_)
  on.exit({for(v in vars)if(is.na(prior[v]))Sys.unsetenv(v) else do.call(Sys.setenv,setNames(list(prior[v]),v))},add=TRUE)
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
  monitor<-sim_monitor_start(cfg,source_root,run,cache);on.exit(sim_monitor_stop(monitor),add=TRUE)
  cl<-parallel::makePSOCKcluster(min(cfg$execution$workers,length(bundles)))
  on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  init<-parallel::clusterCall(cl,function(src,libs,session){.libPaths(libs)
    source(file.path(src,"load.R"),encoding="UTF-8");load_kdrcc_simulation(src)
    options(kdr.resource.session=session);list(pid=Sys.getpid(),R=as.character(getRversion()))
  },source_root,.libPaths(),getOption("kdr.resource.session"))
  kdr_csv(kdr_bind(lapply(seq_along(init),function(i)data.frame(worker=i,pid=init[[i]]$pid,R=init[[i]]$R))),file.path(run,"worker_processes.csv"))
  kdr_log("Active PSOCK workers:",length(cl),"| PIDs:",paste(vapply(init,`[[`,integer(1),"pid"),collapse=", "))
  en<-new.env(parent=emptyenv());en$events<-list();en$times<-numeric();en$run<-run;en$module<-"E"
  en$total<-nrow(plan);en$data_total<-0L;en$initial_fits<-0L;en$initial_data<-0L
  en$started<-en$run_started<-started;en$data_keys<-character()
  sim_scan_events(en,run,FALSE)
  en$initial_fits<-sum(vapply(en$events,function(a)a$kind=="fit",logical(1)))
  sim_dispatch(cl,bundles,"sim_E_worker",cfg,run,cache,en,"fit")
  sim_monitor_stop(monitor);monitor<-NULL
  kdr_log("Module E finished | elapsed",sim_seconds(as.numeric(difftime(Sys.time(),started,units="secs"))))
  kdr_log("Results:",run,"| Run report E separately.");invisible(run)
}
