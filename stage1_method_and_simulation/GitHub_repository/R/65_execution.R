sim_science <- function(cfg) {
  x<-cfg;x$modules<-x$configuration<-x$execution<-x$resources<-x$report<-NULL
  x$block_groups<-cfg$execution$block_groups;x
}
sim_compatibility <- function(cfg,parameters,code_hash) {
  x<-sim_science(cfg);x$A<-x$B<-x$C<-NULL
  object_id(list(science=x,parameters=parameters,code=code_hash))
}
sim_snapshot <- function(root,run) {
  dest<-file.path(run,"source_snapshot");if(dir.exists(dest))return(invisible(NULL))
  dir.create(dest,recursive=TRUE)
  for(d in c("R","config","docs","tests","monitor")) {
    dir.create(file.path(dest,d),showWarnings=FALSE)
    fs<-list.files(file.path(root,d),full.names=TRUE)
    if(length(fs)&&!all(file.copy(fs,file.path(dest,d),recursive=TRUE)))stop("Source snapshot failed")
  }
  fs<-list.files(root,pattern="[.](R|md|txt)$|^DESCRIPTION$|^LICENSE$",full.names=TRUE)
  if(length(fs))file.copy(fs,dest)
}
sim_unlock <- function(path,run) {
  if(!dir.exists(path))return(invisible(NULL))
  a<-normalizePath(path,winslash="/",mustWork=TRUE);b<-normalizePath(run,winslash="/",mustWork=TRUE)
  if(dirname(a)!=b || !basename(a)%in%c("run_lock","report_lock","contribution_lock"))stop("Invalid lock path")
  unlink(a,recursive=TRUE);invisible(NULL)
}
sim_pointer <- function(root,module) {
  f<-file.path(root,"results",paste0("LATEST_",module,".txt"))
  if(!file.exists(f))return(NULL)
  run<-file.path(root,"results",basename(readLines(f,warn=FALSE)[1]))
  if(dir.exists(run))run else NULL
}
sim_reuse_B <- function(root,run,cfg,plan,manifest) {
  if(!cfg$C$reuse_B_baseline || !any(plan$null_count==0L))return(0L)
  source<-sim_pointer(root,"B");if(is.null(source))return(0L)
  old<-sim_read(file.path(source,"manifest.rds"))
  if(is.null(old)||!identical(old$compatibility_id,manifest$compatibility_id)) {
    kdr_log("Module B reuse skipped: effective scientific configurations differ.");return(0L)
  }
  n<-0L
  for(i in which(plan$null_count==0L)) {
    job<-as.list(plan[i,,drop=FALSE]);dest<-file.path(run,"fits",paste0(job$key,".rds"))
    if(!is.null(sim_read(dest)))next
    a<-sim_read(file.path(source,"fits",paste0(job$key,".rds")))
    if(is.null(a)||!isTRUE(a$ok)||!isTRUE(a$interval_ok)||!isTRUE(a$selection$valid))next
    a$job<-job;a$reused<-TRUE;a$reused_from<-file.path(basename(source),"fits",paste0(job$key,".rds"))
    a$elapsed_seconds<-0;sim_save_fit(run,a);n<-n+1L
    for(d in c("kernels","kernel_records")) {
      src<-file.path(source,d,paste0(job$kernel_key,".rds"));to<-file.path(run,d,basename(src))
      if(file.exists(src)&&!file.exists(to))file.copy(src,to)
    }
  }
  if(n)kdr_log("Reused",n,"identical module B G1 fits for the zero-null baseline; no independent repetition is claimed.")
  n
}

sim_scan_events <- function(en,run,notify=TRUE,force=FALSE) {
  directory<-file.path(run,"events");stamp<-as.numeric(file.info(directory)$mtime);now<-Sys.time()
  if(!force && !is.null(en$event_stamp) && identical(stamp,en$event_stamp) &&
     as.numeric(difftime(now,en$last_scan,units="secs"))<30)return(invisible(NULL))
  en$event_stamp<-stamp;en$last_scan<-now
  fs<-list.files(directory,pattern="[.]rds$",full.names=TRUE)
  if(!length(fs))return(invisible(NULL))
  times<-as.numeric(file.info(fs)$mtime)
  for(i in seq_along(fs)) {
    id<-basename(fs[i]);old<-en$times[id]
    if(length(old)&&!is.na(old)&&old==times[i])next
    a<-sim_read(fs[i]);if(is.null(a))next
    en$times[id]<-times[i];en$events[[id]]<-a
    if(notify && a$kind=="fit") {
      kdr_log(if(a$ok)"Completed" else "Failed",a$key,"| fit elapsed",sim_seconds(a$elapsed_seconds))
      sim_progress(en,list(),FALSE)
    }
  }
}
sim_progress <- function(en,active,detail=TRUE,kind="fit") {
  a<-en$events;selected<-a[vapply(a,function(x)x$kind==kind,logical(1))]
  if(kind=="data")selected<-selected[vapply(selected,function(x)x$key%in%en$data_keys,logical(1))]
  done<-length(selected);total<-if(kind=="fit")en$total else en$data_total
  elapsed<-as.numeric(difftime(Sys.time(),en$run_started,units="secs"))
  phase_elapsed<-as.numeric(difftime(Sys.time(),en$started,units="secs"))
  # ETA uses observed completed work in the current phase. It is not a time limit.
  fresh<-max(0,done-if(kind=="fit")en$initial_fits else en$initial_data)
  eta<-if(fresh>0)phase_elapsed/fresh*max(0,total-done) else NA_real_
  if(done>=total)eta<-0
  good<-sum(vapply(selected,function(x)isTRUE(x$ok),logical(1)))
  kdr_log(if(kind=="fit")paste("Module",en$module,"overall") else "Data preparation",
    sim_bar(done,total),"| success",good,"failed",done-good,
    "| elapsed",sim_seconds(elapsed),"| estimated remaining",sim_seconds(eta))
  if(detail)for(w in names(active)) {
    st<-sim_read(file.path(en$run,"status",paste0(active[[w]]$key,".rds")))
    if(is.null(st))kdr_log("Worker",w,active[[w]]$key,"| starting") else
      kdr_log("Worker",w,st$key,"|",sim_bar(st$stage_number,st$stage_total,8),st$stage,
        "|",sim_bar(st$i,st$n,8),st$detail,
        if(!is.null(st$validation_total))paste("| validation",st$validation_done,"/",st$validation_total) else "")
  }
  invisible(NULL)
}
sim_dispatch <- function(cl,bundles,fun,cfg,run,cache,en,kind) {
  if(!length(bundles))return(invisible(NULL))
  active<-list();next_job<-1L
  send<-function(w) {
    if(next_job>length(bundles))return(FALSE)
    b<-bundles[[next_job]];next_job<<-next_job+1L
    key<-if(kind=="data")sim_data_key(b) else b$kernel_key[1]
    parallel:::sendCall(cl[[w]],function(b,fn,cfg,run,cache,slot)get(fn)(b,cfg,run,cache,slot),
      list(b,fun,cfg,run,cache,w),tag=next_job-1L)
    active[[as.character(w)]]<<-list(key=key,bundle=b);TRUE
  }
  for(w in seq_len(min(length(cl),length(bundles))))send(w)
  last<-Sys.time();sim_progress(en,active,TRUE,kind)
  while(length(active)) {
    sim_scan_events(en,run,TRUE)
    slots<-as.integer(names(active))
    ready<-socketSelect(lapply(cl[slots],function(x)x$con),timeout=cfg$execution$poll_seconds)
    if(any(ready)) {
      answer<-tryCatch(parallel:::recvOneResult(cl),error=function(e)stop("Worker connection interrupted; checkpoints can be resumed: ",conditionMessage(e)))
      w<-answer$node;current<-active[[as.character(w)]];value<-answer$value
      if(is.null(current))stop("Unregistered worker result")
      if(inherits(value,"try-error")||!is.list(value)||!isTRUE(value$ok)) {
        err<-if(is.list(value)&&!is.null(value$error))value$error else paste(value,collapse=" ")
        if(kind=="data")stop("Data preparation failed: ",err)
        for(i in seq_len(nrow(current$bundle))) {
          job<-as.list(current$bundle[i,,drop=FALSE])
          if(is.null(sim_read(file.path(run,"fits",paste0(job$key,".rds")))))sim_failed_fit(run,job,err,"worker")
        }
      }
      active[[as.character(w)]]<-NULL;send(w);sim_scan_events(en,run,TRUE)
    }
    if(as.numeric(difftime(Sys.time(),last,units="secs"))>=cfg$execution$progress_seconds) {
      sim_progress(en,active,TRUE,kind);last<-Sys.time()
    }
  }
  sim_scan_events(en,run,TRUE,force=TRUE);sim_progress(en,list(),FALSE,kind);invisible(NULL)
}

run_simulation_module <- function(cfg,root,module) {
  if(module=="E")return(run_module_E(cfg,root))
  module_started<-Sys.time()
  # root is the writable work directory; source_root may be read-only.
  source_root<-getOption("kdr.root",root)
  if(!file.exists(file.path(source_root,"load.R")))stop("Simulation source directory is unavailable.")
  resolved<-if(module=="A")NULL else sim_resolve_fixed(cfg,root,module)
  if(module=="D")cfg$D$source_B<-sim_D_source_B(cfg,root,resolved)
  pp<-sim_plan(cfg,module,resolved);plan<-pp$tasks;sim_validate(cfg,plan)
  cfg$batch<-if(module=="A")"configuration" else if(module=="D")"correlation" else "validation"
  code_files<-list.files(file.path(source_root,"R"),pattern="[.]R$",full.names=TRUE)
  code_files<-code_files[as.integer(substr(basename(code_files),1,2))<80L]
  hashes<-tools::md5sum(code_files);names(hashes)<-basename(code_files)
  id<-substr(object_id(list(science=sim_science(cfg),module=module,plan=plan,
    parameters=resolved$parameters,code=hashes)),1,12)
  folder<-paste0("module_",module,"_",id)
  if(!cfg$execution$resume)folder<-paste0(folder,"_",format(Sys.time(),"%Y%m%d_%H%M%S"))
  run<-file.path(root,"results",folder)
  for(d in c("fits","kernels","kernel_records","comparators","auxiliaries","points","events","status","logs","evaluation"))
    dir.create(file.path(run,d),recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(run,"run_lock")
  if(dir.exists(lock)&&cfg$execution$clear_stale_lock)sim_unlock(lock,run)
  if(!dir.create(lock,showWarnings=FALSE))stop("Run lock exists. Confirm no process is active before enabling clear_stale_lock.")
  on.exit(sim_unlock(lock,run),add=TRUE)
  writeLines(as.character(Sys.getpid()),file.path(lock,"owner_pid.txt"))
  logfile<-file(file.path(run,"logs",paste0("run_",format(Sys.time(),"%Y%m%d_%H%M%S"),".txt")),"at",encoding="UTF-8")
  sink(logfile,split=TRUE);on.exit({sink();close(logfile)},add=TRUE)
  if(module=="D") {
    sim_D_import_B(cfg,run)
    cfg<-sim_D_prepare_generator(cfg,run)
    pp<-sim_plan(cfg,module,resolved);plan<-pp$tasks
  }
  manifest<-list(module=module,run_id=id,created_at=Sys.time(),effective_configuration=resolved,
    compatibility_id=if(module=="A")NA_character_ else sim_compatibility(cfg,resolved$parameters,hashes),
    code_hash=hashes,planned_fits=nrow(plan),batch=cfg$batch)
  kdr_save(cfg,file.path(run,"config.rds"));kdr_save(pp,file.path(run,"plan.rds"));kdr_save(manifest,file.path(run,"manifest.rds"))
  kdr_csv(plan,file.path(run,"task_manifest.csv"));kdr_csv(pp$design$configurations,file.path(run,"parameter_configurations.csv"))
  if(nrow(pp$design$membership))kdr_csv(pp$design$membership,file.path(run,"configuration_views.csv"))
  if(!is.null(resolved))dput(resolved,file=file.path(run,"effective_configuration.R"))
  sim_snapshot(source_root,run);capture.output(sessionInfo(),file=file.path(run,"sessionInfo.txt"))
  kdr_csv(data.frame(package="CompQuadForm",version=as.character(packageVersion("CompQuadForm"))),file.path(run,"package_versions.csv"))
  writeLines(folder,file.path(root,"results",paste0("LATEST_",module,".txt")))
  cache<-if(is.null(cfg$execution$cache_parent))file.path(run,"cache") else file.path(cfg$execution$cache_parent,folder)
  dir.create(cache,recursive=TRUE,showWarnings=FALSE);writeLines("owned",file.path(cache,".owned_simulation_cache"))
  kdr_log("Module",module,"| planned model fits",nrow(plan),"| requested workers",cfg$execution$workers)
  kdr_log("Batch:",cfg$batch,"| generating mechanism and semantic random seeds are recorded.")
  if(module=="D") {
    kdr_log("Fixed configuration:");print(as.data.frame(resolved$parameters),row.names=FALSE)
  }
  if(module=="C")sim_reuse_B(root,run,cfg,plan,manifest)
  if(module=="D")sim_prepare_D_evaluations(cfg,run) else for(p in sort(unique(plan$p))) {
    kdr_log("Preparing fixed independent evaluation points for",p,"exposures")
    sim_evaluation(cfg,p,file.path(run,"evaluation",paste0("p",p,".rds")))
  }
  needed<-vapply(plan$key,function(k)!sim_complete(sim_read(file.path(run,"fits",paste0(k,".rds"))),cfg$execution$retry_failed),logical(1))
  if(!any(needed)){kdr_log("All planned records are available. Run 02_prepare_results.R separately.");return(invisible(run))}
  todo<-plan[needed,,drop=FALSE]
  # Retried records must not be counted as completed in the new progress tally.
  for(k in todo$key) {
    event_path<-file.path(run,"events",paste0("fit_",k,".rds"))
    if(file.exists(event_path))unlink(event_path)
  }
  data_jobs<-if(module=="D")sim_D_data_jobs(todo) else
    lapply(sort(unique(todo$replicate)),function(r)list(replicate=r,L=max(todo$L[todo$replicate==r]),N=max(todo$N[todo$replicate==r])))
  monitor<-sim_monitor_start(cfg,source_root,run,cache);on.exit(sim_monitor_stop(monitor),add=TRUE)
  # Avoid nested BLAS/OpenMP oversubscription in PSOCK workers.
  thread_vars<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
  previous_threads<-Sys.getenv(thread_vars,unset=NA_character_)
  on.exit({for(v in thread_vars)if(is.na(previous_threads[v]))Sys.unsetenv(v) else do.call(Sys.setenv,setNames(list(previous_threads[v]),v))},add=TRUE)
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
  cl<-parallel::makePSOCKcluster(min(cfg$execution$workers,max(length(data_jobs),length(unique(todo$kernel_key)))))
  on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  libpaths<-.libPaths();resource_session<-getOption("kdr.resource.session")
  init<-parallel::clusterCall(cl,function(root,libs,resource_session){
    .libPaths(libs);source(file.path(root,"load.R"),encoding="UTF-8");load_kdrcc_simulation(root)
    options(kdr.resource.session=resource_session);list(pid=Sys.getpid(),R=as.character(getRversion()))
  },source_root,libpaths,resource_session)
  kdr_csv(kdr_bind(lapply(seq_along(init),function(i)data.frame(worker=i,pid=init[[i]]$pid,R=init[[i]]$R))),file.path(run,"worker_processes.csv"))
  kdr_log("Active PSOCK workers:",length(cl),"| PIDs:",paste(vapply(init,`[[`,integer(1),"pid"),collapse=", "))
  en<-new.env(parent=emptyenv());en$events<-list();en$times<-numeric();en$run<-run;en$module<-module
  en$total<-nrow(plan);en$data_total<-length(data_jobs);en$initial_fits<-0;en$initial_data<-0;en$started<-Sys.time()
  en$run_started<-module_started
  en$data_keys<-vapply(data_jobs,sim_data_key,character(1))
  sim_scan_events(en,run,FALSE)
  en$initial_fits<-sum(vapply(en$events,function(x)x$kind=="fit",logical(1)))
  en$initial_data<-0
  sim_dispatch(cl,data_jobs,"sim_prepare_worker",cfg,run,cache,en,"data")
  groups<-split(seq_len(nrow(todo)),todo$kernel_key)
  bundles<-lapply(groups,function(ix)todo[ix,,drop=FALSE])
  # Scheduling changes execution order only; it does not enter random seeds.
  cost<-vapply(bundles,function(d)d$L[1]*d$M[1]^2*d$p[1]+sum(d$N*d$M^2),numeric(1))
  bundles<-bundles[order(cost,decreasing=TRUE)]
  en$started<-Sys.time();sim_dispatch(cl,bundles,"sim_worker",cfg,run,cache,en,"fit")
  sim_monitor_stop(monitor);monitor<-NULL
  kdr_log("Module",module,"finished. Completed and failed records are retained separately in the fit metadata.")
  kdr_log("Results:",run)
  kdr_log("No figures or analysis tables were generated. Run 02_prepare_results.R when ready.")
  invisible(run)
}
run_simulation <- function(cfg,root,modules=cfg$modules) {
  root<-normalizePath(root,winslash="/",mustWork=TRUE)
  if(!all(modules%in%c("A","B","C","D","E")))stop("modules must contain A, B, C, D and/or E")
  out<-setNames(vector("list",length(modules)),modules)
  for(m in modules)out[[m]]<-run_simulation_module(cfg,root,m)
  invisible(out)
}
