# Contributions evaluated from saved fitted functions at a specified background.
sim_contribution_files <- function(run,a) {
  j<-a$job
  kernel<-if(j$method=="KDR")file.path(run,"kernels",paste0(a$kernel_key,".rds")) else
    file.path(run,"comparators",paste0(j$kernel_key,"_",j$method,".rds"))
  evaluation<-if(identical(j$module,"E"))file.path(run,"evaluation",paste0(j$scenario_id,".rds")) else
    sim_evaluation_path(run,j$p,j$correlation_id)
  c(fit=file.path(run,"fits",paste0(j$key,".rds")),kernel=kernel,evaluation=evaluation,
    config=file.path(run,"config.rds"))
}

sim_contribution_signature <- function(files,a,source_root=getOption("kdr.root")) {
  if(any(!file.exists(files)))stop("Missing saved contribution input: ",paste(files[!file.exists(files)],collapse=", "))
  code<-file.path(source_root,"R",c("42_attribution.R","45_interpretation.R"))
  inputs<-files[names(files)!="fit"]
  list(rule=sim_contribution_rule(),fit=object_id(list(job=a$job,beta=a$fit$beta)),
    inputs=setNames(unname(tools::md5sum(inputs)),names(inputs)),
    code=setNames(unname(tools::md5sum(code)),basename(code)))
}

sim_contribution_task <- function(key,run,resume=TRUE) {
  started<-proc.time()[3]
  answer<-tryCatch({
    a<-readRDS(file.path(run,"fits",paste0(key,".rds")))
    if(!isTRUE(a$ok))return(list(key=key,status="unavailable",error="Fitted function unavailable",seconds=0))
    files<-sim_contribution_files(run,a);signature<-sim_contribution_signature(files,a)
    if(resume&&identical(a$interpretation$contribution_scope,sim_contribution_rule())&&
       identical(a$interpretation$contribution_input_signature,signature))
      return(list(key=key,status="cached",error="",seconds=unname(proc.time()[3]-started)))
    s<-readRDS(files["kernel"]);ev<-readRDS(files["evaluation"]);cfg<-readRDS(files["config"])
    value<-sim_active_contributions(a$fit,s,ev,cfg)
    for(n in names(value))a$interpretation[[n]]<-value[[n]]
    a$interpretation$contribution_input_signature<-signature
    a$interpretation$contribution_evaluated_at<-Sys.time()
    kdr_save(a,files["fit"])
    list(key=key,status="completed",error="",seconds=unname(proc.time()[3]-started))
  },error=function(e)list(key=key,status="failed",error=conditionMessage(e),seconds=unname(proc.time()[3]-started)))
  answer
}

sim_contribution_attach <- function(run,records) {
  for(i in seq_along(records)) {
    a<-records[[i]]
    if(!isTRUE(a$interpretation_ok))next
    if(identical(a$interpretation$contribution_scope,sim_contribution_rule()))next
    # With eight candidate coordinates the full and active games coincide.
    if(a$job$p==length(exposure_names)) {
      z<-a$interpretation$contributions
      if(nrow(z)!=length(exposure_names)||!setequal(z$variable,exposure_names))
        stop("Saved active-exposure contribution rows are incomplete: ",a$job$key)
      a$interpretation$contributions<-z[match(exposure_names,z$variable),,drop=FALSE]
      a$interpretation$contribution_scope<-sim_contribution_rule()
    } else {
      stop("Active-exposure contributions are required for ",a$job$key,
        ". Run: Rscript --vanilla 05_submit.R contributions ",a$job$module)
    }
    records[[i]]<-a
  }
  records
}

prepare_contributions <- function(run,workers=32L,resume=TRUE,source_root=getOption("kdr.root")) {
  run<-normalizePath(run,winslash="/",mustWork=TRUE)
  if(dir.exists(file.path(run,"run_lock"))||dir.exists(file.path(run,"report_lock")))
    stop("Wait for fitting and reporting to finish before preparing contributions.")
  lock<-file.path(run,"contribution_lock")
  if(!dir.create(lock,showWarnings=FALSE))stop("Contribution lock exists. Confirm no contribution process is active before removing it.")
  on.exit(sim_unlock(lock,run),add=TRUE)
  plan<-readRDS(file.path(run,"plan.rds"))$tasks
  present<-file.exists(file.path(run,"fits",paste0(plan$key,".rds")))
  keys<-plan$key[present];n<-length(keys)
  if(!n)stop("No saved fitted records were found.")
  if(length(workers)!=1L||!is.finite(workers)||workers<1||workers!=as.integer(workers))stop("Invalid worker count.")
  workers<-min(as.integer(workers),n);started<-Sys.time();answers<-vector("list",n)
  kdr_log("Active-exposure contribution evaluation |",basename(run),"| saved records",n,"| pending records",sum(!present))
  finish<-function(ans,i,completed) {
    answers[[i]]<<-as.data.frame(ans,stringsAsFactors=FALSE)
    elapsed<-as.numeric(difftime(Sys.time(),started,units="secs"))
    kdr_csv(kdr_bind(answers),file.path(run,"contribution_processing.csv"))
    kdr_log("Contributions",sim_bar(completed,n),"|",ans$key,"|",ans$status,
      "| elapsed",sim_seconds(elapsed),"| remaining",sim_seconds(elapsed/completed*(n-completed)))
  }
  if(workers==1L) {
    for(i in seq_len(n))finish(sim_contribution_task(keys[i],run,resume),i,i)
  } else {
    vars<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    prior<-Sys.getenv(vars,unset=NA_character_)
    on.exit({for(v in vars)if(is.na(prior[v]))Sys.unsetenv(v) else do.call(Sys.setenv,setNames(list(prior[v]),v))},add=TRUE)
    for(v in vars)do.call(Sys.setenv,setNames(list("1"),v))
    cl<-parallel::makePSOCKcluster(workers)
    on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
    pids<-parallel::clusterCall(cl,function(src,libs){.libPaths(libs)
      source(file.path(src,"load.R"),encoding="UTF-8");load_kdrcc_simulation(src);Sys.getpid()
    },source_root,.libPaths())
    kdr_log("Active PSOCK workers:",workers,"| PIDs:",paste(unlist(pids),collapse=", "))
    send<-getFromNamespace("sendCall","parallel");receive<-getFromNamespace("recvOneResult","parallel")
    dispatch<-function(node,i)send(cl[[node]],sim_contribution_task,list(keys[i],run,resume),tag=i)
    for(i in seq_len(workers))dispatch(i,i)
    next_index<-workers+1L
    for(completed in seq_len(n)) {
      a<-receive(cl)
      if(inherits(a$value,"try-error"))stop("Contribution worker error: ",a$value)
      finish(a$value,a$tag,completed)
      if(next_index<=n){dispatch(a$node,next_index);next_index<-next_index+1L}
    }
  }
  status<-kdr_bind(answers)
  kdr_log("Contribution records saved in:",file.path(run,"fits"),"| elapsed",
    sim_seconds(as.numeric(difftime(Sys.time(),started,units="secs"))))
  if(any(status$status=="failed"))stop("Some contribution tasks failed; inspect contribution_processing.csv and rerun the contribution action.")
  # The archived evaluation source describes the calculations in the fitted records.
  files<-c("R/45_interpretation.R","R/65_execution.R","R/80_summary.R","R/81_interpretation_summary.R",
    "R/82_contribution_processing.R","R/88_interpretation_publication.R","R/90_report.R",
    "R/92_E_publication.R","R/95_submission.R","05_submit.R","06_prepare_contributions.R",
    "README.md","运行说明.txt","docs/METHOD_SPECIFICATION.md","docs/SIMULATION_PROTOCOL.md",
    "docs/MODULE_E.md","docs/OUTPUTS.md","docs/CONTRIBUTIONS.md","docs/SERVER_RUN.md",
    "tests/contribution_checks.R","tests/saved_contribution_checks.R")
  snapshot<-file.path(run,"source_snapshot")
  if(dir.exists(snapshot))for(f in files) {
    src<-file.path(source_root,f);dest<-file.path(snapshot,f)
    if(file.exists(src)){dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
      if(!file.copy(src,dest,overwrite=TRUE))stop("Cannot save evaluation source: ",f)}
  }
  hashes<-tools::md5sum(file.path(source_root,files[file.exists(file.path(source_root,files))]))
  names(hashes)<-substring(names(hashes),nchar(source_root)+2L)
  kdr_save(list(scope=sim_contribution_rule(),active_exposures=exposure_names,
    background="Marginal-median reference",completed_at=Sys.time(),source_md5=hashes,
    available_records=n,pending_records=sum(!present)),file.path(run,"contribution_processing.rds"))
  invisible(list(run=run,completed=sum(status$status=="completed"),cached=sum(status$status=="cached"),
    unavailable=sum(status$status=="unavailable"),pending=sum(!present)))
}
