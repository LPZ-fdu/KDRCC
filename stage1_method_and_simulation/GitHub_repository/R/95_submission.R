# Submission orchestration; code and writable work directories are separate.
sim_submission_config <- function(settings,source_root,work_root) {
  e<-new.env(parent=.GlobalEnv)
  sys.source(file.path(source_root,"00_parameters.R"),e,keep.source=FALSE)
  cfg<-e$simulation_config(source_root)
  cfg<-sim_merge(cfg,settings$advanced_overrides)
  cfg$seed<-settings$seed;cfg$modules<-settings$modules
  cfg$A<-sim_merge(cfg$A,settings$A_overrides)
  es<-settings$E_settings;ek<-es$kernel_candidates;es$kernel_candidates<-NULL
  cfg$E<-sim_merge(cfg$E,es)
  if(!is.null(ek))cfg$E$kernel_candidates<-ek
  for(m in c("A","B","C","D","E"))if(!is.null(settings$repetitions[[m]]))
    cfg[[m]]$repetitions<-settings$repetitions[[m]]
  if(!is.null(settings$D_correlation_levels))cfg$D$correlation_levels<-settings$D_correlation_levels
  cfg$D$baseline_B_run<-settings$D_baseline_run
  cfg$C$null_counts<-settings$null_counts
  cfg$configuration$fixed_file<-file.path("config","fixed_parameters.R")
  cfg$configuration$overrides<-list(
    B=sim_merge(settings$common_fixed,settings$B_overrides),
    C=sim_merge(settings$common_fixed,settings$C_overrides),
    D=sim_merge(settings$common_fixed,settings$D_overrides),
    E=sim_merge(settings$common_fixed,settings$E_overrides))
  cfg$execution$workers<-settings$workers
  cfg$execution$resume<-settings$resume
  cfg$execution$retry_failed<-settings$retry_failed
  cfg$execution$clear_stale_lock<-settings$clear_stale_lock
  cfg$report<-sim_merge(cfg$report,settings$report_overrides)
  cfg
}

execute_submission <- function(settings,source_root,work_root,submission_file=NULL) {
  source_root<-normalizePath(path.expand(source_root),winslash="/",mustWork=TRUE)
  work_root<-path.expand(work_root)
  dir.create(work_root,recursive=TRUE,showWarnings=FALSE)
  work_root<-normalizePath(work_root,winslash="/",mustWork=TRUE)
  compare<-function(x)if(.Platform$OS.type=="windows")tolower(x) else x
  a<-compare(source_root);b<-compare(work_root)
  if(a==b||startsWith(b,paste0(a,"/"))||startsWith(a,paste0(b,"/")))
    stop("Use separate, non-nested source and work directories.")
  allowed<-c("install","check","plan","run","report","freeze","contributions")
  if(length(settings$action)!=1L||!settings$action%in%allowed)stop("Unknown submission action.")
  if(!length(settings$modules)||anyDuplicated(settings$modules)||!all(settings$modules%in%c("A","B","C","D","E")))
    stop("modules must contain unique entries from A, B, C, D, E.")
  if(length(settings$workers)!=1L||!is.finite(settings$workers)||settings$workers<1||settings$workers!=floor(settings$workers))
    stop("workers must be a positive integer.")
  options(kdr.root=source_root)
  cfg<-sim_submission_config(settings,source_root,work_root)
  for(d in c("results","config","submissions"))dir.create(file.path(work_root,d),recursive=TRUE,showWarnings=FALSE)
  stamp<-paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid(),"_",settings$action)
  record_dir<-file.path(work_root,"submissions",stamp)
  if(dir.exists(record_dir))record_dir<-paste0(record_dir,"_",length(list.dirs(file.path(work_root,"submissions"),recursive=FALSE))+1L)
  dir.create(record_dir,recursive=TRUE)
  kdr_save(list(settings=settings,source_root=source_root,work_root=work_root,created_at=Sys.time()),
    file.path(record_dir,"submission_settings.rds"))
  dput(settings,file=file.path(record_dir,"submission_settings.R"))
  kdr_save(cfg,file.path(record_dir,"effective_input_config.rds"))
  if(!is.null(submission_file))file.copy(submission_file,file.path(record_dir,"submitted_script.R"))
  fs<-list.files(source_root,pattern="[.](R|md|txt|ps1|csv)$|^DESCRIPTION$",recursive=TRUE,full.names=TRUE)
  fs<-fs[!grepl("/(results|validation)/",gsub("\\\\","/",fs))]
  hash<-tools::md5sum(fs)
  kdr_csv(data.frame(file=substring(fs,nchar(source_root)+2L),md5=unname(hash)),file.path(record_dir,"source_md5.csv"))
  capture.output(sessionInfo(),file=file.path(record_dir,"sessionInfo.txt"))
  log<-file(file.path(record_dir,"submission.log"),open="at",encoding="UTF-8")
  sink(log,split=TRUE);on.exit({sink();close(log)},add=TRUE)
  kdr_log("Submission action:",settings$action,"| modules:",paste(settings$modules,collapse=", "))
  kdr_log("Read-only source:",source_root)
  kdr_log("Writable work directory:",work_root)
  kdr_log("Requested PSOCK workers:",settings$workers)
  started<-Sys.time()
  answer<-tryCatch({
    if(settings$action=="install") {
      source(file.path(source_root,"00_install_dependencies.R"),encoding="UTF-8",local=new.env(parent=.GlobalEnv))
      result<-record_dir
    } else if(settings$action=="check") {
      sim_check_compatibility(source_root,work_root,parallel_check=TRUE)
      source(file.path(source_root,"tests","unit_checks.R"),encoding="UTF-8")
      result<-run_unit_checks(source_root,output_root=work_root)
      source(file.path(source_root,"tests","module_D_checks.R"),encoding="UTF-8")
      run_module_D_checks(source_root,work_root)
      source(file.path(source_root,"tests","module_E_checks.R"),encoding="UTF-8")
      run_module_E_checks(source_root,work_root)
      source(file.path(source_root,"tests","method_contract_checks.R"),encoding="UTF-8")
      run_method_contract_checks(source_root,work_root)
      source(file.path(source_root,"tests","contribution_checks.R"),encoding="UTF-8")
      run_contribution_checks(source_root,work_root)
    } else if(settings$action%in%c("plan","run")) {
      # Resolve every requested module before starting expensive computations.
      plans<-lapply(settings$modules,function(m){
        resolved<-if(m=="A")NULL else sim_resolve_fixed(cfg,work_root,m)
        p<-sim_plan(cfg,m,resolved);sim_validate(cfg,p$tasks)
        if(m=="D") {
          reference<-sim_D_source_B(cfg,work_root,resolved)
          kdr_log("Existing B reference:",reference$path,"| saved fits",reference$fit_records)
          kdr_log("Module D fixed configuration:");print(as.data.frame(resolved$parameters),row.names=FALSE)
          kdr_log("Module D overall-correlation settings:");print(sim_D_levels(cfg),row.names=FALSE)
        }
        if(m=="E") {
          kdr_log("E fixed comparison configuration:");print(as.data.frame(resolved$parameters),row.names=FALSE)
          kdr_log("E scenarios:");print(p$design$scenarios,row.names=FALSE)
          kdr_log("E candidate configurations:");print(p$design$configurations,row.names=FALSE)
          kdr_log("E validation fits:",nrow(p$tasks)/2*p$cv_fits_per_repetition,
            "| independent paired repetitions:",nrow(p$tasks)/2,
            "| independent learning sets:",resolved$parameters$L,
            "| configuration pool sets:",cfg$E$tuning_sets,"| final estimation sets:",resolved$parameters$N)
        }
        kdr_csv(p$tasks,file.path(record_dir,paste0("plan_",m,".csv")))
        if(!is.null(resolved))dput(resolved,file=file.path(record_dir,paste0("resolved_",m,".R")))
        data.frame(module=m,planned_fits=nrow(p$tasks),repetitions=cfg[[m]]$repetitions,
          configurations=nrow(p$design$configurations),null_count_settings=length(unique(p$tasks$null_count)),
          distinct_learning_tasks=length(unique(p$tasks$kernel_key))*
            if(m=="E")length(unique(p$design$configurations$kernel_id)) else 1L)
      })
      overview<-kdr_bind(plans);print(overview,row.names=FALSE)
      kdr_csv(overview,file.path(record_dir,"plan_overview.csv"))
      kdr_log("Planned fit records:",sum(overview$planned_fits),
        if("C"%in%settings$modules)"| compatible B/G1 reuse may reduce new fits in C." else "")
      if(settings$action=="run")sim_check_compatibility(source_root,work_root,parallel_check=FALSE)
      result<-if(settings$action=="plan")overview else run_simulation(cfg,work_root)
    } else if(settings$action%in%c("report","contributions")) {
      runs<-settings$report_run_directories
      if(is.null(runs)) {
        runs<-lapply(settings$modules,function(m)sim_pointer(work_root,m))
        if(any(vapply(runs,is.null,logical(1))))stop("A requested module run was not found. Set report_run_directories to its complete path.")
      }
      if(!length(runs))stop("No requested module results were found in the work directory.")
      result<-if(settings$action=="report")lapply(runs,prepare_results,report_overrides=settings$report_overrides) else
        lapply(runs,prepare_contributions,workers=settings$workers,resume=settings$resume,source_root=source_root)
    } else if(settings$action=="freeze") {
      run<-settings$freeze$source_run
      if(is.null(run))run<-sim_pointer(work_root,"A")
      result<-freeze_configuration(work_root,run,settings$freeze$config_id,settings$freeze$manual_overrides)
    }
    result
  },error=function(e){
    kdr_log("Submission stopped:",conditionMessage(e))
    kdr_save(list(ok=FALSE,error=conditionMessage(e),finished_at=Sys.time()),file.path(record_dir,"outcome.rds"))
    stop(e)
  })
  kdr_save(list(ok=TRUE,result=answer,finished_at=Sys.time()),file.path(record_dir,"outcome.rds"))
  kdr_log("Submission finished | elapsed",sim_seconds(as.numeric(difftime(Sys.time(),started,units="secs"))))
  if(settings$action=="run")kdr_log("Run reporting separately with action='report'; no figures were generated during fitting.")
  invisible(answer)
}
