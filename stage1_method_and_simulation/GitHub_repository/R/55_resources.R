# Resource sampling is observational and never changes the estimator or RNG.
sim_resource_context <- function(run,kind,key,job=NULL,work="",phase="") {
  session<-getOption("kdr.resource.session");if(is.null(session))return(invisible(NULL))
  context<-list(run=run,kind=kind,key=key,job=job,work=work,phase=phase)
  options(kdr.resource.context=context)
  val<-function(n)if(is.null(job)||is.null(job[[n]]))NA_real_ else job[[n]]
  row<-data.frame(pid=Sys.getpid(),role=if(kind=="main")"main" else "worker",kind=kind,key=key,
    phase=phase,N=val("N"),L=val("L"),M=val("M"),kappa=val("kappa"),lambda=val("lambda"),work=work,
    updated=as.numeric(Sys.time()),stringsAsFactors=FALSE)
  folder<-file.path(run,"resources",session,"contexts");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  path<-file.path(folder,paste0(Sys.getpid(),".csv"));tmp<-paste0(path,".tmp")
  kdr_csv(row,tmp)
  if(file.exists(path))unlink(path)
  if(!file.rename(tmp,path))warning("Resource context write failed; estimation continues")
  invisible(NULL)
}
sim_resource_phase <- function(phase) {
  c<-getOption("kdr.resource.context");if(!is.null(c))sim_resource_context(c$run,c$kind,c$key,c$job,c$work,phase)
}
sim_disk_bytes <- function(path) {
  if(!dir.exists(path))return(0)
  f<-list.files(path,recursive=TRUE,full.names=TRUE,all.files=TRUE)
  if(!length(f))return(0)
  a<-file.info(f);sum(a$size[!a$isdir],na.rm=TRUE)
}
sim_resource_checkpoint <- function(run,kind,key,work) {
  folder<-file.path(run,"resource_checkpoints");dir.create(folder,showWarnings=FALSE)
  path<-file.path(folder,paste0(kind,"_",key,".rds"));old<-sim_read(path)
  n<-sim_disk_bytes(work)
  kdr_save(list(kind=kind,key=key,cache_bytes=max(n,if(is.null(old))0 else old$cache_bytes),
    measured_at=Sys.time(),work=work),path);invisible(n)
}
sim_monitor_start <- function(cfg,root,run,cache) {
  if(!isTRUE(cfg$resources$enabled))return(NULL)
  if(.Platform$OS.type!="windows"){
    kdr_log("The process memory sampler requires Windows; memory fields are unavailable on this platform.");return(NULL)}
  session<-paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())
  folder<-file.path(run,"resources",session);dir.create(file.path(folder,"contexts"),recursive=TRUE,showWarnings=FALSE)
  options(kdr.resource.session=session)
  sim_resource_context(run,"main","main",work=cache,phase="Preparing independent data pools")
  exe<-Sys.which("powershell.exe");if(!nzchar(exe))exe<-Sys.which("powershell")
  if(!nzchar(exe)){options(kdr.resource.session=NULL);kdr_log("PowerShell unavailable; memory sampling disabled.");return(NULL)}
  script<-file.path(root,"monitor","sample_resources.ps1")
  args<-c("-NoProfile","-NonInteractive","-ExecutionPolicy","Bypass","-File",shQuote(script),
    "-Folder",shQuote(folder),"-OwnerPid",Sys.getpid(),"-SampleSeconds",cfg$resources$sample_seconds,
    "-DiskSeconds",cfg$resources$disk_seconds)
  launch<-tryCatch({system2(exe,args,wait=FALSE,invisible=TRUE,stdout=FALSE,stderr=FALSE);TRUE},error=function(e)FALSE)
  for(i in 1:50){if(file.exists(file.path(folder,"ready.txt"))||file.exists(file.path(folder,"error.txt")))break;Sys.sleep(.1)}
  if(!launch||!file.exists(file.path(folder,"ready.txt")))kdr_log("Resource sampler startup was not confirmed; available samples will be reported.") else
    kdr_log("Resource sampler started: memory interval",cfg$resources$sample_seconds,"seconds; cache interval",cfg$resources$disk_seconds,"seconds.")
  folder
}
sim_monitor_stop <- function(folder) {
  if(is.null(folder))return(invisible(NULL))
  writeLines("stop",file.path(folder,"stop.txt"))
  for(i in 1:40){if(file.exists(file.path(folder,"finished.txt")))break;Sys.sleep(.1)}
  options(kdr.resource.session=NULL,kdr.resource.context=NULL);invisible(NULL)
}
sim_read_resource_csv <- function(path) {
  if(!file.exists(path))return(data.frame())
  tryCatch(suppressWarnings(read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)),error=function(e)data.frame())
}
sim_resource_tables <- function(run) {
  base<-file.path(run,"resources");folders<-if(dir.exists(base))list.dirs(base,recursive=FALSE,full.names=TRUE) else character()
  rows<-total<-list()
  for(f in folders){
    x<-sim_read_resource_csv(file.path(f,"process_samples.csv"));if(nrow(x)){x$session<-basename(f);rows[[length(rows)+1L]]<-x}
    x<-sim_read_resource_csv(file.path(f,"total_samples.csv"));if(nrow(x)){x$session<-basename(f);total[[length(total)+1L]]<-x}
  }
  d<-kdr_bind(rows);tot<-kdr_bind(total);peak<-list()
  mx<-function(x)if(any(is.finite(x)))max(x[is.finite(x)]) else NA_real_
  if(nrow(d))for(ix in sim_groups(d,c("kind","key"))){a<-d[ix,];
    peak[[length(peak)+1L]]<-data.frame(kind=a$kind[1],key=a$key[1],memory_samples=sum(is.finite(a$working_set_bytes)),
      sampled_working_set_gb=mx(a$working_set_bytes)/1024^3,sampled_private_bytes_gb=mx(a$private_bytes)/1024^3,
      sampled_cache_gb=mx(a$cache_bytes)/1024^3,stringsAsFactors=FALSE)}
  list(peaks=kdr_bind(peak),totals=tot)
}
