exposure_names <- c("PM25","PMcoarse","O3","NO2","SO2","CO","temperature","humidity")
kdr_names <- function(p, profile="eight") c(exposure_names, if(p>8L) paste0("noise",seq_len(p-8L)))
kdr_seed <- function(cfg, ...) stage_seed(cfg$seed, cfg$schema, cfg$batch, ...)
kdr_log <- function(...) {
  cat(sprintf("[%s] ",format(Sys.time(),"%H:%M:%S")),paste(...,collapse=" "),"\n",sep="")
  flush.console()
}
sim_read <- function(path) if(file.exists(path)) suppressWarnings(tryCatch(readRDS(path),error=function(e)NULL)) else NULL
sim_seconds <- function(x) {
  if(length(x)!=1L || !is.finite(x) || x<0) return("estimating")
  x<-round(x);sprintf("%02d:%02d:%02d",x%/%3600,(x%%3600)%/%60,x%%60)
}
sim_bar <- function(i,n,width=20L) {
  k<-if(n>0)max(0,min(width,floor(width*i/n))) else 0L
  paste0("[",paste(rep("=",k),collapse=""),paste(rep(".",width-k),collapse=""),"] ",i,"/",n)
}
sim_merge <- function(base, override) {
  for(n in names(override)) base[[n]] <- if(is.list(base[[n]]) && is.list(override[[n]]))
    sim_merge(base[[n]],override[[n]]) else override[[n]]
  base
}
sim_safe_remove <- function(path, root) {
  if(!dir.exists(path))return(invisible(NULL))
  a<-normalizePath(path,winslash="/",mustWork=TRUE);b<-normalizePath(root,winslash="/",mustWork=TRUE)
  if(!file.exists(file.path(b,".owned_simulation_cache")) || !startsWith(tolower(a),paste0(tolower(b),"/work/")))
    stop("Refusing to remove a directory outside the owned work cache.")
  unlink(a,recursive=TRUE);invisible(NULL)
}
sim_status <- function(detail=NULL,i=NULL,n=NULL,force=FALSE) {
  e<-getOption("kdr.worker");if(is.null(e))return(invisible(NULL))
  if(!is.null(detail))e$detail<-detail
  if(!is.null(i)){e$i<-i;e$n<-n}
  now<-Sys.time()
  if(force || as.numeric(difftime(now,e$updated,units="secs"))>=e$refresh) {
    e$updated<-now
    kdr_save(list(key=e$key,stage=e$stage,stage_number=e$number,stage_total=6L,
      detail=e$detail,i=e$i,n=e$n,updated=now,started=e$started,pid=Sys.getpid(),
      validation_done=e$validation_done,validation_total=e$validation_total),e$path)
  }
  invisible(NULL)
}
kdr_bar <- function(label,i,n,detail="") sim_status(paste(label,detail),i,n,TRUE)
kdr_heartbeat <- function(label) sim_status(label)
sim_major <- function(i,label) {
  e<-getOption("kdr.worker");if(!is.null(e)){e$stage<-label;e$number<-i}
  sim_status(label,0,1,TRUE);sim_resource_phase(label)
}
sim_groups <- function(d,cols) {
  if(!nrow(d))return(list())
  split(seq_len(nrow(d)),do.call(paste,c(d[cols],sep="\r")))
}
sim_mean <- function(x) if(any(is.finite(x)))mean(x[is.finite(x)]) else NA_real_
sim_sd <- function(x) if(sum(is.finite(x))>1L)sd(x[is.finite(x)]) else NA_real_
sim_mcse <- function(x) sim_sd(x)/sqrt(sum(is.finite(x)))
sim_label <- function(x) {
  labels<-c(PM25="PM₂.₅",PMcoarse="PM₂.₅–₁₀",O3="O₃",NO2="NO₂",SO2="SO₂",CO="CO",
    temperature="Temperature",humidity="Relative humidity",KDR="KDR-CC",SplineCLR="Spline-CLR",LinearCLR="Linear-CLR",
    single="Single-exposure effect",joint="Joint exposure effect",individual="Observed-combination effect",
    G1="G1: g and f correct",G2="G2: g correct, f incorrect",G3="G3: g incorrect, f correct",G4="G4: g and f incorrect")
  out<-unname(labels[x]);out[is.na(out)]<-x[is.na(out)];out
}
