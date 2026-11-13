# Import an immutable module B reference without fitting a baseline model again.
sim_D_source_B <- function(cfg,root,resolved) {
  path<-cfg$D$baseline_B_run
  if(is.null(path))path<-sim_pointer(root,"B")
  if(is.null(path)) {
    candidates<-list.dirs(file.path(root,"results"),recursive=FALSE,full.names=TRUE)
    candidates<-candidates[grepl("^module_B_",basename(candidates))]
    if(length(candidates)==1L)path<-candidates
  }
  if(is.null(path))stop("Set D_baseline_run to the completed module B directory.")
  if(!dir.exists(path)&&dir.exists(file.path(root,"results",path)))path<-file.path(root,"results",path)
  path<-normalizePath(path,winslash="/",mustWork=TRUE)
  required<-file.path(path,c("config.rds","plan.rds","manifest.rds","evaluation/p8.rds"))
  if(!all(file.exists(required)))stop("Module B reference is incomplete: ",path)
  b<-readRDS(required[1]);pp<-readRDS(required[2]);m<-readRDS(required[3])
  if(!identical(m$module,"B")||any(pp$tasks$p!=8L)||nrow(pp$design$configurations)!=1L)
    stop("D requires a completed eight-exposure, single-configuration module B run.")
  pars<-as.list(pp$design$configurations[1,sim_parameter_names,drop=FALSE])
  if(!isTRUE(all.equal(pars,resolved$parameters,check.attributes=FALSE)))
    stop("D parameters must match the reused B configuration; review D_overrides and fixed_parameters.R.")
  for(n in c("baseline","kernel","nuisance","solver","comparator","selection","interpretation","evaluation"))
    if(!isTRUE(all.equal(cfg[[n]],b[[n]],check.attributes=FALSE)))
      stop("D and reused B have different scientific settings: ",n)
  fits<-file.path(path,"fits",paste0(pp$tasks$key,".rds"))
  if(!all(file.exists(fits)))stop("Module B still has pending fit records; complete B before reusing it.")
  key_files<-c(required,fits)
  hashes<-tools::md5sum(key_files);names(hashes)<-substring(key_files,nchar(path)+2L)
  list(path=path,run_id=basename(path),signature=object_id(hashes),hashes=hashes,
    fit_records=length(fits),repetitions=length(unique(pp$tasks$replicate)))
}
sim_D_import_B <- function(cfg,run) {
  info<-cfg$D$source_B;dest<-file.path(run,"baseline_B")
  marker<-sim_read(file.path(dest,"import_complete.rds"))
  if(!is.null(marker)) {
    if(!identical(marker$signature,info$signature))stop("Imported B reference identity differs.")
    return(invisible(dest))
  }
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  kdr_log("Importing completed B reference:",info$run_id,"| saved fits",info$fit_records)
  for(d in c("fits","kernels","kernel_records","comparators","evaluation","source_snapshot")) {
    src<-file.path(info$path,d)
    if(dir.exists(src)) {
      fs<-list.files(src,recursive=TRUE,full.names=TRUE,all.files=TRUE)
      fs<-fs[!file.info(fs)$isdir]
      for(f in fs) {
        out<-file.path(dest,d,substring(f,nchar(src)+2L));dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
        if(!file.copy(f,out,overwrite=TRUE))stop("Failed to import B file: ",f)
      }
    }
  }
  fs<-c("config.rds","plan.rds","manifest.rds","effective_configuration.R","sessionInfo.txt")
  fs<-fs[file.exists(file.path(info$path,fs))]
  if(!all(file.copy(file.path(info$path,fs),dest,overwrite=TRUE)))stop("Failed to import B metadata.")
  copied<-tools::md5sum(file.path(dest,names(info$hashes)))
  if(!identical(unname(copied),unname(info$hashes)))stop("B import checksum verification failed.")
  bcfg<-readRDS(file.path(dest,"config.rds"));ev<-readRDS(file.path(dest,"evaluation","p8.rds"))
  d<-kdr_generate(bcfg$evaluation$groups,8L,bcfg,kdr_seed(bcfg,"evaluation"),events=FALSE)
  rows<-match(seq_len(bcfg$evaluation$function_points),d$group)
  if(!isTRUE(all.equal(ev$individual,d$A[rows,,drop=FALSE],tolerance=1e-12,check.attributes=FALSE)))
    stop("Regenerated B evaluation exposures differ from the saved B evaluation sample.")
  ev$overall_correlation<-cor(d$A);ev$correlation_sample_rows<-nrow(d$A);ev$correlation_sample_sets<-d$n
  ev$correlation_setting<-data.frame(correlation_id="D0",correlation_label="Original (module B)",
    target_mean=NA_real_,min_abs=NA_real_,mean_tolerance=NA_real_,correlation_mix=0)
  kdr_save(ev,sim_evaluation_path(run,8L,"D0"))
  kdr_save(info,file.path(dest,"import_complete.rds"))
  kdr_csv(data.frame(file=names(info$hashes),md5=unname(info$hashes)),file.path(dest,"source_md5.csv"))
  invisible(dest)
}
sim_D_bind <- function(a,b) {
  if(!is.data.frame(a)||!nrow(a))return(b)
  if(!is.data.frame(b)||!nrow(b))return(a)
  cols<-union(names(a),names(b))
  for(n in setdiff(cols,names(a)))a[[n]]<-NA
  for(n in setdiff(cols,names(b)))b[[n]]<-NA
  rbind(a[,cols,drop=FALSE],b[,cols,drop=FALSE])
}
sim_D_combine_B <- function(run,result) {
  b<-sim_collect(file.path(run,"baseline_B"))
  level<-readRDS(sim_evaluation_path(run,8L,"D0"))$correlation_setting
  b$plan$design$configurations<-cbind(b$plan$design$configurations,level)
  for(n in names(level))b$plan$tasks[[n]]<-level[[n]]
  b$status$reused<-TRUE
  if(nrow(b$cost))b$cost$reused<-TRUE
  for(n in c("status","selection","effects","points","repetitions","cost"))
    result[[n]]<-sim_D_bind(b[[n]],result[[n]])
  for(n in names(result$extended))result$extended[[n]]<-sim_D_bind(b$extended[[n]],result$extended[[n]])
  for(n in names(result$interpretation))result$interpretation[[n]]<-sim_D_bind(b$interpretation[[n]],result$interpretation[[n]])
  result$or_errors<-c(b$or_errors,result$or_errors)
  result$plan$tasks<-sim_D_bind(b$plan$tasks,result$plan$tasks)
  result$plan$design$configurations<-sim_D_bind(b$plan$design$configurations,result$plan$design$configurations)
  result$plan$design$correlation_levels<-sim_D_bind(level,result$plan$design$correlation_levels)
  result
}
