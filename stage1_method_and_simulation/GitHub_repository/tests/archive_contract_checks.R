# Read-only comparison of a release with an archived formal run.
run_archive_contract_checks <- function(source_root, results_root, output_root=tempdir()) {
  load_env <- function(root) {
    e<-new.env(parent=.GlobalEnv)
    for(f in sort(list.files(file.path(root,"R"),pattern="[.]R$",full.names=TRUE)))
      sys.source(f,e,keep.source=FALSE)
    sys.source(file.path(root,"config","mechanism.R"),e,keep.source=FALSE)
    e
  }
  current<-load_env(source_root);rows<-list()
  for(run in list.dirs(results_root,recursive=FALSE,full.names=TRUE)) {
    if(!file.exists(file.path(run,"config.rds")))next
    cfg<-readRDS(file.path(run,"config.rds"));module<-readRDS(file.path(run,"manifest.rds"))$module
    snapshot<-file.path(run,"source_snapshot")
    old<-load_env(snapshot)
    # The common, original-dependence generator and A-D nuisance dictionary.
    cfg$baseline$dgp$copula_eta<-NULL;cfg$nuisance$raw_spec<-NULL
    a<-current$kdr_generate(24L,12L,cfg,712390L)
    b<-old$kdr_generate(24L,12L,cfg,712390L)
    checks<-c(generator=identical(a,b),
      ordered_pairs=identical(current$kdr_pairs(a),old$kdr_pairs(b)),
      nuisance_dictionary=identical(current$kdr_nuisance_raw(a$A,a$W,cfg),old$kdr_nuisance_raw(b$A,b$W,cfg)))
    core<-c("01_numeric.R","02_calendar.R","20_kernel.R","30_estimation.R","31_covariance.R",
            "40_evaluation.R","41_response_grid.R","42_attribution.R","43_interaction_inference.R")
    same<-vapply(core,function(f)identical(readBin(file.path(source_root,"R",f),"raw",n=file.info(file.path(source_root,"R",f))$size),
      readBin(file.path(snapshot,"R",f),"raw",n=file.info(file.path(snapshot,"R",f))$size)),logical(1))
    checks<-c(checks,setNames(same,paste0("core_bytes:",core)))
    stopifnot(all(checks))
    rows[[length(rows)+1L]]<-data.frame(module=module,check=names(checks),passed=unname(checks),stringsAsFactors=FALSE)
  }
  out<-do.call(rbind,rows);dir.create(file.path(output_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  write.csv(out,file.path(output_root,"validation","archive_contract_checks.csv"),row.names=FALSE,fileEncoding="UTF-8")
  message("Archive contract checks passed: ",nrow(out));invisible(out)
}
