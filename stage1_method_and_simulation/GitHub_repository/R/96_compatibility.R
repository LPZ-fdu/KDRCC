# Fast checks run on the interpreter that will actually execute the study.
sim_check_compatibility <- function(source_root,work_root,parallel_check=FALSE) {
  checks<-list();add<-function(name,passed,detail="") {
    checks[[length(checks)+1L]]<<-data.frame(check=name,passed=isTRUE(passed),detail=detail,stringsAsFactors=FALSE)
    if(!isTRUE(passed))stop("Compatibility check failed: ",name,". ",detail)
  }
  add("R >= 3.6.3",getRversion()>="3.6.3",as.character(getRversion()))
  add("64-bit R",.Machine$sizeof.pointer==8L,R.version$arch)
  for(p in c("parallel","splines","CompQuadForm"))
    add(paste("Package",p),requireNamespace(p,quietly=TRUE),"Run the install action if missing.")
  ns<-asNamespace("parallel")
  add("PSOCK asynchronous interfaces",all(vapply(c("sendCall","recvOneResult"),exists,logical(1),envir=ns,inherits=FALSE)))
  for(f in list.files(file.path(source_root,"R"),pattern="[.]R$",full.names=TRUE))parse(f,encoding="UTF-8")
  add("Core syntax on current R",TRUE)
  probe<-data.frame(label="G1")
  add("Character-column convention",is.character(probe$label))
  temp<-tempfile("rds_compatibility_",tmpdir=work_root);on.exit(unlink(temp),add=TRUE)
  value<-list(x=c(1,NA_real_),name="KDR-CC")
  saveRDS(value,temp,version=3)
  add("Serialization v3 read/write",identical(value,readRDS(temp)))
  q<-CompQuadForm::davies(2,lambda=c(1,1),lim=100000L,acc=1e-6)
  add("Numerical quadratic-form tail",q$ifault==0L&&is.finite(q$Qq)&&abs(q$Qq-pchisq(2,2,lower.tail=FALSE))<2e-6)
  add("R 3.6 color palette API",length(grDevices::hcl.colors(4,"YlOrRd"))==4L)
  if(parallel_check) {
    vars<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    old<-Sys.getenv(vars,unset=NA_character_)
    on.exit({for(v in vars)if(is.na(old[v]))Sys.unsetenv(v) else do.call(Sys.setenv,setNames(list(old[v]),v))},add=TRUE)
    Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
    cl<-parallel::makePSOCKcluster(2L);on.exit(parallel::stopCluster(cl),add=TRUE)
    got<-parallel::clusterCall(cl,function(root,libs){
      .libPaths(libs);source(file.path(root,"load.R"),encoding="UTF-8");load_kdrcc_simulation(root)
      list(R=as.character(getRversion()),package=as.character(packageVersion("CompQuadForm")),
        pid=Sys.getpid(),character_columns=is.character(data.frame(x="G1")$x),
        seed=stage_seed(27182818L,"compatibility",1L))
    },source_root,.libPaths())
    add("Two independent PSOCK workers",length(unique(vapply(got,`[[`,integer(1),"pid")))==2L)
    add("Worker R and dependency versions",all(vapply(got,function(x)x$R==as.character(getRversion())&&
      x$package==as.character(packageVersion("CompQuadForm"))&&x$character_columns,logical(1))))
    add("Worker seed parity",all(vapply(got,`[[`,integer(1),"seed")==stage_seed(27182818L,"compatibility",1L)))
  }
  if(!isTRUE(capabilities("cairo")))kdr_log("Cairo PNG is unavailable on this R build. Fitting is supported; prepare figures on a Cairo-enabled computer.")
  if(.Platform$OS.type!="windows")kdr_log("Linux: Windows memory sampling is disabled; no zero-memory values will be substituted.")
  d<-kdr_bind(checks)
  dir.create(file.path(work_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  kind<-if(parallel_check)"full" else "runtime"
  kdr_csv(d,file.path(work_root,"validation",paste0("compatibility_checks_",kind,".csv")))
  capture.output(sessionInfo(),file=file.path(work_root,"validation",paste0("compatibility_sessionInfo_",kind,".txt")))
  kdr_log("Compatibility checks passed:",nrow(d),"| R",as.character(getRversion()),"| CompQuadForm",as.character(packageVersion("CompQuadForm")))
  invisible(d)
}
