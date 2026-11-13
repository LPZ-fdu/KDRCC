# Saved-model integration checks run on temporary copies of selected records.
run_saved_contribution_checks <- function(source_root,results_root,output_root=tempdir()) {
  rows<-list();check<-function(name,ok){rows[[length(rows)+1L]]<<-data.frame(check=name,passed=isTRUE(ok));if(!isTRUE(ok))stop(name)}
  for(module in c("C","E")) {
    candidates<-list.dirs(results_root,recursive=FALSE,full.names=TRUE)
    candidates<-candidates[startsWith(basename(candidates),paste0("module_",module,"_"))]
    run<-tail(candidates,1);pp<-readRDS(file.path(run,"plan.rds"));plan<-pp$tasks
    selected<-which(plan$null_count%in%c(8L,80L)&plan$replicate==1L)
    if(module=="E")selected<-selected[plan$strategy[selected]=="Adaptive"]
    selected<-selected[!duplicated(plan$null_count[selected])];pp$tasks<-plan[selected,,drop=FALSE]
    dest<-file.path(output_root,paste0("saved_",module));dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    copied<-character();original<-list()
    for(i in seq_len(nrow(pp$tasks))) {
      key<-pp$tasks$key[i];a<-readRDS(file.path(run,"fits",paste0(key,".rds")))
      original[[key]]<-a;fs<-sim_contribution_files(run,a)
      for(f in fs){relative<-substring(f,nchar(run)+2L);to<-file.path(dest,relative)
        dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(f,to,overwrite=TRUE))}
      copied<-c(copied,unname(fs))
    }
    kdr_save(pp,file.path(dest,"plan.rds"));file.copy(file.path(run,"manifest.rds"),file.path(dest,"manifest.rds"),overwrite=TRUE)
    dir.create(file.path(dest,"source_snapshot","R"),recursive=TRUE,showWarnings=FALSE)
    file.copy(file.path(run,"source_snapshot","R","45_interpretation.R"),
      file.path(dest,"source_snapshot","R","45_interpretation.R"),overwrite=TRUE)
    before_original<-tools::md5sum(unique(copied))
    prepare_contributions(dest,workers=2L,source_root=source_root)
    values<-lapply(pp$tasks$key,function(k)readRDS(file.path(dest,"fits",paste0(k,".rds"))));names(values)<-pp$tasks$key
    strip<-function(a){a$interpretation$contributions<-NULL
      for(n in grep("^contribution_",names(a$interpretation),value=TRUE))a$interpretation[[n]]<-NULL
      a}
    check(paste(module,"preserves all fitted and other inference fields"),all(vapply(names(values),function(k)
      identical(strip(original[[k]]),strip(values[[k]])),logical(1))))
    check(paste(module,"stores eight active contributions"),all(vapply(values,function(a)
      nrow(a$interpretation$contributions)==8L&&all(a$interpretation$contributions$true_active),logical(1))))
    check(paste(module,"normalized active shares"),all(vapply(values,function(a)
      abs(sum(a$interpretation$contributions$estimated_share)-100)<1e-10,logical(1))))
    check(paste(module,"saved model allocation closure"),all(vapply(values,function(a)
      a$interpretation$contribution_identity_error<1e-8,logical(1))))
    check(paste(module,"archived contribution source matches evaluation"),
      unname(tools::md5sum(file.path(dest,"source_snapshot","R","45_interpretation.R")))==
      unname(tools::md5sum(file.path(source_root,"R","45_interpretation.R"))))
    check(paste(module,"does not retain alternate contribution result files"),!dir.exists(file.path(dest,"contributions"))&&
      !length(list.files(dest,pattern="[.]previous$",recursive=TRUE)))
    fits<-file.path(dest,"fits",paste0(pp$tasks$key,".rds"));before<-tools::md5sum(fits)
    resumed<-prepare_contributions(dest,workers=2L,source_root=source_root)
    check(paste(module,"resume preserves completed records"),resumed$cached==length(fits)&&identical(before,tools::md5sum(fits)))
    collected<-sim_collect(dest)
    check(paste(module,"summary uses eight active exposures"),setequal(unique(collected$interpretation$contributions$variable),exposure_names)&&
      !"null_contribution_fraction"%in%names(collected$extended$scores))
    check(paste(module,"read-only source results"),identical(before_original,tools::md5sum(unique(copied))))
  }
  # The eight-dimensional reporting target agrees with the saved formal allocation.
  b<-tail(list.dirs(results_root,recursive=FALSE,full.names=TRUE)[
    startsWith(basename(list.dirs(results_root,recursive=FALSE,full.names=TRUE)),"module_B_")],1)
  plan<-readRDS(file.path(b,"plan.rds"))$tasks;j<-plan[plan$method=="KDR"&plan$scenario_gf=="G1",][1,]
  a<-readRDS(file.path(b,"fits",paste0(j$key,".rds")));fs<-sim_contribution_files(b,a)
  fresh<-sim_active_contributions(a$fit,readRDS(fs["kernel"]),readRDS(fs["evaluation"]),readRDS(fs["config"]))
  check("B eight-exposure shares agree with saved values",max(abs(fresh$contributions$estimated_share-
    a$interpretation$contributions$estimated_share))<1e-10)
  out<-kdr_bind(rows);dir.create(file.path(output_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  kdr_csv(out,file.path(output_root,"validation","saved_contribution_checks.csv"))
  message("Saved contribution checks passed: ",nrow(out));invisible(out)
}
