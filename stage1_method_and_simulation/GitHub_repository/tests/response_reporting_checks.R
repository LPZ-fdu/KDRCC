# Explicit saved-fit reporting check; no simulation or model fitting is run.
run_response_reporting_checks <- function(source_root,run,render=TRUE) {
  source(file.path(source_root,"load.R"),encoding="UTF-8");load_kdrcc_simulation(source_root)
  cfg<-readRDS(file.path(run,"config.rds"));pp<-readRDS(file.path(run,"plan.rds"))
  stopifnot(readRDS(file.path(run,"manifest.rds"))$module=="B")
  raw_files<-unlist(lapply(c("fits","kernels","comparators","evaluation"),function(d)
    list.files(file.path(run,d),full.names=TRUE)))
  raw_files<-c(raw_files,file.path(run,c("config.rds","plan.rds","manifest.rds")))
  before<-tools::md5sum(raw_files)
  settings<-sim_B_response_settings(list(),cfg)
  ev<-readRDS(file.path(run,"evaluation","p8.rds"))
  grid<-sim_B_response_grid(ev,cfg,settings)
  tab<-grid$points;stopifnot(length(unique(tab$component[tab$target=="single"]))==8L,
    length(unique(tab$component[tab$target=="interaction"]))==4L,
    all(tab$truth[!tab$eligible]==0))
  # Bilinear generating interactions have an analytic anchored contrast.
  for(pair in unique(tab$component[tab$target=="interaction"])) {
    a<-tab[tab$component==pair,];j<-match(a$exposure_1[1],exposure_names);k<-match(a$exposure_2[1],exposure_names)
    t<-cfg$baseline$truth$interactions
    coef<-sum(t$coefficient[(t$q==j&t$r==k)|(t$q==k&t$r==j)])
    want<-coef*(a$concentration_1-ev$reference[j])*(a$concentration_2-ev$reference[k])/
      (cfg$baseline$reference$scale[j]*cfg$baseline$reference$scale[k])
    stopifnot(max(abs(a$truth-want))<1e-12)
  }
  # A specified null pair remains in the evaluation, with zero truth off the anchors.
  nullset<-sim_B_response_settings(list(B_interaction_pairs=rbind(c("PM25","NO2"))),cfg)
  nullgrid<-sim_B_response_grid(ev,cfg,nullset)
  nullpoints<-nullgrid$points[nullgrid$points$target=="interaction",]
  stopifnot(any(nullpoints$eligible),max(abs(nullpoints$truth))<1e-12)
  for(m in c("KDR","SplineCLR","LinearCLR")) {
    job<-pp$tasks[pp$tasks$method==m,,drop=FALSE][1,]
    path<-if(m=="KDR")file.path(run,"kernels",paste0(job$kernel_key,".rds")) else
      file.path(run,"comparators",paste0(job$kernel_key,"_",m,".rds"))
    spec<-readRDS(path);a<-readRDS(file.path(run,"fits",paste0(job$key,".rds")))
    D<-sim_B_response_operator(grid,spec)
    pred<-drop(D%*%a$fit$beta)
    stopifnot(length(pred)==nrow(tab),all(is.finite(pred)),all(pred[!tab$eligible]==0))
    if(m!="KDR")stopifnot(max(abs(pred[tab$target=="interaction"]))<1e-10)
    # Direct point predictions reproduce the saved joint-curve values.
    ref<-drop(kdr_basis(matrix(ev$reference,1,8),spec))
    direct<-drop(sweep(kdr_basis(ev$joint,spec),2,ref)%*%a$fit$beta)
    saved<-a$outputs[a$outputs$target=="joint","estimate"]
    stopifnot(max(abs(direct-saved))<1e-10)
  }
  Y<-rbind(tab$truth+.03,tab$truth-.03)
  Y[,!tab$eligible]<-0
  summary<-sim_B_response_summary(Y,grid,c(1L,2L),2L)
  stopifnot(max(summary$metrics$RMS_Bias)<1e-12,max(abs(summary$metrics$RMSE-.03))<1e-12)
  Y[1,tab$component=="PM25"]<-NA_real_
  partial<-sim_B_response_summary(Y,grid,c(1L,2L),2L)$metrics
  stopifnot(partial$usable_repetitions[partial$component=="PM25"]==1L,
    all(partial$usable_repetitions[partial$component!="PM25"]==2L))
  out<-sim_B_response_collect(run,cfg,pp,settings)
  stopifnot(all(out$status$prediction_ok))
  display<-sim_B_display_scenarios(out$points)
  for(m in c("SplineCLR","LinearCLR"))for(pair in list(c("G1","G2"),c("G3","G4"))) {
    a<-display[display$method==m&display$scenario_gf==pair[1],]
    b<-display[display$method==m&display$scenario_gf==pair[2],]
    stopifnot(identical(a$mean_estimate,b$mean_estimate),identical(a$usable_repetitions,b$usable_repetitions))
  }
  empty<-sim_B_display_scenarios(data.frame());stopifnot(nrow(empty)==0L)
  if(render) {
    prepare_results(run,list(dpi=100L,width_inches=12,B_response_summaries=TRUE))
    fig<-file.path(run,"analysis","figures")
    stopifnot(all(file.exists(file.path(fig,paste0(c(paste0("Figure_B_single_G",1:4),"Figure_B_interactions"),".png")))))
    met<-read.csv(file.path(run,"analysis","tables","Table_B_single_accuracy.csv"),stringsAsFactors=FALSE)
    stopifnot(nrow(met)==8L*3L*4L)
    met<-read.csv(file.path(run,"analysis","tables","Table_B_interaction_accuracy.csv"),stringsAsFactors=FALSE)
    stopifnot(nrow(met)==4L*3L*4L)
    captions<-read.csv(file.path(run,"analysis","Figure_and_table_captions.csv"),stringsAsFactors=FALSE)
    stopifnot(all(file.exists(file.path(run,"analysis",captions$file))))
  }
  stopifnot(identical(before,tools::md5sum(raw_files)))
  kdr_log("Response reporting checks passed: analytic contrasts, saved-fit predictions, CLR reuse, aggregation, availability, figures and unchanged fit files.")
  invisible(TRUE)
}
