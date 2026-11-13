# Overall-correlation summaries and complete fixed-scenario reports.
sim_D_panel_layout <- function(n,methods) {
  cols<-min(2L,n);rows<-ceiling(n/cols)
  ids<-seq_len(n);if(length(ids)<rows*cols)ids<-c(ids,rep(n,rows*cols-length(ids)))
  layout(rbind(rep(1L,cols),matrix(ids+1L,ncol=cols,byrow=TRUE)),heights=c(.45,rep(4.2,rows)))
  old<-par("mar");par(mar=c(0,0,0,0));plot.new()
  colors<-c(KDR="#0072B2",SplineCLR="#D55E00",LinearCLR="#009E73")
  legend("center",sim_label(methods),col=colors[methods],lty=1,pch=16,lwd=2,bty="n",horiz=TRUE,cex=.9)
  par(mar=old)
}

sim_D_attach <- function(x,design,correlations) {
  if(!is.data.frame(x)||!nrow(x)||!"config_id"%in%names(x))return(x)
  cols<-c("correlation_id","correlation_label","correlation_mix")
  m<-match(x$config_id,design$config_id)
  for(n in cols)x[[n]]<-design[[n]][m]
  x$mean_absolute_correlation<-correlations$mean_absolute_correlation[match(x$correlation_id,correlations$correlation_id)]
  x
}

sim_D_subset_result <- function(result,id) {
  out<-result
  subset<-function(d)if(is.data.frame(d)&&"config_id"%in%names(d))d[d$config_id==id,,drop=FALSE] else d
  for(n in c("status","selection","effects","points","repetitions","cost"))out[[n]]<-subset(out[[n]])
  out$extended<-lapply(out$extended,subset);out$interpretation<-lapply(out$interpretation,subset)
  out$plan$tasks<-subset(out$plan$tasks);out$plan$design$configurations<-subset(out$plan$design$configurations)
  out$or_errors<-out$or_errors[startsWith(names(out$or_errors),paste0(id,"|"))]
  out
}

sim_D_trend_figure <- function(data,value,label,id,title,cfg,figures,captions,by_target=TRUE,reference=NULL) {
  if(!nrow(data)||!value%in%names(data))return(invisible(NULL))
  panels<-list()
  for(g in paste0("G",1:4))for(t in if(by_target)c("joint","individual") else "all") {
    z<-data[data$scenario_gf==g,,drop=FALSE]
    if(by_target)z<-z[z$target==t,,drop=FALSE]
    z<-z[is.finite(z[[value]])&is.finite(z$mean_absolute_correlation),,drop=FALSE]
    if(nrow(z))panels[[length(panels)+1L]]<-list(data=z,g=g,target=t)
  }
  if(!length(panels))return(invisible(NULL))
  colors<-c(KDR="#0072B2",SplineCLR="#D55E00",LinearCLR="#009E73")
  sim_figure(file.path(figures,paste0(id,".png")),title,ceiling(length(panels)/2),min(2,length(panels)),cfg,function(){
    sim_D_panel_layout(length(panels),names(colors)[names(colors)%in%data$method])
    for(a in panels) {
      z<-a$data;main<-paste(a$g,if(a$target=="all")"" else sim_label(a$target))
      plot(sim_range(z$mean_absolute_correlation),sim_range(z[[value]],reference),type="n",
        xlab="Mean absolute exposure correlation",ylab=label,main=main)
      if(!is.null(reference))abline(h=reference,lty=3,col="grey60")
      methods<-names(colors)[names(colors)%in%z$method]
      for(m in methods) {
        b<-z[z$method==m,,drop=FALSE];b<-b[order(b$correlation_mix),,drop=FALSE]
        lines(b$mean_absolute_correlation,b[[value]],type="b",pch=16,col=colors[m],lwd=2)
      }
    }
  })
  captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
    note=paste("Each point summarizes usable independent repetitions under one overall-correlation setting.",
      "The horizontal coordinate is the mean absolute off-diagonal Pearson correlation in the independent evaluation distribution.",
      "Joint contrasts use a common concentration path; observed-combination contrasts use the distribution of the corresponding setting.",
      "Model configuration and response coefficients are fixed. The matching correct/incorrect CLR adjustment fits supply G1/G2 and G3/G4 displays.",
      "Original results are reused from B; moderate/high results are new D repetitions. B and D batches are independent.",
      "Tables report usable counts and Monte Carlo uncertainty."))
  invisible(NULL)
}

sim_report_D <- function(run,result,tables,figures,captions) {
  result<-sim_D_combine_B(run,result)
  cfg<-result$config;design<-result$plan$design$configurations
  levels<-result$plan$design$correlation_levels;analysis<-dirname(tables)
  ev<-lapply(levels$correlation_id,function(id)readRDS(sim_evaluation_path(run,8L,id)))
  correlations<-kdr_bind(lapply(seq_len(nrow(levels)),function(i) {
    R<-ev[[i]]$overall_correlation;v<-R[upper.tri(R)]
    cbind(levels[i,,drop=FALSE],data.frame(source=if(levels$correlation_id[i]=="D0")"Reused module B" else "New module D",
      mean_absolute_correlation=mean(abs(v)),
      minimum_signed_correlation=min(v),maximum_signed_correlation=max(v),
      minimum_absolute_correlation=min(abs(v)),maximum_absolute_correlation=max(abs(v)),
      evaluation_matched_sets=ev[[i]]$correlation_sample_sets,evaluation_rows=ev[[i]]$correlation_sample_rows))
  }))
  emit<-function(d,id,title,note)sim_export_table(d,id,title,note,tables,captions)
  add<-function(x)sim_D_attach(x,design,correlations)
  import<-readRDS(file.path(run,"baseline_B","import_complete.rds"))
  emit(data.frame(module_B_run=import$run_id,source_signature=import$signature,
    reused_fit_records=import$fit_records,reused_repetitions=import$repetitions),
    "Table_D_baseline_source","Source of the original-correlation reference",
    "Saved B fits are reused after scientific-setting and checksum verification. The D run contains their snapshot and file hashes in baseline_B. B and D use independent simulation batches.")
  emit(cfg$D$copula$levels,"Table_D_correlation_calibration","Independent calibration of exposure dependence",
    "Marginal maps, calibration and validation samples are independent. The calibrated score-mixing parameter is fixed before fitting. Target mean absolute Pearson correlations and pairwise lower bounds refer to final lag-defined exposures.")
  emit(cfg$D$copula$marginals,"Table_D_marginal_distributions","Marginal exposure distributions after recoupling",
    "Means and standard deviations use the independent calibration-validation sample. Fixed empirical marginal transforms approximate the original univariate distributions. Recoupling changes the joint exposure dependence, including its temporal structure and associations with the time-varying confounder.")
  emit(correlations,"Table_D_correlations","Overall exposure correlation settings",
    "Original results reuse the completed module B fits. Moderate and high dependence are calibrated on final lag-defined exposures using fixed marginal maps and a common Gaussian score. Target means and minimum absolute correlations are checked on independent calibration and evaluation samples. All 28 exposure pairs enter the mean absolute coefficient. The recorded mixing parameter is internal to calibration.")
  matrices<-kdr_bind(lapply(seq_along(ev),function(i){
    ix<-which(upper.tri(ev[[i]]$overall_correlation),arr.ind=TRUE)
    data.frame(correlation_id=levels$correlation_id[i],exposure_1=exposure_names[ix[,1]],
      exposure_2=exposure_names[ix[,2]],Pearson_correlation=ev[[i]]$overall_correlation[ix])
  }))
  emit(matrices,"Table_D_pairwise_correlations","Pairwise overall exposure correlations",
    "Pearson correlations use all candidate-date exposure measurements in the independent evaluation sample. Candidate-date rows retain the prescribed fixed-lag exposure definitions.")
  id<-"Figure_D_correlations";title<-"Overall exposure correlation structures"
  cor_cfg<-cfg;cor_cfg$report$width_inches<-max(cfg$report$width_inches,16)
  sim_figure(file.path(figures,paste0(id,".png")),title,ceiling(nrow(levels)/3),min(3,nrow(levels)),cor_cfg,function(){
    par(oma=c(.2,.2,2,.2))
    pal<-colorRampPalette(c("#2166AC","white","#B2182B"))(101)
    for(i in seq_along(ev)) {
      R<-ev[[i]]$overall_correlation
      par(mar=c(6.5,6.2,2.5,.4),pty="s")
      plot(c(.5,8.5),c(.5,8.5),type="n",axes=FALSE,xlab="",ylab="",main=levels$correlation_label[i],asp=1)
      for(j in 1:8)for(k in 1:8) {
        value<-R[j,k];index<-max(1,min(101,round(50*(value+1))+1))
        rect(j-.5,8-k+.5,j+.5,8-k+1.5,col=pal[index],border="white")
        text(j,9-k,formatC(value,digits=2,format="f"),cex=.85,col=if(abs(value)>.65)"white" else "black")
      }
      axis(1,1:8,sim_label(exposure_names),las=2,cex.axis=.88,tick=FALSE)
      axis(2,1:8,rev(sim_label(exposure_names)),las=1,cex.axis=.88,tick=FALSE)
    }
  })
  captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
    note="Cell labels are observed Pearson coefficients in the independent evaluation distribution. The shared color scale spans -1 (blue) to +1 (red). Original data and fits come from B. D recouples the final lag-defined exposures using the original marginal distributions; response coefficients and the fixed model configuration are shared. The temporal dependence of exposure trajectories also changes under this recoupling.")
  effects<-sim_B_display_scenarios(add(result$effects))
  emit(effects,"Table_D_performance","Effect estimation across overall correlations",
    "Joint paths and their true contrasts are common across settings. Actual-combination contrasts use the distribution of each setting with the common median reference. Bias and RMSE use log OR. OR errors are percentages. CLR adjustment fits are reused for their corresponding nuisance scenarios. Repetition counts, uncertainty and Monte Carlo errors accompany each metric.")
  settings<-list(
    c("RMSE","RMSE (log OR)","RMSE","Effect estimation error across exposure correlations"),
    c("RMS_Bias","RMS bias (log OR)","bias","Systematic error across exposure correlations"),
    c("mean_OR_absolute_relative_error","Mean absolute OR relative error (%)","OR_error","Odds-ratio error across exposure correlations"),
    c("SE_ratio","RMS model SE / RMS empirical SD","SE_calibration","Standard-error calibration across exposure correlations"),
    c("coverage_percent","Pointwise interval coverage (%)","coverage","Interval coverage across exposure correlations"),
    c("correct_significant_percent","Correct significant direction (%)","effect_detection","Effect detection across exposure correlations"))
  for(s in settings)sim_D_trend_figure(effects,s[1],s[2],paste0("Figure_D_",s[3]),s[4],cfg,figures,captions,
    reference=if(s[1]=="SE_ratio")1 else if(s[1]=="coverage_percent")95 else NULL)
  scores<-sim_B_display_scenarios(add(result$extended$scores))
  if(nrow(scores)) {
    scores$contribution_error_pp<-100*scores$contribution_share_error
    sim_D_trend_figure(scores,"contribution_error_pp","Contribution share RMSE (pp)",
      "Figure_D_contribution_accuracy","Contribution accuracy across exposure correlations",cfg,figures,captions,FALSE)
  }
  for(n in names(result$extended))emit(sim_B_display_scenarios(add(result$extended[[n]])),paste0("Table_D_correlation_",n),
    paste("Correlation-specific",gsub("_"," ",n)),
    "Each setting uses the same model configuration and its own independent evaluation distribution for observed combinations. Constructed contrasts share the reference concentrations. Counts and Monte Carlo errors are based on complete independent repetitions.")
  pair<-sim_B_display_scenarios(add(result$extended$pairs));direction<-sim_B_display_scenarios(add(result$extended$directions))
  for(g in paste0("G",1:4)) {
    panels<-list()
    for(component in unique(pair$component))for(mode in c("overall","direction")) {
      z<-if(mode=="overall")pair else direction
      z<-z[z$scenario_gf==g&z$component==component,,drop=FALSE]
      field<-if(mode=="overall")"corrected_detection_percent" else "corrected_correct_percent"
      if(nrow(z)&&any(is.finite(z[[field]])))panels[[length(panels)+1L]]<-list(z=z,component=component,field=field,mode=mode)
    }
    if(!length(panels))next
    id<-paste0("Figure_D_interaction_detection_",g);title<-paste("Interaction detection across exposure correlations:",g)
    sim_figure(file.path(figures,paste0(id,".png")),title,ceiling(length(panels)/2),min(2,length(panels)),cfg,function(){
      colors<-c(KDR="#0072B2",SplineCLR="#D55E00",LinearCLR="#009E73")
      sim_D_panel_layout(length(panels),names(colors)[names(colors)%in%pair$method])
      for(a in panels) {
        z<-a$z
        plot(sim_range(z$mean_absolute_correlation),c(0,105),type="n",xlab="Mean absolute exposure correlation",
          ylab=if(a$mode=="overall")"Overall detection (%)" else "Correct-direction detection (%)",
          main=paste(sim_label(strsplit(a$component,"__",fixed=TRUE)[[1]]),collapse=" × "))
        methods<-names(colors)[names(colors)%in%z$method]
        for(m in methods){b<-z[z$method==m,,drop=FALSE];b<-b[order(b$correlation_mix),]
          lines(b$mean_absolute_correlation,b[[a$field]],type="b",pch=16,col=colors[m],lwd=2)}
      }
    })
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note="Left panels summarize Holm-corrected overall pairwise detection. Right panels summarize correct significant pointwise directions with Holm correction across the prespecified interaction grid. The same four pairs and constructed grids are evaluated at every correlation setting. Additive comparators have zero four-point interaction contrasts. Tables report usable tests, Monte Carlo uncertainty and wrong-direction rates.")
  }
  selection<-sim_selection_summary(result$selection,result$status)
  variables<-add(selection$variables)
  emit(add(selection$overall),"Table_D_correlation_selection","Active-exposure identification across correlations",
    "Mean true-positive rates average the selected fraction of eight active exposures across complete usable repetitions. Variable tests use Holm adjustment. This module contains eight active exposures.")
  emit(variables,"Table_D_correlation_active_variables","Exposure-specific identification across correlations",
    "Rates are computed across usable KDR-CC G1 repetitions. Limits are 95% Wilson Monte Carlo intervals. Humidity enters through its interaction with temperature.")
  active<-exposure_names[exposure_names%in%variables$variable]
  if(length(active)) {
    id<-"Figure_D_active_exposures";title<-"Active-exposure identification across correlations"
    sim_figure(file.path(figures,paste0(id,".png")),title,ceiling(length(active)/2),min(2,length(active)),cfg,function(){
      sim_panel_layout(length(active))
      for(v in active) {
        z<-variables[variables$variable==v,,drop=FALSE];z<-z[order(z$correlation_mix),]
        plot(z$mean_absolute_correlation,100*z$rate,type="b",pch=16,lwd=2,col="#0072B2",ylim=c(0,105),
          xlab="Mean absolute exposure correlation",ylab="Identification rate (%)",main=sim_label(v))
        good<-is.finite(z$lower)&is.finite(z$upper)&z$upper>z$lower
        if(any(good))arrows(z$mean_absolute_correlation[good],100*z$lower[good],z$mean_absolute_correlation[good],
          100*z$upper[good],angle=90,code=3,length=.035,col="#0072B2")
      }
    })
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note="Points are identification rates across KDR-CC G1 repetitions; bars are Wilson 95% Monte Carlo intervals. Selection uses Holm correction across the eight candidates. Overall correlations are measured in the fixed independent evaluation distribution. Tables give planned and usable counts.")
  }
  # Complete B-style assessment within each correlation setting, sharing calculations.
  for(i in seq_len(nrow(design))) {
    id<-design$config_id[i];cid<-design$correlation_id[i];part<-sim_D_subset_result(result,id)
    source_run<-if(cid=="D0")file.path(run,"baseline_B") else run
    part$config$report$correlation_id<-if(cid=="D0")NULL else cid
    folder<-file.path(analysis,"correlation_scenarios",cid)
    ft<-file.path(folder,"tables");ff<-file.path(folder,"figures")
    dir.create(ft,recursive=TRUE,showWarnings=FALSE);dir.create(ff,recursive=TRUE,showWarnings=FALSE)
    subcap<-new.env(parent=emptyenv())
    sim_plot_B(part,ff,subcap)
    sim_report_interpretation(source_run,part,ft,ff,subcap,"B")
    sim_report_B_responses(source_run,part,ft,ff,subcap)
    rename<-function(x)sub("(^|/)(Figure_|Table_)?B_",paste0("\\1\\2D_",cid,"_"),x)
    for(f in list.files(folder,recursive=TRUE,full.names=TRUE)) {
      to<-file.path(dirname(f),rename(basename(f)))
      if(to!=f&&!file.rename(f,to))stop("Unable to name the correlation-specific output: ",f)
    }
    for(k in ls(subcap)) {
      a<-subcap[[k]];a$file<-paste0("correlation_scenarios/",cid,"/",rename(a$file))
      a$title<-paste0(a$title," | ",design$correlation_label[i])
      a$note<-paste(a$note,"Overall correlation setting:",design$correlation_label[i],
        "; mean absolute Pearson coefficient =",format(correlations$mean_absolute_correlation[match(cid,correlations$correlation_id)],digits=4),".")
      captions[[rename(k)]]<-a
    }
    kdr_log("Correlation-specific report prepared:",cid)
  }
  invisible(NULL)
}
