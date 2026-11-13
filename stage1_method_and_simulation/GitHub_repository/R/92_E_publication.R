# Publication reports for a selection procedure, averaged over independent repetitions.
sim_E_attach <- function(d,scenarios) {
  if(is.null(d)||!nrow(d))return(data.frame())
  keys<-kdr_bind(lapply(c("Fixed","Adaptive"),function(s){a<-scenarios;a$strategy<-s
    a$config_id<-paste(a$scenario_id,s,sep="_");a}))
  ix<-match(d$config_id,keys$config_id)
  for(n in c("scenario_id","scenario_label","strategy","dependence"))d[[n]]<-keys[[n]][ix]
  d
}
sim_E_paired <- function(d,metric,groups=c("scenario_id","target")) {
  if(!nrow(d))return(data.frame())
  a<-d[d$strategy=="Fixed",c(groups,"replicate",metric),drop=FALSE]
  b<-d[d$strategy=="Adaptive",c(groups,"replicate",metric),drop=FALSE]
  names(a)[ncol(a)]<-"fixed";names(b)[ncol(b)]<-"adaptive"
  z<-merge(a,b,by=c(groups,"replicate"));z<-z[is.finite(z$fixed)&is.finite(z$adaptive),,drop=FALSE]
  if(!nrow(z))return(data.frame())
  z$difference<-z$adaptive-z$fixed;z$metric<-metric;z
}
sim_E_comparison_panel <- function(d,metric,label,scenarios,reference=NULL,error=NULL,axis_label=NULL) {
  if(!nrow(d)||!metric%in%names(d)||!any(is.finite(d[[metric]])))return(invisible(NULL))
  levels<-scenarios$scenario_id;cols<-c(Fixed="#0072B2",Adaptive="#D55E00")
  limits<-d[[metric]]
  if(!is.null(error)&&error%in%names(d))limits<-c(limits,d[[metric]]-1.96*d[[error]],d[[metric]]+1.96*d[[error]])
  yrange<-if(metric%in%c("rate","mean_TPR"))c(0,1) else if(metric=="mean_FPR")
    c(0,min(1,max(.01,limits[is.finite(limits)])*1.1)) else sim_range(limits,reference)
  plot(c(.6,length(levels)+.4),yrange,type="n",xaxt="n",
    xlab="Exposure setting",ylab=if(is.null(axis_label))tail(strsplit(label,"\n",fixed=TRUE)[[1]],1) else axis_label,main=label)
  axis(1,seq_along(levels),labels=FALSE);text(seq_along(levels),par("usr")[3],
    labels=scenarios$scenario_label,srt=30,adj=c(1,1.2),xpd=TRUE,cex=.78)
  if(!is.null(reference))abline(h=reference,lty=3,col="grey60")
  for(s in names(cols)) {
    z<-d[d$strategy==s & is.finite(d[[metric]]),,drop=FALSE];x<-match(z$scenario_id,levels)
    ord<-order(x);x<-x[ord];z<-z[ord,,drop=FALSE];off<-if(s=="Fixed")-.07 else .07
    if(!nrow(z))next
    lines(x+off,z[[metric]],type="b",pch=if(s=="Fixed")16 else 17,col=cols[s],lwd=2)
    if(!is.null(error)&&error%in%names(z)) {
      good<-is.finite(z[[error]])&z[[error]]>0
      if(any(good))arrows(x[good]+off,z[[metric]][good]-1.96*z[[error]][good],
        x[good]+off,z[[metric]][good]+1.96*z[[error]][good],angle=90,code=3,length=.035,col=cols[s])
    }
  }
  legend("topright",c("Fixed configuration","Data-selected configuration"),col=cols,pch=c(16,17),lty=1,bty="n",cex=.76)
}
sim_report_E <- function(run,report_overrides=list()) {
  lock<-file.path(run,"report_lock");if(!dir.create(lock,showWarnings=FALSE))stop("Report lock exists.")
  on.exit(sim_unlock(lock,run),add=TRUE)
  analysis<-file.path(run,"analysis");tables<-file.path(analysis,"tables");figures<-file.path(analysis,"figures")
  dir.create(tables,recursive=TRUE,showWarnings=FALSE);dir.create(figures,showWarnings=FALSE)
  x<-sim_collect(run);cfg<-x$config;cfg$report<-sim_report_font(sim_merge(cfg$report,report_overrides));x$config<-cfg
  if(!isTRUE(capabilities("cairo")))stop("Reporting requires Cairo PNG support.")
  sc<-x$plan$design$scenarios;cap<-new.env(parent=emptyenv())
  emit<-function(d,id,title,note)sim_export_table(d,id,title,note,tables,cap)
  figure<-function(id,title,note,nr,nc,draw) {
    sim_figure(file.path(figures,paste0(id,".png")),title,nr,nc,cfg,draw)
    cap[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,note=note)
  }
  attach<-function(d)sim_E_attach(d,sc)
  status<-attach(x$status);kdr_csv(status,file.path(tables,"task_status.csv"))
  counts<-kdr_bind(lapply(sim_groups(status,c("scenario_id","strategy")),function(ix){z<-status[ix,]
    data.frame(scenario=z$scenario_label[1],strategy=z$strategy[1],planned=nrow(z),
      completed=sum(z$completed),point_estimates=sum(z$point_ok),intervals=sum(z$interval_ok),
      variable_tests=sum(z$selection_ok),interaction_tests=sum(z$interaction_ok))}))
  emit(counts,"Table_E_availability","Availability of final paired estimates",
    "Each scenario uses new independent repetitions. Both strategies use identical final-estimation data. Completed records include recorded failures; usable estimates, intervals and tests have separate denominators.")
  emit(x$plan$design$configurations,"Table_E_candidates","Prespecified configuration candidates",
    "N and L are final-estimation and independent kernel-learning matched-set counts. M, kappa and lambda are selected by the minimum pooled held-out conditional negative log likelihood. The frozen configuration is always included. Three validation folds each use two-thirds of the configuration pool for estimation. The independent kernel-learning pool is common to all folds.")
  effects<-attach(x$effects);points<-attach(x$points);rr<-attach(x$repetitions)
  emit(effects,"Table_E_effect_performance","Effect recovery after data-driven parameter selection",
    "RMS bias is the root mean squared pointwise bias across independent repetitions; RMSE averages squared errors over repetitions and evaluation points. Errors use log OR; absolute relative OR error is 100*abs(exp(estimate-truth)-1). Empirical and model SE use RMS aggregation. Intervals condition on the independently learned kernel and selected configuration. The empirical evaluations include configuration selection variability. Selection never uses these truth-based metrics.")
  kdr_csv(rr,file.path(tables,"effect_metrics_by_repetition.csv"));kdr_csv(points,file.path(tables,"pointwise_effect_summary.csv"))
  for(block in c("accuracy","uncertainty")) {
    measures<-if(block=="accuracy")c("RMS_Bias","RMSE","mean_OR_absolute_relative_error") else c("coverage_percent","SE_ratio","mean_interval_width")
    eligible<-unlist(lapply(measures,function(m)lapply(c("joint","individual"),function(t) {
      z<-effects[effects$target==t,,drop=FALSE]
      if(nrow(z)&&any(is.finite(z[[m]])))list(list(metric=m,target=t,data=z)) else list()
    })),recursive=FALSE)
    if(length(eligible)) {
      # Explicit panel lists keep an unavailable target out of the figure.
      panels<-list();for(m in measures)for(t in c("joint","individual")) {
        z<-effects[effects$target==t,,drop=FALSE]
        if(nrow(z)&&any(is.finite(z[[m]])))panels[[length(panels)+1L]]<-list(m=m,t=t,z=z)
      }
      if(length(panels))figure(paste0("Figure_E_",block),if(block=="accuracy")"Effect accuracy across exposure settings" else "Uncertainty calibration across exposure settings",
        "Fixed and data-selected strategies use matched final samples. Error bars, where available, are 95% Monte Carlo intervals across independent repetitions. Coverage refers to nominal pointwise 95% intervals.",
        ceiling(length(panels)/2),min(2,length(panels)),function(){for(a in panels) {
          error<-switch(a$m,RMS_Bias="MCSE_RMS_Bias",RMSE="MCSE_RMSE",mean_OR_absolute_relative_error="MCSE_mean_OR_error",coverage_percent="MCSE_coverage_percent",NULL)
          sim_E_comparison_panel(a$z,a$m,paste(sim_label(a$t),sim_metric_label(a$m),sep="\n"),sc,
            reference=if(a$m=="coverage_percent")95 else if(a$m=="SE_ratio")1 else NULL,error=error)
        }})
    }
  }
  paired<-kdr_bind(lapply(c("RMSE","mean_OR_error","coverage","correct_significant"),function(m)sim_E_paired(rr,m)))
  kdr_csv(paired,file.path(tables,"paired_effect_differences.csv"))
  summarize<-function(d,groups)kdr_bind(lapply(sim_groups(d,groups),function(ix){z<-d[ix,];n<-nrow(z);se<-sim_mcse(z$difference)
    cbind(z[1,groups,drop=FALSE],data.frame(pairs=n,fixed_mean=mean(z$fixed),adaptive_mean=mean(z$adaptive),
      mean_difference=mean(z$difference),MCSE=se,lower_MC=mean(z$difference)-if(n>1)qt(.975,n-1)*se else NA_real_,
      upper_MC=mean(z$difference)+if(n>1)qt(.975,n-1)*se else NA_real_))}))
  emit(summarize(paired,c("scenario_id","target","metric")),"Table_E_paired_changes","Paired changes in final performance",
    "Differences are adaptive minus fixed within the same repetition, restricted to pairs with both available values. Negative RMSE and OR-error changes favor adaptive selection. Positive coverage and correct-direction detection changes indicate increases. The table reports mean repetition-level changes; aggregate RMS bias and pooled RMSE appear in the effect-performance table.")
  sel<-sim_selection_summary(x$selection,x$status);sel$overall<-attach(sel$overall);sel$variables<-attach(sel$variables)
  emit(sel$overall,"Table_E_variable_selection","Effective-variable identification and null-variable inclusion",
    "TPR and FPR are averaged over complete usable repetitions after Holm correction. FPR is undefined when the null count is zero. Monte Carlo intervals use whole-repetition rates.")
  emit(sel$variables,"Table_E_active_rates","Identification of each active exposure",
    "Rates use complete usable variable-test repetitions. Humidity acts through interaction. Wilson intervals quantify Monte Carlo uncertainty.")
  kdr_csv(attach(x$selection),file.path(tables,"variable_tests_by_repetition.csv"))
  sm<-c("mean_TPR","mean_FPR")
  sm<-sm[vapply(sm,function(m)nrow(sel$overall)>0 && any(is.finite(sel$overall[[m]])),logical(1))]
  if(length(sm))figure("Figure_E_variable_selection","Variable identification across exposure settings",
    "Mean true-positive and false-positive rates compare the two strategies. Only settings with defined false-positive rates appear in that panel.",1,length(sm),function(){
      for(m in sm)sim_E_comparison_panel(sel$overall,m,if(m=="mean_TPR")"Mean true-positive rate" else "Mean false-positive rate",sc)
    })
  if(nrow(sel$variables))figure("Figure_E_active_exposures","Identification of active exposures",
    "Identification rates compare the fixed and adaptive procedures for each of the eight active exposures. Variable-level Wilson intervals and available repetitions are tabulated.",4,2,function(){
      for(v in exposure_names)sim_E_comparison_panel(sel$variables[sel$variables$variable==v,,drop=FALSE],
        "rate",sim_label(v),sc,axis_label="Identification rate")
    })
  # Interpretation tables retain all quantitative ORs, allocation shares and interaction tests.
  for(n in names(x$extended)) {
    d<-attach(x$extended[[n]]);d<-d[,setdiff(names(d),c("M","kappa","lambda")),drop=FALSE]
    emit(d,paste0("Table_E_",n),paste("Configuration selection:",gsub("_"," ",n)),
      "Fixed and adaptive results use common evaluation combinations and references. Contribution shares normalize mean absolute baseline-Shapley values within the eight prespecified active exposures, with all additional coordinates fixed at their reference values. Overall interactions use Holm across four pairs; corrected directions use Holm across all evaluated grid contrasts. Quantitative contrasts use the log-OR scale for error and exponentiation for OR limits. Repetition-specific selected parameters are in Table_E_choices.")
  }
  for(n in names(x$interpretation))kdr_csv(attach(x$interpretation[[n]]),file.path(tables,paste0(n,"_by_repetition.csv")))
  scores<-attach(x$extended$scores)
  if(nrow(scores)) {
    scores$contribution_error_pp<-100*scores$contribution_share_error
    scores$interaction_detection<-100*(1-scores$interaction_miss);scores$direction_detection<-100*(1-scores$direction_miss)
    scores$wrong_direction<-100*scores$direction_wrong
    figure("Figure_E_interpretation","Contribution and interaction recovery",
      "Contribution error uses the mean repetition-level RMSE of percentage shares across the eight prespecified active exposures at the reference background. Detection rates use the specified Holm corrections.",2,2,function(){
        sim_E_comparison_panel(scores,"contribution_error_pp","Contribution share error (pp)",sc)
        sim_E_comparison_panel(scores,"interaction_detection","Overall interaction detection (%)",sc)
        sim_E_comparison_panel(scores,"direction_detection","Correct-direction detection (%)",sc)
        sim_E_comparison_panel(scores,"wrong_direction","Wrong-direction detection (%)",sc)
      })
  }
  choices<-list();rankings<-list()
  for(f in list.files(file.path(run,"choices"),pattern="[.]rds$",full.names=TRUE)) {
    a<-sim_read(f);if(is.null(a))next;c<-a$candidate
    val<-function(n)if(is.null(a[[n]]))NA_real_ else a[[n]]
    choices[[length(choices)+1L]]<-data.frame(scenario_id=a$scenario_id,replicate=a$replicate,
      candidate_id=c$candidate_id,M=c$M,kappa=c$kappa,lambda=c$lambda,
      selection_ok=a$selection_ok,fallback=a$fallback,selection_seconds=val("selection_elapsed_seconds"),
      final_fits_seconds=val("final_fits_elapsed_seconds"),total_seconds=val("total_elapsed_seconds"))
    rankings[[length(rankings)+1L]]<-cbind(scenario_id=a$scenario_id,replicate=a$replicate,a$ranking)
  }
  choices<-kdr_bind(choices);ranking<-kdr_bind(rankings)
  emit(choices,"Table_E_choices","Selected configuration in each independent repetition",
    "Selection success requires at least one candidate with finite validation scores in every fold. Fallback denotes use of the frozen configuration after all candidates failed. The paired final performance includes this prespecified fallback. Elapsed component times reflect checkpoint reuse when a task is resumed.")
  kdr_csv(ranking,file.path(tables,"validation_rankings.csv"))
  if(nrow(choices)) {
    costs<-kdr_bind(lapply(sim_groups(choices,"scenario_id"),function(ix){z<-choices[ix,]
      data.frame(scenario_id=z$scenario_id[1],records=nrow(z),selection_successes=sum(z$selection_ok),
        fallback_count=sum(z$fallback),selection_mean_seconds=sim_mean(z$selection_seconds),
        selection_MCSE_seconds=sim_mcse(z$selection_seconds),final_pair_mean_seconds=sim_mean(z$final_fits_seconds),
        total_pair_mean_seconds=sim_mean(z$total_seconds))}))
    emit(costs,"Table_E_selection_costs","Configuration-selection availability and elapsed time",
      "Elapsed times are within-worker wall times. Selection includes learning and validation; final-pair time covers the paired final computations. Shared work and concurrent execution make their sum different from whole-study wall time. Resumed tasks can have reduced times.")
    if(any(is.finite(choices$selection_seconds)))figure("Figure_E_selection_cost","Configuration-selection computing time",
      "Boxplots show selection elapsed time in minutes across completed independent repetitions. Checkpoint reuse can shorten resumed tasks.",1,1,function(){
        vals<-lapply(sc$scenario_id,function(id)choices$selection_seconds[choices$scenario_id==id]/60)
        good<-vapply(vals,function(v)any(is.finite(v)),logical(1));par(mar=c(8,6,3,1))
        boxplot(vals[good],names=sc$scenario_label[good],las=2,col="#B3DDF2",ylab="Selection time (minutes)")
      })
    freq<-as.data.frame(table(factor(choices$scenario_id,levels=sc$scenario_id),
      factor(choices$candidate_id,levels=x$plan$design$configurations$candidate_id)))
    names(freq)<-c("scenario_id","candidate_id","selected_count")
    denom<-table(choices$scenario_id);freq$percent<-100*freq$selected_count/as.numeric(denom[as.character(freq$scenario_id)])
    emit(freq,"Table_E_selection_frequency","Configuration selection frequencies",
      "Frequencies use available selection records in each scenario and include the explicitly recorded fallback rule. Planned repetitions and availability are tabulated separately.")
    figure("Figure_E_configuration_frequency","Configuration selection frequencies",
      "Cell values are percentages of available selection records. Candidate definitions are in Table_E_candidates.",1,1,function(){
        ids<-x$plan$design$configurations$candidate_id;pal<-colorRampPalette(c("white","#0072B2"))(101)
        par(mar=c(8,5,3,2));plot(c(.5,nrow(sc)+.5),c(.5,length(ids)+.5),type="n",axes=FALSE,xlab="",ylab="Candidate configuration")
        axis(2,seq_along(ids),ids,las=1,cex.axis=.75);axis(1,seq_len(nrow(sc)),sc$scenario_label,las=2,cex.axis=.75)
        for(i in seq_len(nrow(freq)))if(is.finite(freq$percent[i])) {
          xx<-match(freq$scenario_id[i],sc$scenario_id);yy<-match(freq$candidate_id[i],ids);v<-freq$percent[i]
          rect(xx-.5,yy-.5,xx+.5,yy+.5,col=pal[1+round(v)],border="grey90")
          if(v>0)text(xx,yy,format(round(v,1),trim=TRUE),col=if(v>50)"white" else "black",cex=.8)
        }
      })
  }
  for(id in sc$scenario_id) {
    d<-points[points$scenario_id==id,,drop=FALSE];if(!nrow(d))next
    name<-sc$scenario_label[match(id,sc$scenario_id)]
    figure(paste0("Figure_E_recovery_",id),paste("Effect recovery:",name),
      "Dashed black denotes truth or the identity line. Colored curves/points show repetition-mean estimates. These display recovery under the chosen procedure; pointwise interval calibration and repetition variability are tabulated separately.",1,2,function(){
        z<-d[d$target=="joint",];plot(range(z$q),sim_range(c(z$truth,z$mean_estimate)),type="n",
          xlab="Common exposure quantile",ylab="Log odds ratio",main="Joint exposure effect")
        u<-z[z$strategy==unique(z$strategy)[1],];u<-u[order(u$q),];lines(u$q,u$truth,lty=2,lwd=2)
        cols<-c(Fixed="#0072B2",Adaptive="#D55E00")
        for(s in names(cols)){u<-z[z$strategy==s,];u<-u[order(u$q),];lines(u$q,u$mean_estimate,col=cols[s],lwd=2)}
        legend("topleft",c("Truth",names(cols)),col=c("black",cols),lty=c(2,1,1),bty="n",cex=.8)
        z<-d[d$target=="individual",];plot(sim_range(z$truth),sim_range(c(z$mean_estimate,z$truth)),type="n",
          xlab="True log odds ratio",ylab="Mean estimated log odds ratio",main="Observed-combination effect");abline(0,1,lty=2)
        for(s in names(cols)){u<-z[z$strategy==s,];points(u$truth,u$mean_estimate,pch=16,cex=.45,col=adjustcolor(cols[s],alpha.f=.5))}
      })
  }
  supports<-kdr_bind(lapply(sc$scenario_id,function(id){e<-sim_read(file.path(run,"evaluation",paste0(id,".rds")))
    if(is.null(e))NULL else cbind(scenario_id=id,e$support)}))
  emit(supports,"Table_E_joint_support","Geometric support of joint evaluation combinations",
    "Squared Mahalanobis distances use the eight active-exposure coordinates and an independent scenario evaluation sample. The reference is its empirical 95th percentile. This is a geometric description of extrapolation along the common joint path.")
  emit(attach(x$cost),"Table_E_final_fit_resources","Computational resources for final estimation",
    "Per-fit times cover learning, nuisance estimation and final inference components. Configuration selection costs are recorded separately in Table_E_choices. Shared kernel costs are reused. Memory sampling availability depends on the platform.")
  sim_report_model_diagnostics(run,x,cap,"E")
  captions<-kdr_bind(lapply(sort(ls(cap)),function(id){a<-cap[[id]];data.frame(id=id,type=a$type,file=a$file,title=a$title,note=a$note)}))
  kdr_csv(captions,file.path(analysis,"Figure_and_table_captions.csv"))
  writeLines(unlist(lapply(seq_len(nrow(captions)),function(i)c(paste(captions$id[i],captions$title[i]),captions$note[i],""))),
    file.path(analysis,"Figure_and_table_captions.txt"),useBytes=TRUE)
  html<-c("<!doctype html><html lang='en'><meta charset='utf-8'><title>Module E</title><style>body{font-family:serif;margin:30px;max-width:1300px}img{max-width:100%}</style><h1>Data-driven configuration selection</h1>")
  for(i in seq_len(nrow(captions)))html<-c(html,paste0("<h2>",sim_escape(captions$title[i]),"</h2>"),
    if(captions$type[i]=="Figure")paste0("<img src='",captions$file[i],"'>") else paste0("<a href='",captions$file[i],"'>Formatted table</a>"),
    paste0("<p>",sim_escape(captions$note[i]),"</p>"))
  writeLines(c(html,"</html>"),file.path(analysis,"index.html"),useBytes=TRUE)
  dir.create(file.path(analysis,"report_source"),showWarnings=FALSE)
  fs<-list.files(file.path(getOption("kdr.root"),"R"),pattern="^(80|81|85|89|92).*R$",full.names=TRUE)
  file.copy(fs,file.path(analysis,"report_source"),overwrite=TRUE)
  kdr_save(list(report_settings=cfg$report,generated_at=Sys.time()),file.path(analysis,"report_metadata.rds"))
  kdr_log("Module E report saved:",file.path(analysis,"index.html"));invisible(analysis)
}
