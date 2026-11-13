sim_B_response_summary <- function(estimates,grid,replicates,planned) {
  metrics<-points<-reps<-list();tab<-grid$points
  for(ix in sim_groups(tab,c("target","component"))) {
    z<-tab[ix,,drop=FALSE];Y<-estimates[,ix,drop=FALSE]
    usable<-rowSums(!is.finite(Y))==0L;Y<-Y[usable,,drop=FALSE];R<-nrow(Y)
    if(!R)next
    truth<-z$truth;E<-sweep(Y,2,truth);eligible<-z$eligible
    if(!any(eligible))next
    mse<-rowMeans(E[,eligible,drop=FALSE]^2);bias<-colMeans(E)
    jack<-function(v)if(length(v)>1L)sqrt((length(v)-1)/length(v)*sum((v-mean(v))^2)) else NA_real_
    jb<-if(R>1L)vapply(seq_len(R),function(r){
      b<-(colSums(E)-E[r,])/(R-1L);sqrt(mean(b[eligible]^2))
    },numeric(1)) else NA_real_
    jr<-if(R>1L)sqrt(pmax(0,(sum(mse)-mse)/(R-1L))) else NA_real_
    meta<-z[1,c("target","component","exposure_1","exposure_2"),drop=FALSE]
    metrics[[length(metrics)+1L]]<-cbind(meta,data.frame(planned_repetitions=planned,
      usable_repetitions=R,evaluation_points=sum(eligible),
      RMS_Bias=sqrt(mean(bias[eligible]^2)),RMSE=sqrt(mean(mse)),
      MCSE_RMS_Bias=jack(jb),MCSE_RMSE=jack(jr),
      RMS_empirical_SD=sqrt(sim_mean(apply(Y[,eligible,drop=FALSE],2,sim_sd)^2))))
    points[[length(points)+1L]]<-cbind(z,data.frame(usable_repetitions=R,
      mean_estimate=colMeans(Y),bias=bias,pointwise_RMSE=sqrt(colMeans(E^2)),
      empirical_SD=apply(Y,2,sim_sd),
      lower_repetition_quantile=apply(Y,2,quantile,.025),upper_repetition_quantile=apply(Y,2,quantile,.975)))
    reps[[length(reps)+1L]]<-cbind(meta[rep(1,R),,drop=FALSE],data.frame(
      replicate=replicates[usable],RMSE=sqrt(mse),mean_signed_error=rowMeans(E[,eligible,drop=FALSE])))
  }
  list(metrics=kdr_bind(metrics),points=kdr_bind(points),repetitions=kdr_bind(reps))
}

sim_B_response_collect <- function(run,cfg,pp,settings) {
  plan<-pp$tasks;p<-unique(plan$p);stopifnot(length(p)==1L)
  ev<-readRDS(sim_evaluation_path(run,p,cfg$report$correlation_id))
  grid<-ev$interpretation$response
  records<-lapply(plan$key,function(k)sim_read(file.path(run,"fits",paste0(k,".rds"))))
  status<-plan[c("key","config_id","method","scenario_gf","replicate")]
  status$prediction_ok<-vapply(records,function(a)isTRUE(a$interpretation_ok),logical(1))
  status$error<-vapply(records,function(a)if(is.null(a))"Pending" else if(is.null(a$interpretation_error))a$error else a$interpretation_error,character(1))
  metrics<-points<-reps<-list()
  keys<-c("config_id","method","scenario_gf")
  for(ix in sim_groups(plan,keys)) {
    ids<-ix[vapply(records[ix],function(a)isTRUE(a$interpretation_ok),logical(1))]
    if(!length(ids))next
    Y<-do.call(rbind,lapply(records[ids],function(a)a$interpretation$response$estimate))
    q<-sim_B_response_summary(Y,grid,plan$replicate[ids],length(ix))
    lower<-do.call(rbind,lapply(records[ids],function(a)a$interpretation$response$lower))
    upper<-do.call(rbind,lapply(records[ids],function(a)a$interpretation$response$upper))
    q$points$mean_CI_lower<-apply(lower,2,sim_mean)[q$points$index]
    q$points$mean_CI_upper<-apply(upper,2,sim_mean)[q$points$index]
    q$points$interval_repetitions<-colSums(is.finite(lower)&is.finite(upper))[q$points$index]
    meta<-plan[ix[1],keys,drop=FALSE]
    add<-function(z)if(nrow(z))cbind(meta[rep(1,nrow(z)),,drop=FALSE],z) else NULL
    metrics[[length(metrics)+1L]]<-add(q$metrics);points[[length(points)+1L]]<-add(q$points);reps[[length(reps)+1L]]<-add(q$repetitions)
  }
  list(grid=grid,plan=plan,status=status,metrics=kdr_bind(metrics),points=kdr_bind(points),repetitions=kdr_bind(reps))
}

sim_B_display_scenarios <- function(x) {
  if(!nrow(x))return(x)
  kdr_bind(lapply(paste0("G",1:4),function(g){
    sg<-if(g%in%c("G1","G2"))"correct" else "wrong"
    z<-x[(x$method=="KDR"&x$scenario_gf==g)|(x$method!="KDR"&x$scenario_gf==sg),,drop=FALSE]
    z$fit_scenario<-z$scenario_gf;z$scenario_gf<-rep(g,nrow(z));z
  }))
}

sim_B_exposure_axis <- function(variable) {
  switch(variable,PM25=expression(PM[2.5]~(mu*g/m^3)),PMcoarse=expression(PM[2.5-10]~(mu*g/m^3)),
    O3=expression(O[3]~(mu*g/m^3)),NO2=expression(NO[2]~(mu*g/m^3)),SO2=expression(SO[2]~(mu*g/m^3)),
    CO=expression(CO~(mg/m^3)),temperature=expression(Temperature~(degree*C)),humidity="Relative humidity (%)")
}

sim_B_response_figure <- function(path,title,nr,nc,np,cfg,methods,draw,interaction=FALSE) {
  cf<-cfg;cf$report$width_inches<-max(cf$report$width_inches,4.8*nc)
  colors<-c(KDR="#0072B2",SplineCLR="#D55E00",LinearCLR="#009E73")
  sim_figure(path,title,nr,nc,cf,function(){
    # A shared legend has its own strip; the last incomplete row spans its cells.
    panels<-matrix(c(seq_len(np)+1L,rep(0L,nr*nc-np)),nr,nc,byrow=TRUE)
    if(any(panels[nr,]==0L))panels[nr,panels[nr,]==0L]<-max(panels)
    layout(rbind(rep(1L,nc),panels),heights=c(.65,rep(4.2,nr)))
    par(cex=cf$report$base_cex,cex.axis=1.02,cex.lab=1.15,cex.main=1.14)
    old<-par("mar");par(mar=c(0,0,0,0));plot.new()
    legend("center",c(if(interaction)"Identity line" else "Truth",sim_label(methods)),
      col=c("black",colors[methods]),lty=c(2,if(interaction)rep(NA,length(methods)) else rep(1,length(methods))),
      pch=c(NA,if(interaction)c(KDR=16,SplineCLR=17,LinearCLR=1)[methods] else rep(NA,length(methods))),
      lwd=2,horiz=TRUE,bty="n",cex=1.02)
    if(!interaction)mtext("Light shading: mean pointwise 95% CI limits    Dark shading: 95% repetition range",side=1,line=-.7,cex=.78)
    par(mar=old);draw(colors)
  })
}

sim_plot_B_responses <- function(out,cfg,folder,captions) {
  d<-sim_B_display_scenarios(out$points);if(!nrow(d))return(invisible(NULL))
  methods<-c("KDR","SplineCLR","LinearCLR");settings<-out$grid$settings
  scenarios<-paste0("G",1:4)
  base_note<-paste("All other exposures are fixed at the saved marginal-median reference.",
    "The evaluation range for each exposure is the range of its saved joint-curve concentrations.",
    "These constructed combinations define fixed-background response targets;",
    "some combinations can have limited data support.",
    "The same saved fits and reference vector are used in all comparisons.",
    "CLR uses the outcome adjustment: correct-adjustment fits supply G1/G2 and misspecified-adjustment fits supply G3/G4.")
  for(g in scenarios) {
    z<-d[d$target=="single"&d$scenario_gf==g,,drop=FALSE]
    active<-exposure_names[exposure_names%in%z$component];if(!length(active))next
    prefix<-if(is.null(cfg$report$response_prefix))"B" else cfg$report$response_prefix
    id<-paste0("Figure_",prefix,"_single_",g);title<-paste("Single-exposure effect recovery:",g)
    ms<-methods[methods%in%z$method];nc<-min(4L,length(active));nr<-ceiling(length(active)/nc)
    sim_B_response_figure(file.path(folder,paste0(id,".png")),title,nr,nc,length(active),cfg,ms,function(cols){
      for(j in seq_along(active)) {
        a<-z[z$component==active[j],,drop=FALSE];truth<-a[a$method==a$method[1],];truth<-truth[order(truth$concentration_1),]
        # Consistent scales for a given exposure across G1-G4.
        all<-d[d$target=="single"&d$component==active[j],,drop=FALSE]
        yy<-c(all$truth,all$mean_estimate,all$mean_CI_lower,all$mean_CI_upper)
        if(settings$B_repetition_ranges)yy<-c(yy,all$lower_repetition_quantile,all$upper_repetition_quantile)
        plot(range(truth$concentration_1),sim_range(yy,0),type="n",xlab=sim_B_exposure_axis(active[j]),
          ylab="Log odds ratio",main=paste0(LETTERS[j],"  ",sim_label(active[j])))
        abline(h=0,lty=3,col="grey75")
        for(m in ms) {
          b<-a[a$method==m,,drop=FALSE];b<-b[order(b$concentration_1),]
          good<-is.finite(b$mean_CI_lower)&is.finite(b$mean_CI_upper)
          if(sum(good)>1L)polygon(c(b$concentration_1[good],rev(b$concentration_1[good])),
            c(b$mean_CI_lower[good],rev(b$mean_CI_upper[good])),col=adjustcolor(cols[m],.10),border=NA)
        }
        if(settings$B_repetition_ranges)for(m in ms) {
          b<-a[a$method==m,,drop=FALSE];b<-b[order(b$concentration_1),]
          good<-is.finite(b$lower_repetition_quantile)&is.finite(b$upper_repetition_quantile)
          if(sum(good)>1L&&max(b$usable_repetitions)>1L)polygon(c(b$concentration_1[good],rev(b$concentration_1[good])),
            c(b$lower_repetition_quantile[good],rev(b$upper_repetition_quantile[good])),col=adjustcolor(cols[m],.22),border=NA)
        }
        for(m in ms) {
          b<-a[a$method==m,,drop=FALSE];b<-b[order(b$concentration_1),]
          if(nrow(b))lines(b$concentration_1,b$mean_estimate,col=cols[m],lwd=2.2)
        }
        lines(truth$concentration_1,truth$truth,lty=2,col="black",lwd=2.1)
      }
    })
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note=paste(base_note,"Each contrast is h(a, reference[-j]) - h(reference), on the log OR scale.",
        "This background-specific effect includes interaction contributions at the specified background.",
        "Colored lines are means over usable independent repetitions. Interval availability is tabulated at each point.",
        if(settings$B_repetition_ranges)"Light shading averages the pointwise 95% interval lower and upper limits separately; darker shading is the empirical 2.5th to 97.5th percentile range of point estimates. These show typical model uncertainty and between-repetition variability, respectively." else "Repetition variability is reported in the tables.",
        "Axes are identical for each exposure across the four scenario figures. Reference contrasts equal zero by construction and are excluded from aggregate errors.",
        "Humidity has no main effect; its median-background response can be small despite its temperature interaction."))
  }
  z<-d[d$target=="interaction",,drop=FALSE];pairs<-unique(z$component)
  if(!length(pairs))return(invisible(NULL))
  gs<-scenarios[scenarios%in%z$scenario_gf]
  panels<-list()
  for(pair in pairs)for(g in gs) {
    a<-z[z$component==pair&z$scenario_gf==g&z$eligible,,drop=FALSE]
    if(nrow(a))panels[[length(panels)+1L]]<-a
  }
  if(!length(panels))return(invisible(NULL))
  prefix<-if(is.null(cfg$report$response_prefix))"B" else cfg$report$response_prefix
  id<-paste0("Figure_",prefix,"_interactions");title<-"Interaction-effect recovery"
  ms<-methods[methods%in%z$method];nc<-length(gs);nr<-ceiling(length(panels)/nc)
  sim_B_response_figure(file.path(folder,paste0(id,".png")),title,nr,nc,length(panels),cfg,ms,function(cols){
    for(i in seq_along(panels)) {
      a<-panels[[i]];all<-z[z$component==a$component[1]&z$eligible,,drop=FALSE]
      yy<-c(all$truth,all$mean_estimate)
      if(settings$B_repetition_ranges)yy<-c(yy,all$lower_repetition_quantile,all$upper_repetition_quantile)
      lim<-sim_range(yy,0)
      plot(lim,lim,type="n",xlab="True interaction (log OR)",ylab="Estimated interaction (log OR)",
        main=paste(sim_label(a$exposure_1[1]),"x",sim_label(a$exposure_2[1]),"|",a$scenario_gf[1]),asp=1)
      abline(h=0,v=0,lty=3,col="grey85");abline(0,1,lty=2,col="black",lwd=1.6)
      for(m in ms) {
        b<-a[a$method==m,,drop=FALSE];if(!nrow(b))next
        if(settings$B_repetition_ranges) {
          good<-is.finite(b$lower_repetition_quantile)&is.finite(b$upper_repetition_quantile)
          segments(b$truth[good],b$lower_repetition_quantile[good],b$truth[good],b$upper_repetition_quantile[good],
            col=adjustcolor(cols[m],.12))
        }
        points(b$truth,b$mean_estimate,pch=c(KDR=16,SplineCLR=17,LinearCLR=1)[m],col=adjustcolor(cols[m],.68),cex=.72)
      }
    }
  },interaction=TRUE)
  captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
    note=paste(base_note,"Rows identify exposure pairs and columns identify nuisance scenarios.",
      "Interaction is h(a,b,z0) - h(a,b0,z0) - h(a0,b,z0) + h(a0,b0,z0): departure from additivity on the log OR scale.",
      "Both concentrations vary independently over a rectangular grid. Points show true contrasts versus mean predictions across repetitions.",
      if(settings$B_repetition_ranges)"Faint vertical segments show empirical central 95% repetition ranges, showing between-repetition variability." else "Repetition variability is reported in the tables.",
      "Both CLR exposure models are additive and predict zero interaction up to numerical precision; their points can overlap.",
      "This structural difference is distinct from nuisance-model robustness. Anchor rows and columns, which have identically zero contrast, are excluded from the scatter and aggregate errors.",
      "All remaining points, including genuine zero interactions, enter error summaries. Scales are common across scenarios for each pair."))
  invisible(NULL)
}

sim_report_B_responses <- function(run,result,tables,figures,captions) {
  cfg<-result$config;s<-sim_B_response_settings(cfg$report,cfg)
  if(!s$B_response_summaries)return(invisible(NULL))
  out<-sim_B_response_collect(run,cfg,result$plan,s)
  analysis<-dirname(tables)
  kdr_save(out,file.path(analysis,"B_response_predictions.rds"))
  kdr_csv(out$grid$points,file.path(tables,"B_response_evaluation_points.csv"))
  kdr_csv(out$grid$concentration_ranges,file.path(tables,"B_response_concentration_ranges.csv"))
  emit<-function(d,id,title,note)sim_export_table(d,id,title,note,tables,captions)
  emit(out$status,"Table_B_response_availability","Availability of fixed-background predictions",
    "Each row identifies one planned saved fit. Predictions require its saved coefficients and kernel or comparator basis. Pending fits, failed fits and missing model files remain explicit. Comparator fits are evaluated once per adjustment specification and repetition, then displayed in the corresponding G scenarios.")
  m<-sim_B_display_scenarios(out$metrics)
  for(target in c("single","interaction")) {
    a<-m[m$target==target,,drop=FALSE];if(!nrow(a))next
    emit(a,paste0("Table_B_",target,"_accuracy"),
      if(target=="single")"Accuracy of single-exposure effect curves" else "Accuracy of interaction-effect contrasts",
      paste("Contrasts are on the log OR scale at a common marginal-median background.",
        "RMS bias is sqrt(mean over evaluation points of squared repetition-mean error). RMSE is sqrt(mean squared error over all usable repetitions and points).",
        "MCSE uses deletion of an entire repetition; the empirical SD is summarized on its root-mean-square scale.",
        "Identically zero anchor contrasts are excluded; other zero true effects enter the absolute-error summaries.",
        "Usable/planned counts are reported separately for each component. Independent repetitions determine Monte Carlo uncertainty.",
        "Fit scenario records the CLR adjustment fit and its original independent-repetition count.",
        if(target=="interaction")"Interaction is defined by four function evaluations; additive CLR models have zero log-OR interaction by construction." else
          "A single-exposure curve describes a fixed-background contrast, including any interactions at that background, rather than a unique main-effect decomposition."))
  }
  if(nrow(out$repetitions))kdr_csv(sim_B_display_scenarios(out$repetitions),file.path(tables,"B_response_metrics_by_repetition.csv"))
  if(isTRUE(cfg$report$save_pointwise_tables)&&nrow(out$points))
    kdr_csv(sim_B_display_scenarios(out$points),file.path(tables,"B_response_pointwise_summary.csv"))
  sim_plot_B_responses(out,cfg,figures,captions)
  kdr_log("Fixed-background response summaries saved.")
  invisible(out)
}
