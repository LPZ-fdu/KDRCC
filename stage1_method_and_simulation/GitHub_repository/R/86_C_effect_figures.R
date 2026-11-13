# Reporting only: all quantities are read from saved fits; no model is refitted.
sim_panel_layout <- function(n) {
  if(n>1L && n%%2L==1L)
    layout(matrix(c(seq_len(n-1L),n,n),ncol=2L,byrow=TRUE))
}
sim_C_box_key <- function(config_id,method,scenario,target,null_count)
  paste(config_id,method,scenario,target,null_count,sep="|")

sim_plot_C_effects <- function(result,folder,captions) {
  cfg<-result$config;d<-result$points;m<-result$effects
  if(!nrow(d)||!nrow(m))return(invisible(NULL))
  wanted<-cfg$report$C_null_counts
  if(!is.null(wanted)) {
    d<-d[d$null_count%in%wanted,,drop=FALSE]
    m<-m[m$null_count%in%wanted,,drop=FALSE]
  }
  if(!nrow(d)||!nrow(m))return(invisible(NULL))
  per_page<-cfg$report$C_panels_per_page
  if(is.null(per_page))per_page<-6L
  if(length(per_page)!=1L||!is.finite(per_page)||per_page<1||per_page!=floor(per_page))
    stop("C_panels_per_page must be a positive integer.")
  blue<-"#0072B2";orange<-"#D55E00"
  caption<-function(id,title,note) {
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,note=note)
  }
  base_note<-paste("KDR-CC with both nuisance working families correctly specified (G1).",
    "Each model includes all eight active exposures and the indicated number of null exposures;",
    "effects are evaluated from the full candidate model.",
    "Usable/planned repetitions and effective configurations are reported in the tables.")
  for(target in c("joint","individual")) {
    z<-d[d$target==target,,drop=FALSE]
    z<-z[is.finite(z$truth)&is.finite(z$mean_estimate),,drop=FALSE]
    if(!nrow(z))next
    qs<-sort(unique(z$null_count))
    sets<-lapply(qs,function(q){
      a<-z[z$null_count==q,,drop=FALSE]
      if(target=="joint") {
        if(!any(abs(a$q-.5)<1e-12)) {
          ref<-a[1,,drop=FALSE];ref$q<-.5
          for(v in c("truth","mean_estimate","lower_repetition_quantile","upper_repetition_quantile"))ref[[v]]<-0
          a<-rbind(a,ref)
        }
        a<-a[order(a$q),,drop=FALSE]
      } else a<-a[order(a$index),,drop=FALSE]
      a
    })
    all_y<-if(target=="joint")unlist(lapply(sets,function(a)c(a$truth,a$mean_estimate,
      a$lower_repetition_quantile,a$upper_repetition_quantile))) else c(z$truth,z$mean_estimate)
    lim<-sim_range(all_y)
    if(target=="joint")lim[2]<-lim[2]+.15*diff(lim)
    pages<-split(seq_along(qs),ceiling(seq_along(qs)/per_page))
    for(pg in seq_along(pages)) {
      ids<-pages[[pg]];nc<-min(2L,length(ids));nr<-ceiling(length(ids)/nc)
      id<-sprintf("Figure_C_%s_page%02d",target,pg)
      title<-if(target=="joint")"Joint-effect recovery as null exposures increase" else
        "Observed-combination effect recovery as null exposures increase"
      sim_figure(file.path(folder,paste0(id,".png")),title,nr,nc,cfg,function(){
        sim_panel_layout(length(ids))
        for(ii in seq_along(ids)) {
          k<-ids[ii];a<-sets[[k]]
          main<-paste0(LETTERS[ii],"  ",qs[k]," null exposures")
          if(target=="joint") {
            plot(range(a$q),lim,type="n",xlab="Common exposure quantile",ylab="Log odds ratio",main=main)
            abline(h=0,lty=3,col="grey75")
            good<-is.finite(a$lower_repetition_quantile)&is.finite(a$upper_repetition_quantile)
            if(sum(good)>1L)polygon(c(a$q[good],rev(a$q[good])),
              c(a$lower_repetition_quantile[good],rev(a$upper_repetition_quantile[good])),
              col=adjustcolor(blue,alpha.f=.13),border=NA)
            lines(a$q,a$mean_estimate,col=blue,lwd=2.3)
            lines(a$q,a$truth,col="black",lty=2,lwd=2.1)
            legend("topleft",c("Truth","Mean KDR-CC estimate"),col=c("black",blue),
              lty=c(2,1),lwd=2,bty="n",cex=.86)
          } else {
            plot(lim,lim,type="n",xlab="True log odds ratio",ylab="Mean estimated log odds ratio",main=main)
            abline(0,1,lty=2,lwd=2,col="black")
            keep<-unique(round(seq(1,nrow(a),length.out=min(cfg$report$scatter_points,nrow(a)))))
            points(a$truth[keep],a$mean_estimate[keep],pch=16,cex=.65,col=adjustcolor(blue,alpha.f=.5))
            legend("topleft",c("Identity line","Mean KDR-CC estimate"),col=c("black",blue),
              lty=c(2,NA),pch=c(NA,16),bty="n",cex=.86)
          }
        }
      })
      note<-if(target=="joint")paste("All candidate exposures follow a common marginal quantile relative to their median reference.",
        "Null exposures have zero true response contribution and are included in estimation.",
        "Lines are repetition means; shading is the empirical central 95% range across repetitions, showing between-repetition variability.",
        "The median-reference point is exactly zero. Identical axes are used across pages.") else
        paste("Each point is a fixed independent exposure combination relative to the median reference.",
        "The vertical coordinate is its mean estimate across repetitions. The dashed identity line denotes perfect agreement.",
        "A deterministic subset of points is displayed, with identical axes across pages.",
        "Repetition-level RMSE, empirical SD and interval calibration describe individual-fit reliability alongside the mean curve.")
      caption(id,title,paste(base_note,note,"Page",pg,"of",length(pages),"."))
    }
  }
  # All null counts remain in tables even if a plotting subset was requested.
  metric_names<-c("RMS_Bias","RMSE","mean_OR_absolute_relative_error","correct_significant_percent")
  panels<-list()
  for(metric in metric_names)for(target in c("joint","individual")) {
    a<-m[m$target==target&is.finite(m[[metric]]),,drop=FALSE]
    if(nrow(a))panels[[length(panels)+1L]]<-list(a=a,metric=metric,target=target)
  }
  if(length(panels)) {
    id<-"Figure_C_effect_accuracy";title<-"Effect-estimation accuracy as null exposures increase"
    sim_figure(file.path(folder,paste0(id,".png")),title,ceiling(length(panels)/2),min(2,length(panels)),cfg,function(){
      sim_panel_layout(length(panels))
      for(i in seq_along(panels)) {
        pp<-panels[[i]];a<-pp$a;a<-a[order(a$null_count),];metric<-pp$metric
        mc<-switch(metric,RMS_Bias="MCSE_RMS_Bias",RMSE="MCSE_RMSE",
          mean_OR_absolute_relative_error="MCSE_mean_OR_error",correct_significant_percent="MCSE_correct_significant_percent")
        y<-a[[metric]];lo<-pmax(0,y-1.96*a[[mc]]);hi<-y+1.96*a[[mc]]
        if(metric=="correct_significant_percent")hi<-pmin(100,hi)
        plot(a$null_count,y,type="b",pch=16,col=blue,lwd=2,xlab="Number of null exposures",
          ylab=sim_metric_label(metric),ylim=if(metric=="correct_significant_percent")c(0,100) else sim_range(c(y,lo,hi),0,TRUE),
          main=paste0(LETTERS[i],"  ",sim_label(pp$target)))
        good<-is.finite(lo)&is.finite(hi)&hi>lo
        if(any(good))arrows(a$null_count[good],lo[good],a$null_count[good],hi[good],angle=90,code=3,length=.035,col=blue)
      }
    })
    caption(id,title,paste(base_note,"Rows show RMS bias, RMSE, mean absolute OR relative error and correct-direction significant detection.",
      "Bars are approximate 95% Monte Carlo intervals (estimate plus/minus 1.96 MCSE); bounds respect the metric range.",
      "Bias and RMSE MCSE use deletion of an entire repetition. Other MCSE use repetition-level summaries.",
      "Exact zero contrasts are excluded from detection denominators. Bars describe Monte Carlo uncertainty across independent repetitions."))
  }
  panels<-list()
  for(metric in c("dispersion","SE_ratio","coverage_percent"))for(target in c("joint","individual")) {
    a<-m[m$target==target,,drop=FALSE]
    good<-if(metric=="dispersion")is.finite(a$RMS_empirical_SD)|is.finite(a$RMS_model_SE) else is.finite(a[[metric]])
    if(any(good))panels[[length(panels)+1L]]<-list(a=a[good,,drop=FALSE],metric=metric,target=target)
  }
  if(length(panels)) {
    id<-"Figure_C_uncertainty";title<-"Uncertainty calibration as null exposures increase"
    sim_figure(file.path(folder,paste0(id,".png")),title,ceiling(length(panels)/2),min(2,length(panels)),cfg,function(){
      sim_panel_layout(length(panels))
      for(i in seq_along(panels)) {
        pp<-panels[[i]];a<-pp$a;a<-a[order(a$null_count),];metric<-pp$metric
        if(metric=="dispersion") {
          lim<-sim_range(c(a$RMS_empirical_SD,a$RMS_model_SE),0,TRUE);lim[2]<-lim[2]+.2*diff(lim)
          plot(range(a$null_count),lim,type="n",xlab="Number of null exposures",ylab="SD / SE (log OR)",
            main=paste0(LETTERS[i],"  ",sim_label(pp$target)))
          lines(a$null_count,a$RMS_empirical_SD,col=blue,lwd=2,type="b",pch=16)
          lines(a$null_count,a$RMS_model_SE,col=orange,lwd=2,type="b",pch=17)
          legend("topleft",c("Empirical SD (RMS)","Model SE (RMS)"),col=c(blue,orange),pch=c(16,17),lty=1,bty="n",cex=.85)
        } else {
          ref<-if(metric=="SE_ratio")1 else 100*cfg$evaluation$confidence
          plot(a$null_count,a[[metric]],type="b",pch=16,col=blue,lwd=2,xlab="Number of null exposures",
            ylab=sim_metric_label(metric),ylim=if(metric=="coverage_percent")c(0,100) else sim_range(a[[metric]],ref,TRUE),
            main=paste0(LETTERS[i],"  ",sim_label(pp$target)))
          abline(h=ref,lty=3,col="grey45")
        }
      }
    })
    caption(id,title,paste(base_note,"Empirical SD and model SE are summarized by their root mean square at common interval-evaluation points.",
      "The reference SE ratio is one. Coverage is a secondary assessment of the nominal",100*cfg$evaluation$confidence,"percent pointwise intervals.",
      "Intervals condition on the learned kernel. Their empirical calibration is assessed alongside finite-representation error and regularization bias; widths are supplied in the tables."))
  }
  panels<-list()
  for(target in c("joint","individual")) {
    a<-m[m$target==target,,drop=FALSE];a<-a[order(a$null_count),];series<-list()
    for(i in seq_len(nrow(a))) {
      key<-sim_C_box_key(a$config_id[i],a$method[i],a$scenario_gf[i],target,a$null_count[i])
      v<-result$or_errors[[key]];v<-v[is.finite(v)]
      if(length(v))series[[as.character(a$null_count[i])]]<-v
    }
    if(length(series))panels[[length(panels)+1L]]<-list(target=target,series=series)
  }
  if(length(panels)) {
    id<-"Figure_C_OR_errors";title<-"Absolute OR relative errors as null exposures increase"
    sim_figure(file.path(folder,paste0(id,".png")),title,length(panels),1,cfg,function(){
      for(i in seq_along(panels)) {
        pp<-panels[[i]]
        boxplot(pp$series,outline=FALSE,col=adjustcolor(blue,.22),cex.axis=.92,
          xlab="Number of null exposures",ylab="Absolute OR relative error (%)",
          main=paste0(LETTERS[i],"  ",sim_label(pp$target)))
      }
    })
    caption(id,title,paste(base_note,"Error is 100 * abs(exp(estimated log OR - true log OR) - 1).",
      "Boxes pool finite point-by-repetition errors; whiskers use 1.5 IQR. Outliers are omitted from display only.",
      "All errors, including infinite overflow values, remain in numeric summaries. Means, medians, quartiles and 95th percentiles are tabulated.",
      "Monte Carlo uncertainty uses independent repetitions as the sampling unit."))
  }
  invisible(NULL)
}
