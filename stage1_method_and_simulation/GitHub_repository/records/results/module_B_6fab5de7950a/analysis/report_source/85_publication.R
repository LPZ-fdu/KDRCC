sim_escape <- function(x) {
  x<-gsub("&","&amp;",as.character(x),fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE)
  gsub(">","&gt;",x,fixed=TRUE)
}
sim_html_table <- function(d,title,note,path) {
  if(!nrow(d))return(invisible(NULL))
  format_cell<-function(x) {
    if(!is.numeric(x))return(sim_escape(x))
    vapply(x,function(a)if(is.na(a))"NA" else if(!is.finite(a))as.character(a) else
      if(a!=0&&abs(a)<.0001)format(a,digits=3,scientific=TRUE) else format(round(a,4),trim=TRUE,scientific=FALSE),character(1))
  }
  body<-as.data.frame(lapply(d,format_cell),stringsAsFactors=FALSE,check.names=FALSE)
  headings<-gsub("_"," ",names(d),fixed=TRUE)
  labels<-c(config_id="Configuration",scenario_gf="Scenario",method="Method",target="Target",
    available_effect_repetitions="Usable repetitions",available_interval_repetitions="Usable intervals",
    RMS_Bias="RMS bias",RMSE="RMSE",MCSE_RMS_Bias="MCSE: RMS bias",MCSE_RMSE="MCSE: RMSE",
    correct_significant_percent="Correct-direction detection (%)",wrong_significant_percent="Wrong-direction detection (%)",
    MCSE_correct_significant_percent="Detection MCSE (%)",mean_OR_absolute_relative_error="Mean absolute OR error (%)",
    median_OR_absolute_relative_error="Median absolute OR error (%)",OR_error_Q25="OR error Q25 (%)",
    OR_error_Q75="OR error Q75 (%)",OR_error_P95="OR error P95 (%)",MCSE_mean_OR_error="Mean OR-error MCSE (%)",
    RMS_empirical_SD="Empirical SD (RMS)",RMS_model_SE="Model SE (RMS)",SE_ratio="SE / empirical SD",
    coverage_percent="Coverage (%)",MCSE_coverage_percent="Coverage MCSE (%)",mean_interval_width="Mean interval width",
    null_count="Null exposures",p="Total exposures",mean_TPR="Mean TPR",mean_FPR="Mean FPR")
  ix<-names(d)%in%names(labels);headings[ix]<-unname(labels[names(d)[ix]])
  if("method"%in%names(body))body$method<-sim_label(body$method)
  if("target"%in%names(body))body$target<-sim_label(body$target)
  rows<-apply(body,1,function(z)paste0("<tr>",paste0("<td>",z,"</td>",collapse=""),"</tr>"))
  text<-c("<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<style>body{font-family:'Times New Roman',serif;font-size:12pt;margin:28px}table{border-collapse:collapse;font-size:11pt}th,td{padding:7px 10px;text-align:right}th{border-top:2px solid #111;border-bottom:1px solid #111}tr:last-child td{border-bottom:2px solid #111}th:first-child,td:first-child{text-align:left}.note{max-width:1100px;line-height:1.4}h2{font-size:14pt}</style></head><body>",
    paste0("<h2>",sim_escape(title),"</h2><table><thead><tr>",paste0("<th>",sim_escape(headings),"</th>",collapse=""),"</tr></thead><tbody>"),
    rows,paste0("</tbody></table><p class='note'><b>Note.</b> ",sim_escape(note),"</p></body></html>"))
  writeLines(text,path,useBytes=TRUE);invisible(path)
}
sim_export_table <- function(d,id,title,note,folder,captions) {
  if(!nrow(d))return(invisible(NULL))
  kdr_csv(d,file.path(folder,paste0(id,".csv")))
  sim_html_table(d,title,note,file.path(folder,paste0(id,".html")))
  captions[[id]]<-list(type="Table",file=paste0("tables/",id,".html"),title=title,note=note)
  invisible(NULL)
}
sim_report_font <- function(report) {
  requested<-report$font;effective<-requested
  if(.Platform$OS.type!="windows" && !requested%in%c("serif","sans","mono")) {
    fc<-Sys.which("fc-match")
    family<-if(nzchar(fc))tryCatch(system2(fc,c("-f",shQuote("%{family}"),shQuote(requested)),
      stdout=TRUE,stderr=FALSE),error=function(e)character()) else character()
    names<-trimws(unlist(strsplit(paste(family,collapse=","),",",fixed=TRUE)))
    if(!tolower(requested)%in%tolower(names)) {
      effective<-"serif"
      kdr_log("Requested font",requested,"is unavailable; using serif. The effective font is recorded in report metadata.")
    }
  }
  report$font_requested<-requested;report$font<-effective;report
}
sim_figure <- function(path,title,nrow=1,ncol=1,cfg,draw) {
  height<-max(5.6,4.2*nrow+1.3)
  png(path,width=cfg$report$width_inches,height=height,units="in",res=cfg$report$dpi,type="cairo")
  device<-dev.cur();on.exit(if(device%in%dev.list())dev.off(device),add=TRUE)
  font<-cfg$report$font
  par(mfrow=c(nrow,ncol),mar=c(5.8,6.6,3.2,1.4),oma=c(.8,.6,3.4,.6),
      family=font,cex=cfg$report$base_cex,cex.axis=1.02,cex.lab=1.15,cex.main=1.14,
      las=1,bty="l",mgp=c(3.8,1,0))
  draw();mtext(title,outer=TRUE,line=.9,font=2,cex=1.45)
  invisible(path)
}
sim_range <- function(x,reference=NULL,nonnegative=FALSE) {
  x<-c(x,reference);x<-x[is.finite(x)]
  if(!length(x))return(c(0,1))
  r<-range(x);w<-diff(r);if(w==0)w<-max(abs(r),1)*.1
  r<-r+c(-.06,.10)*w;if(nonnegative)r[1]<-max(0,r[1]);r
}
sim_metric_label <- function(n) {
  map<-c(RMS_Bias="RMS bias (log OR)",RMSE="RMSE (log OR)",
    mean_OR_absolute_relative_error="Mean absolute OR error (%)",
    correct_significant_percent="Correct-direction detection (%)",SE_ratio="Model SE / empirical SD",
    coverage_percent="Pointwise coverage (%)",component_total_seconds="Component time (seconds)",
    sampled_peak_memory_GB="Sampled peak memory (GB)",cache_GB="Working cache (GB)")
  v<-unname(map[n]);if(is.na(v))n else v
}
sim_lines_panel <- function(d,x,y,group,title,xlab,reference=NULL) {
  d<-d[is.finite(d[[x]])&is.finite(d[[y]]),,drop=FALSE];if(!nrow(d))return(invisible(NULL))
  groups<-unique(d[[group]]);cols<-c("#0072B2","#D55E00","#009E73","#CC79A7")
  limits<-sim_range(d[[y]],reference)
  if(length(groups)>1)limits[2]<-limits[2]+.22*diff(limits)
  plot(range(d[[x]]),limits,type="n",xlab=xlab,ylab=sim_metric_label(y),main=title)
  if(!is.null(reference))abline(h=reference,lty=3,col="grey55")
  for(i in seq_along(groups)) {
    z<-d[d[[group]]==groups[i],];z<-z[order(z[[x]]),]
    lines(z[[x]],z[[y]],col=cols[i],lwd=2.3);points(z[[x]],z[[y]],col=cols[i],pch=16,cex=.7)
  }
  if(length(groups)>1)legend("topright",legend=sim_label(groups),col=cols[seq_along(groups)],lty=1,lwd=2,bty="n",cex=.86)
}

sim_plot_A <- function(result,folder,captions) {
  cfg<-result$config;d<-result$effects;views<-result$plan$design$membership
  if(!nrow(d)||!nrow(views))return(invisible(NULL))
  cost<-result$cost;resources<-data.frame()
  if(nrow(cost))resources<-aggregate(cost[c("component_total_seconds","sampled_peak_memory_GB","cache_GB")],
    by=cost[c("config_id",sim_parameter_names)],FUN=sim_mean)
  metric_names<-c("RMS_Bias","RMSE","mean_OR_absolute_relative_error","correct_significant_percent","SE_ratio")
  for(v in intersect(sim_parameter_names,unique(views$view))) {
    ids<-views$config_id[views$view==v];z<-d[d$config_id%in%ids,]
    r<-resources[resources$config_id%in%ids,,drop=FALSE]
    panels<-lapply(metric_names,function(m)list(data=z,metric=m,group="target"))
    for(m in c("component_total_seconds","sampled_peak_memory_GB","cache_GB")) {
      rr<-r;rr$target<-"Resource use";panels[[length(panels)+1L]]<-list(data=rr,metric=m,group="target")
    }
    panels<-Filter(function(p)nrow(p$data)&&any(is.finite(p$data[[p$metric]])),panels)
    if(!length(panels))next
    id<-paste0("Figure_A_",v);title<-paste("Performance across",if(v=="N")"estimation sample sizes" else if(v=="L")"learning sample sizes" else if(v=="M")"numbers of kernel centers" else paste(v,"values"))
    ncol<-if(length(panels)>1)2L else 1L;nrow<-ceiling(length(panels)/ncol)
    sim_figure(file.path(folder,paste0(id,".png")),title,nrow,ncol,cfg,function(){
      # A custom layout avoids allocating empty panels when a metric is unavailable.
      if(length(panels)%%2L==1L && length(panels)>1L)layout(matrix(c(seq_len(length(panels)-1L),length(panels),length(panels)),ncol=2,byrow=TRUE))
      for(i in seq_along(panels)) {
        p<-panels[[i]];ref<-if(p$metric=="SE_ratio")1 else NULL
        xlab<-switch(v,N="Estimation matched sets (N)",L="Learning matched sets (L)",M="Kernel centers (M)",kappa="Kernel-learning penalty (kappa)",lambda="Interaction penalty (lambda)")
        if(v=="kappa") {
          p$data$log10_kappa<-log10(p$data$kappa)
          sim_lines_panel(p$data,"log10_kappa",p$metric,p$group,paste0(LETTERS[i],"  ",sim_metric_label(p$metric)),"log10(kappa)",ref)
        } else sim_lines_panel(p$data,v,p$metric,p$group,paste0(LETTERS[i],"  ",sim_metric_label(p$metric)),xlab,ref)
      }
    })
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note="All other parameters are fixed at the module A reference configuration. Errors concern log odds contrasts unless expressed as OR relative errors. SE ratio is the root mean squared model SE divided by the root mean squared empirical SD. Component time includes recorded learning, nuisance estimation and final computation; shared components are counted once in actual elapsed run time. Memory is a sampled process working-set peak and may retain R allocator memory; it describes measured process memory. Figures use the available resource measurements. Repetitions and effective settings are listed in the accompanying tables.")
  }
  for(view in intersect(c("kappa_lambda","L_M"),unique(views$view))) {
    z<-d[d$config_id%in%views$config_id[views$view==view],];if(!nrow(z))next
    xx<-if(view=="L_M")"L" else "kappa";yy<-if(view=="L_M")"M" else "lambda"
    id<-paste0("Figure_A_",view);title<-if(view=="L_M")"Learning sample size and kernel centers" else "Kernel-learning and interaction penalties"
    grids<-list()
    for(t in c("joint","individual"))for(m in c("RMSE","mean_OR_absolute_relative_error"))
      if(any(z$target==t & is.finite(z[[m]])))grids[[length(grids)+1L]]<-list(t=t,m=m)
    if(!length(grids))next
    sim_figure(file.path(folder,paste0(id,".png")),title,ceiling(length(grids)/2),min(2,length(grids)),cfg,function(){
      if(length(grids)==3)layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
      for(panel in grids) {
        t<-panel$t;m<-panel$m;a<-z[z$target==t,];X<-sort(unique(a[[xx]]));Y<-sort(unique(a[[yy]]))
        mat<-matrix(NA_real_,length(X),length(Y));mat[cbind(match(a[[xx]],X),match(a[[yy]],Y))]<-a[[m]]
        plot(c(.5,length(X)+.5),c(.5,length(Y)+.5),type="n",axes=FALSE,xlab=xx,ylab=yy,
          main=paste(sim_label(t),"\n",sim_metric_label(m)))
        pal<-hcl.colors(48,"YlOrRd",rev=TRUE);lim<-range(mat[is.finite(mat)])
        for(i in seq_along(X))for(j in seq_along(Y))if(is.finite(mat[i,j])) {
          ci<-if(diff(lim)>0)1L+floor(47*(mat[i,j]-lim[1])/diff(lim)) else 24L
          rect(i-.5,j-.5,i+.5,j+.5,col=pal[ci],border="white")
        }
        axis(1,seq_along(X),format(X,trim=TRUE));axis(2,seq_along(Y),format(Y,trim=TRUE),las=1)
        for(i in seq_along(X))for(j in seq_along(Y))if(is.finite(mat[i,j]))text(i,j,formatC(mat[i,j],digits=3,format="fg"),cex=.94)
      }
    })
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note="Cell values are aggregate errors over independent repetitions and common evaluation points. Lower values indicate smaller errors. Each panel has its own color scale; values, rather than colors across panels, should be compared. The complete configuration grid and success counts are supplied in the tables.")
  }
}

sim_plot_B <- function(result,folder,captions) {
  cfg<-result$config;d<-result$points;if(!nrow(d))return(invisible(NULL))
  methods<-c("KDR","SplineCLR","LinearCLR");cols<-c("#0072B2","#D55E00","#009E73")
  obtain<-function(target,g,m) {
    sg<-if(m=="KDR")g else if(g%in%c("G1","G2"))"correct" else "wrong"
    a<-d[d$target==target&d$method==m&d$scenario_gf==sg,,drop=FALSE]
    if(target=="joint"&&nrow(a)&&!any(abs(a$q-.5)<1e-12)) {
      ref<-a[1,,drop=FALSE];ref$q<-.5
      for(v in c("truth","mean_estimate","lower_repetition_quantile","upper_repetition_quantile"))ref[[v]]<-0
      a<-rbind(a,ref);a<-a[order(a$q),,drop=FALSE]
    }
    a
  }
  for(target in c("joint","individual")) {
    available<-Filter(function(g)any(vapply(methods,function(m)nrow(obtain(target,g,m))>0,logical(1))),paste0("G",1:4))
    if(!length(available))next
    id<-paste0("Figure_B_",target);title<-if(target=="joint")"Recovery of the joint exposure effect" else "Recovery of effects at observed exposure combinations"
    sim_figure(file.path(folder,paste0(id,".png")),title,ceiling(length(available)/2),min(2,length(available)),cfg,function(){
      if(length(available)==3)layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
      for(g in available) {
        sets<-lapply(methods,function(m)obtain(target,g,m));first<-Filter(function(x)nrow(x)>0,sets)[[1]]
        if(target=="joint") {
          yy<-c(first$truth,unlist(lapply(sets,function(a)c(a$lower_repetition_quantile,a$upper_repetition_quantile))))
          limits<-sim_range(yy);limits[2]<-limits[2]+.38*diff(limits)
          plot(range(first$q),limits,type="n",xlab="Common exposure quantile",ylab="Log odds ratio",main=sim_label(g))
          abline(h=0,lty=3,col="grey70")
          for(i in seq_along(sets)) {
            a<-sets[[i]];if(!nrow(a))next;a<-a[order(a$q),]
            polygon(c(a$q,rev(a$q)),c(a$lower_repetition_quantile,rev(a$upper_repetition_quantile)),
              col=adjustcolor(cols[i],alpha.f=.09),border=NA)
            lines(a$q,a$mean_estimate,col=cols[i],lwd=2.2)
          }
          lines(first$q,first$truth,col="black",lty=2,lwd=2.2)
        } else {
          yy<-unlist(lapply(sets,function(a)a$mean_estimate));lim<-sim_range(c(first$truth,yy))
          plot(lim,lim,type="n",xlab="True log odds ratio",ylab="Mean estimated log odds ratio",main=sim_label(g))
          abline(0,1,lty=2,lwd=2,col="black")
          ids<-unique(round(seq(1,nrow(first),length.out=min(cfg$report$scatter_points,nrow(first)))))
          for(i in seq_along(sets)) {
            a<-sets[[i]];if(!nrow(a))next;a<-a[ids,]
            points(a$truth,a$mean_estimate,pch=c(16,17,1)[i],col=adjustcolor(cols[i],alpha.f=.48),cex=.58)
          }
        }
        good<-which(vapply(sets,nrow,integer(1))>0)
        legend("topleft",c("Truth",sim_label(methods[good])),col=c("black",cols[good]),
          lty=c(2,rep(1,length(good))),lwd=2,bty="n",cex=.85)
      }
    })
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note=if(target=="joint")"All exposures follow their common marginal quantile, relative to the vector of marginal medians. Curves are means over usable repetitions. Shaded regions are empirical 2.5th–97.5th percentiles across repetitions, showing between-repetition variability. The reference point is zero by construction. CLR uses outcome adjustment, with each correct-adjustment fit supplying G1/G2 and each misspecified-adjustment fit supplying G3/G4. Failure counts and uncertainty metrics are reported separately." else
        "Each point represents one fixed independently sampled exposure combination, relative to the common reference vector. The vertical coordinate is its mean estimate across usable repetitions; the diagonal denotes the true contrast. A deterministic subset of evaluation points is displayed. Variability, RMSE and interval performance are reported in the tables, to quantify individual-fit stability alongside mean prediction agreement. CLR comparisons use the same adjustment settings as the joint-effect figure.")
  }
  id<-"Figure_B_OR_errors";title<-"Absolute relative errors on the odds-ratio scale"
  panels<-list()
  for(t in c("joint","individual"))for(g in paste0("G",1:4)) {
    series<-lapply(methods,function(m){sg<-if(m=="KDR")g else if(g%in%c("G1","G2"))"correct" else "wrong"
      key<-paste(unique(d$config_id)[1],m,sg,t,sep="|");v<-result$or_errors[[key]];v[is.finite(v)]})
    names(series)<-sim_label(methods);series<-Filter(function(x)length(x)&&any(is.finite(x)),series)
    if(length(series))panels[[length(panels)+1L]]<-list(t=t,g=g,series=series)
  }
  if(length(panels)) {
    sim_figure(file.path(folder,paste0(id,".png")),title,ceiling(length(panels)/2),2,cfg,function(){
      if(length(panels)%%2L==1L)layout(matrix(c(seq_len(length(panels)-1L),length(panels),length(panels)),ncol=2,byrow=TRUE))
      for(a in panels)boxplot(a$series,outline=FALSE,col=adjustcolor(cols[match(names(a$series),sim_label(methods))],alpha.f=.25),
        ylab="Absolute OR relative error (%)",xlab="Method",main=paste(a$g,sim_label(a$t)),cex.axis=.88)
    })
    captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note="Errors are 100 times abs(exp(estimated log OR minus true log OR) minus 1). Boxes summarize finite evaluation-point-by-repetition errors. Whiskers follow the 1.5-IQR rule; individual outliers are omitted from display only. Overflow errors remain infinite in numeric summaries and are never clipped to improve performance. Means, medians and upper quantiles are reported in the tables. Monte Carlo uncertainty uses independent repetitions.")
  }
}

sim_plot_C <- function(result,selection,folder,captions) {
  cfg<-result$config;a<-selection$overall;v<-selection$variables;if(!nrow(a))return(invisible(NULL))
  a<-a[is.finite(a$mean_TPR)|is.finite(a$mean_FPR),,drop=FALSE]
  if(!nrow(a))return(invisible(NULL))
  a<-a[order(a$null_count),]
  id<-"Figure_C_selection_rates";title<-"Variable identification as null exposures increase"
  panels<-c(if(any(is.finite(a$mean_TPR)))"TPR",if(any(is.finite(a$mean_FPR)))"FPR")
  sim_figure(file.path(folder,paste0(id,".png")),title,length(panels),1,cfg,function(){
    for(i in seq_along(panels)) {
      m<-panels[i];y<-100*a[[paste0("mean_",m)]];lo<-100*a[[paste0(m,"_lower")]];hi<-100*a[[paste0(m,"_upper")]]
      ok<-is.finite(y);col<-if(m=="TPR")"#0072B2" else "#D55E00"
      plot(a$null_count[ok],y[ok],type="b",pch=16,lwd=2,col=col,xlim=range(a$null_count),
        xlab="Number of null exposures",ylab=if(m=="TPR")"Mean true-positive rate (%)" else "Mean false-positive rate (%)",
        main=paste0(LETTERS[i],"  ",if(m=="TPR")"Active exposures identified" else "Null exposures selected"),
        ylim=if(m=="TPR")c(0,100) else sim_range(c(y,lo,hi),0,TRUE))
      jj<-ok&is.finite(lo)&is.finite(hi)&hi>lo
      if(any(jj))arrows(a$null_count[jj],lo[jj],a$null_count[jj],hi[jj],angle=90,code=3,length=.03,col=col)
    }
  })
  captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
    note="The true-positive rate is the fraction of eight active exposures selected in each repetition; the false-positive rate is the fraction of null exposures selected. Points are means across usable repetitions. Bars are 95% Monte Carlo t intervals computed from repetition-level rates, truncated to [0,1]. FPR is reported when null exposures are present. Variable selection uses numerical quadratic-form tail probabilities and Holm adjustment. Planned and usable repetitions are reported in the tables; unavailable tests retain an explicit availability status.")
  if(!nrow(v))return(invisible(NULL))
  v<-v[is.finite(v$rate),,drop=FALSE]
  active<-exposure_names[exposure_names%in%v$variable]
  if(!length(active))return(invisible(NULL))
  id<-"Figure_C_active_exposures";title<-"Identification rates of the active exposures"
  sim_figure(file.path(folder,paste0(id,".png")),title,ceiling(length(active)/2),min(2,length(active)),cfg,function(){
    sim_panel_layout(length(active))
    for(j in seq_along(active)) {
      z<-v[v$variable==active[j],];z<-z[order(z$null_count),]
      plot(z$null_count,100*z$rate,type="b",pch=16,lwd=2,col="#0072B2",ylim=c(0,100),
        xlab="Number of null exposures",ylab="Identification rate (%)",main=paste0(LETTERS[j],"  ",sim_label(active[j])))
      good<-is.finite(z$lower)&is.finite(z$upper)&z$upper>z$lower
      if(any(good))arrows(z$null_count[good],100*z$lower[good],z$null_count[good],100*z$upper[good],angle=90,code=3,length=.03,col="#0072B2")
    }
  })
  captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
    note="Each point is the proportion of usable independent repetitions selecting the indicated active exposure. Bars are 95% Wilson Monte Carlo intervals. Relative humidity participates through its interaction with temperature and has no main effect. Within a repetition, the active-exposure data are shared across candidate dimensions, permitting paired comparisons; uncertainty is computed across repetitions.")
}
