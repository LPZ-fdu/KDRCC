# Publication graphics drawn from matched-cohort summaries.
rp_open_plot<-function(path,width,height,title,settings) {
  font<-"serif"
  if(.Platform$OS.type=="windows") {
    grDevices::windowsFonts(KDRCCSerif=grDevices::windowsFont(settings$report$font));font<-"KDRCCSerif"
  }
  grDevices::png(path,width=width,height=height,units="in",res=settings$report$dpi,pointsize=settings$report$base_size)
  par(family=font,ps=settings$report$base_size,mar=c(5,5.1,3.1,1.1),oma=c(0,0,3,0),las=1,
    cex.axis=.95,cex.lab=1.05,cex.main=1.05)
  invisible(title)
}
rp_plot_title<-function(title)mtext(title,side=3,outer=TRUE,line=.8,cex=1.2,font=2)
rp_range<-function(x) {
  x<-x[is.finite(x)];if(!length(x))stop("A plotted quantity has no finite values.")
  r<-range(x);d<-diff(r);if(d==0)d<-max(.01,abs(r[1])*.05)
  r+c(-1,1)*.06*d
}
rp_legend<-function(truth=TRUE) {
  par(mar=c(0,0,0,0));plot.new()
  labs<-c(if(truth)"Generating truth","Archived KDR-CC G1","R package KDR-CC G1")
  cols<-c(if(truth)"black","#D55E00","#0072B2");lty<-c(if(truth)2,1,3)
  legend("center",labs,col=cols,lty=lty,lwd=2,horiz=TRUE,bty="n",cex=.9)
}
rp_curve<-function(x,truth,a,b,xlab,main="") {
  plot(x,truth,type="n",ylim=rp_range(c(truth,a$mean_estimate,b$mean_estimate,b$mean_CI_lower,b$mean_CI_upper,
    b$lower_repetition_quantile,b$upper_repetition_quantile)),xlab=xlab,ylab="Log odds ratio",main=main)
  band<-function(lo,hi,alpha) {
    i<-which(is.finite(x)&is.finite(lo)&is.finite(hi))
    if(length(i)>1)polygon(c(x[i],rev(x[i])),c(lo[i],rev(hi[i])),col=adjustcolor("#0072B2",alpha.f=alpha),border=NA)
  }
  band(b$mean_CI_lower,b$mean_CI_upper,.11);band(b$lower_repetition_quantile,b$upper_repetition_quantile,.22)
  abline(h=0,lty=3,col="grey70")
  lines(x,a$mean_estimate,col="#D55E00",lwd=2.2)
  lines(x,b$mean_estimate,col="#0072B2",lwd=2.2,lty=3)
  lines(x,truth,col="black",lwd=2,lty=2)
}
rp_figures<-function(s,analysis,settings,e) {
  folder<-file.path(analysis,"figures");dir.create(folder,showWarnings=FALSE)
  r<-s$reference;p<-s$package;vars<-e$exposure_names
  # Joint curve.
  rp_open_plot(file.path(folder,"Figure_G1_joint.png"),11,7,"",settings)
  layout(matrix(c(1,2),2,1),heights=c(.65,5.3));rp_legend();par(mar=c(5,5.1,2.1,1.1))
  a<-r$points[r$points$target=="joint",];b<-p$points[p$points$target=="joint",]
  o<-order(a$q);a<-a[o,];b<-b[match(a$index,b$index),]
  rp_curve(a$q,a$truth,a,b,"Common exposure quantile")
  rp_plot_title("Recovery of the joint exposure effect: G1");dev.off()
  # Observed combinations.
  rp_open_plot(file.path(folder,"Figure_G1_individual.png"),10,8,"",settings)
  a<-r$points[r$points$target=="individual",];b<-p$points[p$points$target=="individual",];b<-b[match(a$index,b$index),]
  lim<-rp_range(c(a$truth,a$mean_estimate,b$mean_estimate))
  plot(a$truth,a$mean_estimate,type="n",xlim=lim,ylim=lim,xlab="True log odds ratio",ylab="Mean estimated log odds ratio",asp=1)
  abline(0,1,lty=2,lwd=1.7)
  points(a$truth,a$mean_estimate,pch=1,col=adjustcolor("#D55E00",.65),cex=.7)
  points(a$truth,b$mean_estimate,pch=16,col=adjustcolor("#0072B2",.55),cex=.45)
  legend("topleft",c("Identity line","Archived KDR-CC G1","R package KDR-CC G1"),
    lty=c(2,NA,NA),pch=c(NA,1,16),col=c("black","#D55E00","#0072B2"),bty="n",cex=.9)
  rp_plot_title("Recovery at observed exposure combinations: G1");dev.off()
  # Eight single-coordinate curves.
  rp_open_plot(file.path(folder,"Figure_G1_single.png"),18,10,"",settings)
  layout(rbind(rep(1,4),matrix(2:9,2,4,byrow=TRUE)),heights=c(.6,4.2,4.2));rp_legend()
  for(v in vars) {
    a<-r$response$points[r$response$points$target=="single"&r$response$points$component==v,]
    b<-p$response$points[p$response$points$target=="single"&p$response$points$component==v,]
    a<-a[order(a$concentration_1),];b<-b[match(a$index,b$index),]
    par(mar=c(5.1,4.8,2.8,.9))
    rp_curve(a$concentration_1,a$truth,a,b,e$sim_B_exposure_axis(v),e$sim_label(v))
  }
  rp_plot_title("Single-exposure effect recovery: G1");dev.off()
  # Interaction contrasts.
  rp_open_plot(file.path(folder,"Figure_G1_interactions.png"),13,11,"",settings)
  layout(rbind(c(1,1),matrix(2:5,2,2,byrow=TRUE)),heights=c(.55,4.2,4.2));rp_legend(FALSE)
  pairs<-unique(r$response$points$component[r$response$points$target=="interaction"])
  for(pair in pairs) {
    a<-r$response$points[r$response$points$component==pair&r$response$points$target=="interaction"&r$response$points$eligible,]
    b<-p$response$points[p$response$points$component==pair&p$response$points$target=="interaction"&p$response$points$eligible,]
    b<-b[match(a$index,b$index),];lim<-rp_range(c(a$truth,a$mean_estimate,b$mean_estimate,b$lower_repetition_quantile,b$upper_repetition_quantile))
    title<-paste(e$sim_label(a$exposure_1[1]),e$sim_label(a$exposure_2[1]),sep=" × ")
    if(nchar(title)>30)title<-paste(e$sim_label(a$exposure_1[1]),e$sim_label(a$exposure_2[1]),sep=" ×\n")
    par(mar=c(4.7,4.9,3.7,.9))
    plot(lim,lim,type="n",xlab="True interaction (log OR)",ylab="Mean estimated interaction (log OR)",main=title,asp=1)
    abline(0,1,lty=2,lwd=1.5);abline(h=0,v=0,lty=3,col="grey80")
    segments(b$truth,b$lower_repetition_quantile,b$truth,b$upper_repetition_quantile,col=adjustcolor("#0072B2",.12))
    points(a$truth,a$mean_estimate,pch=1,col="#D55E00",cex=.7)
    points(b$truth,b$mean_estimate,pch=16,col=adjustcolor("#0072B2",.65),cex=.5)
  }
  rp_plot_title("Interaction-effect recovery: G1");dev.off()
  # Contribution shares, three stacked bars.
  rp_open_plot(file.path(folder,"Figure_G1_contribution_shares.png"),12,8,"",settings)
  palette<-c("#0072B2","#56B4E9","#D55E00","#E69F00","#009E73","#CC79A7","#882255","#44AA99")
  layout(matrix(c(1,2),2,1),heights=c(1.1,4.8));par(mar=c(0,0,0,0));plot.new()
  legend("center",e$sim_label(vars),fill=palette,ncol=4,bty="n",cex=.9)
  a<-r$extended$contributions[match(vars,r$extended$contributions$variable),]
  b<-p$extended$contributions[match(vars,p$extended$contributions$variable),]
  par(mar=c(5,5.2,1.6,1));barplot(cbind(a$true_share,a$estimated_share,b$estimated_share),col=palette,border=NA,
    names.arg=c("Generating truth","Archived KDR-CC G1","R package KDR-CC G1"),ylim=c(0,100),ylab="Contribution share (%)",cex.names=.9)
  rp_plot_title("Exposure contribution shares: G1");dev.off()
  # Overall detection and pointwise direction rates.
  rp_open_plot(file.path(folder,"Figure_G1_interaction_detection.png"),13,7,"",settings)
  layout(rbind(c(1,1),c(2,3)),heights=c(.6,5));rp_legend(FALSE)
  labels<-vapply(strsplit(pairs,"__",fixed=TRUE),function(z)paste(e$sim_label(z),collapse=" ×\n"),character(1))
  for(mode in c("overall","direction")) {
    a<-if(mode=="overall")r$extended$pairs else r$extended$directions
    b<-if(mode=="overall")p$extended$pairs else p$extended$directions
    a<-a[match(pairs,a$component),];b<-b[match(pairs,b$component),]
    field<-if(mode=="overall")"corrected_detection_percent" else "corrected_correct_percent"
    par(mar=c(7.3,4.7,2.5,1));barplot(rbind(a[[field]],b[[field]]),beside=TRUE,names.arg=labels,
      col=c("#D55E00","#0072B2"),border=NA,ylim=c(0,105),ylab="Detection (%)",cex.names=.75,
      main=if(mode=="overall")"Overall interaction detection" else "Correct significant directions")
    abline(h=c(50,100),lty=3,col="grey80")
  }
  rp_plot_title("Interaction significance and direction recovery: G1");dev.off()
  # Quantitative odds-ratio contrasts.
  rp_open_plot(file.path(folder,"Figure_G1_quantitative_OR.png"),15,10,"",settings)
  layout(rbind(c(1,1),c(2,3)),heights=c(.55,6.5));rp_legend()
  for(target in c("joint","single")) {
    a<-r$extended$OR[r$extended$OR$target==target,];b<-p$extended$OR[p$extended$OR$target==target,]
    key<-paste(a$exposure,a$contrast,sep="|");b<-b[match(key,paste(b$exposure,b$contrast,sep="|")),]
    labels<-paste(ifelse(a$exposure=="all","Joint",e$sim_label(a$exposure)),a$contrast,sep=" | ")
    vals<-c(a$true_OR,a$mean_OR,b$mean_OR,b$mean_CI_lower,b$mean_CI_upper)
    vals<-vals[is.finite(vals)&vals>0];lim<-range(vals);lim<-exp(rp_range(log(lim)))
    par(mar=c(5.2,13,2.7,1));J<-nrow(a)
    plot(lim,c(.5,J+.5),type="n",log="x",yaxt="n",xlab="Odds ratio (log scale)",ylab="",main=if(target=="joint")"Joint exposure contrasts" else "Single-exposure contrasts")
    axis(2,seq_len(J),labels,las=1,cex.axis=.72);abline(v=1,lty=3,col="grey70")
    segments(a$mean_CI_lower,seq_len(J)-.09,a$mean_CI_upper,seq_len(J)-.09,col="#D55E00",lwd=1.5)
    segments(b$mean_CI_lower,seq_len(J)+.09,b$mean_CI_upper,seq_len(J)+.09,col="#0072B2",lwd=1.5)
    points(a$mean_OR,seq_len(J)-.09,col="#D55E00",pch=1,cex=.8);points(b$mean_OR,seq_len(J)+.09,col="#0072B2",pch=16,cex=.7)
    points(a$true_OR,seq_len(J),pch=4,lwd=1.8)
  }
  rp_plot_title("Quantitative odds-ratio recovery: G1");dev.off()
  # OR absolute relative error distributions.
  rp_open_plot(file.path(folder,"Figure_G1_OR_errors.png"),12,7,"",settings)
  par(mfrow=c(1,2),mar=c(5.5,5,3,1))
  for(target in c("joint","individual")) {
    z<-list(r$effects$or_errors[[target]],p$effects$or_errors[[target]])
    boxplot(z,names=c("Archived G1","R package"),col=c("#EFC5A8","#A8CFE4"),outline=FALSE,
      ylab="Absolute relative OR error (%)",main=if(target=="joint")"Joint exposure effect" else "Observed-combination effect")
  }
  rp_plot_title("Distribution of absolute relative OR errors: G1");dev.off()
  # Identification of the eight active exposures.
  rp_open_plot(file.path(folder,"Figure_G1_variable_identification.png"),12,7,"",settings)
  a<-r$selection$variables[match(vars,r$selection$variables$variable),];b<-p$selection$variables[match(vars,p$selection$variables$variable),]
  par(mar=c(6,5.1,3,1));x<-barplot(100*rbind(a$rate,b$rate),beside=TRUE,col=c("#D55E00","#0072B2"),
    border=NA,names.arg=e$sim_label(vars),cex.names=.8,ylim=c(0,112),ylab="Identification rate (%)")
  arrows(x[1,],100*a$lower,x[1,],100*a$upper,angle=90,code=3,length=.035,col="#D55E00")
  arrows(x[2,],100*b$lower,x[2,],100*b$upper,angle=90,code=3,length=.035,col="#0072B2")
  legend("top",c("Archived KDR-CC G1","R package KDR-CC G1"),fill=c("#D55E00","#0072B2"),horiz=TRUE,bty="n",cex=.85)
  rp_plot_title("Identification of active exposures: G1");dev.off()
  # Numerical diagnostics have their own output directory.
  d<-s$diagnostics;diagnostics<-file.path(analysis,"model_diagnostics")
  fields<-c("iterations","maximum_gradient","objective","curvature_condition")
  da<-d[d$implementation=="Archived G1",];db<-d[d$implementation=="R package",]
  db<-db[match(da$replicate,db$replicate),]
  fields<-fields[vapply(fields,function(n) {
    good<-is.finite(da[[n]])&is.finite(db[[n]])
    if(n%in%c("maximum_gradient","curvature_condition"))good<-good&da[[n]]>0&db[[n]]>0
    any(good)
  },logical(1))]
  if(length(fields)) {
    rp_open_plot(file.path(diagnostics,"Figure_G1_model_diagnostics.png"),12,9,"",settings)
    if(length(fields)==3L)layout(matrix(c(1,2,3,3),2,2,byrow=TRUE)) else
      par(mfrow=c(ceiling(length(fields)/2),min(2,length(fields))))
    par(mar=c(4.5,5,3,1))
    for(field in fields) {
      a<-d[d$implementation=="Archived G1",];b<-d[d$implementation=="R package",]
      b<-b[match(a$replicate,b$replicate),];good<-is.finite(a[[field]])&is.finite(b[[field]])
      islog<-field%in%c("maximum_gradient","curvature_condition")
      if(islog)good<-good&a[[field]]>0&b[[field]]>0
      if(!any(good))stop("Diagnostic data changed while plotting: ",field)
      lim<-rp_range(c(a[[field]][good],b[[field]][good]));if(islog)lim<-exp(rp_range(log(c(a[[field]][good],b[[field]][good]))))
      plot(a[[field]][good],b[[field]][good],pch=16,col="#0072B2",log=if(islog)"xy" else "",xlim=lim,ylim=lim,
        xlab="Archived G1",ylab="R package",main=gsub("_"," ",field));abline(0,1,lty=2)
    }
    rp_plot_title("Numerical estimation diagnostics: G1");dev.off()
  }
}
rp_write_captions<-function(analysis,R) {
  notes<-c(
    "Figure_G1_joint.png|Recovery of the joint exposure effect|All eight exposures follow the archived joint quantile path. Curves show generating truth and means over the same matched independent repetitions. The blue light band averages pointwise 95% confidence limits separately; the darker band is the empirical central 95% range of package point estimates.",
    "Figure_G1_individual.png|Recovery at observed exposure combinations|Each point represents one archived evaluation combination. Orange open circles and blue points give the archived and package repetition-mean log ORs; the diagonal is equality with generating truth.",
    "Figure_G1_single.png|Single-exposure effect recovery|All other exposures remain at the archived marginal-median reference. Concentration values, references, contrasts and repetition IDs are identical in the two implementations. Blue bands have the same definitions as the joint-curve figure.",
    "Figure_G1_interactions.png|Interaction-effect recovery|Four-point log-OR contrasts measure departures from additivity. Points compare truth with matched-cohort mean predictions. Blue vertical segments show central 95% repetition ranges. Anchor contrasts are excluded.",
    "Figure_G1_contribution_shares.png|Exposure contribution shares|Baseline Shapley allocations use the same reference and evaluation combinations. Mean absolute allocations are normalized within the eight active exposures for each fit; bar heights are generating shares and averages of fitted shares.",
    "Figure_G1_interaction_detection.png|Interaction significance and direction recovery|Overall tests use Holm correction over the four prespecified pairs. Correct significant directions use Holm over all eligible grid contrasts and require the true sign. Grid rates are averaged within each repetition before averaging across repetitions.",
    "Figure_G1_quantitative_OR.png|Quantitative odds-ratio recovery|Crosses show true ORs. Orange and blue points are averages of OR estimates. Horizontal segments connect separately averaged pointwise interval endpoints. Joint quantile contrasts, single-exposure interquartile contrasts and temperature-versus-optimum contrasts follow the archived definitions.",
    "Figure_G1_OR_errors.png|Absolute relative OR error distributions|Errors equal 100*abs(exp(estimated log OR minus true log OR)-1), pooled over matched repetitions and eligible evaluation points. Boxplots show quartiles and conventional whiskers; all numerical values enter the summary tables.",
    "Figure_G1_variable_identification.png|Identification of active exposures|Derivative-based participation tests use Holm correction. Bars give the proportion of repetitions selecting each active exposure; error bars give Wilson Monte Carlo intervals.")
  output<-c("# Figure and table captions","",paste("Matched independent repetitions:",R),"")
  for(x in notes){p<-strsplit(x,"|",fixed=TRUE)[[1]];output<-c(output,paste0("**",p[1],". ",p[2],".** ",p[3]),"")}
  output<-c(output,"**Model diagnostics.** The dedicated model_diagnostics directory compares final iteration counts, maximum gradients, regularized objective values and exact-curvature condition numbers.","",
    "**Tables.** CSV files preserve numerical columns for further analysis. Matching HTML files include a title, header and note and use three-line table rules. Availability, matched cohort, pointwise effects, errors, interval performance, contributions, quantitative ORs, participation and interaction detection are reported separately.")
  writeLines(output,file.path(analysis,"Figure_and_table_captions.md"))
}
