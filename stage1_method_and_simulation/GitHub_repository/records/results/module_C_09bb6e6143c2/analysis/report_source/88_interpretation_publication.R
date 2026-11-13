# Publication summaries of attribution, quantitative ORs and interaction tests.
sim_report_interpretation <- function(run,result,tables,figures,captions,module) {
  cfg<-result$config;x<-result$extended;raw<-result$interpretation
  emit<-function(d,name,title,note)sim_export_table(d,paste0("Table_",module,"_",name),title,note,tables,captions)
  display<-function(d)if(module=="B")sim_B_display_scenarios(d) else d
  note_contribution<-paste("Baseline Shapley values allocate h(x_active,b_other)-h(b), with b the marginal-median reference.",
    "The eight prespecified active exposures vary over the evaluation distribution; additional coordinates remain at their reference values.",
    "Importance is the mean absolute Shapley value; shares divide this quantity by its sum over the eight active exposures.",
    "Shares sum to 100% within a fit. Their sign-free magnitude includes main and interaction contributions.",
    "True and fitted allocations use identical combinations and reference levels. Uncertainty across repetitions is summarized by SD and MCSE.")
  emit(display(x$contributions),"contribution_shares","Recovery of exposure contribution shares",note_contribution)
  emit(display(x$OR),"quantitative_OR","Recovery of quantitative exposure-effect odds ratios",
    paste("All contrasts and confidence limits are exponentiated log OR contrasts.",
      "Single-exposure P75-versus-P25 contrasts hold the other exposures at their reference; temperature uses P10/P90 versus its specified optimum.",
      "Joint contrasts change all candidate exposures simultaneously. Mean CI limits average the lower and upper pointwise limits separately.",
      "RMS bias and RMSE use log OR; absolute OR error is 100*abs(exp(estimate-truth)-1)."))
  emit(display(x$pairs),"interaction_detection","Detection of overall pairwise interactions",
    paste("For each pair, the statistic is the mean squared four-point log-OR interaction contrast on the fixed evaluation grid.",
      "Its centered-Gaussian null distribution is a weighted chi-square sum computed numerically with covariance eigenvalues.",
      "Holm correction covers the four prespecified pairs within each fit; alpha is",cfg$interpretation$alpha,
      ". Wilson intervals quantify Monte Carlo uncertainty across independent repetitions. Usable tests are reported separately."))
  emit(display(x$directions),"interaction_direction","Recovery of significant interaction directions",
    paste("A correct direction requires a nonzero true contrast, a significant two-sided contrast test and agreement with its sign.",
      "The corrected version applies Holm across all non-anchor grid contrasts of all four pairs in each fit.",
      "Rates first average within each pair and repetition, then across repetitions; MCSE uses whole repetitions.",
      "The interaction sign can differ across the grid. Overall pair detection and pointwise direction detection are distinct criteria."))
  emit(x$scores,"interpretation_metrics","Configuration-level contribution and interaction performance",
    "Contribution share error is the repetition-average root mean squared percentage-point error across the eight prespecified active exposures, divided by 100. Additional coordinates remain at their reference values. Miss rates equal one minus corrected detection rates; direction errors use the pointwise Holm rule. Counts show available interpretation and interaction repetitions.")
  for(n in names(raw))if(nrow(raw[[n]]))kdr_csv(raw[[n]],file.path(tables,paste0(n,"_by_repetition.csv")))
  caption<-function(id,title,note)captions[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,note=note)
  if(module=="A") {
    z<-x$scores;v<-result$plan$design$membership
    if(nrow(z)) {
      z$contribution_error_pp<-100*z$contribution_share_error
      z$overall_detection<-100*(1-z$interaction_miss);z$direction_detection<-100*(1-z$direction_miss);z$wrong_direction<-100*z$direction_wrong
      measures<-c("contribution_error_pp","overall_detection","direction_detection","wrong_direction")
      labels<-c("Contribution share RMSE (pp)","Overall interaction detection (%)","Correct-direction detection (%)","Wrong-direction detection (%)")
      for(param in intersect(sim_parameter_names,unique(v$view))) {
        d<-z[z$config_id%in%v$config_id[v$view==param],,drop=FALSE];if(!nrow(d))next
        id<-paste0("Figure_A_interpretation_",param);title<-paste("Contribution and interaction performance across",param)
        sim_figure(file.path(figures,paste0(id,".png")),title,2,2,cfg,function(){
          for(k in seq_along(measures)) {
            xx<-d[[param]];if(param=="kappa")xx<-log10(xx);yy<-d[[measures[k]]];good<-is.finite(xx)&is.finite(yy)
            if(!any(good)){plot.new();title(main=labels[k]);text(.5,.5,"Insufficient usable results");next}
            ix<-which(good)[order(xx[good])]
            plot(xx[ix],yy[ix],type="b",pch=16,col="#0072B2",lwd=2,
              xlab=if(param=="kappa")"log10(kappa)" else param,ylab=labels[k],main=paste(LETTERS[k],labels[k]))
          }
        });caption(id,title,"Points summarize independent repetitions under G1. Detection uses Holm correction; contribution error uses all eight exposure shares. Underlying values and usable counts are tabulated.")
      }
      for(view in intersect(c("kappa_lambda","L_M"),unique(v$view))) {
        d<-z[z$config_id%in%v$config_id[v$view==view],,drop=FALSE];if(!nrow(d))next
        xx<-if(view=="L_M")"L" else "kappa";yy<-if(view=="L_M")"M" else "lambda"
        X<-sort(unique(d[[xx]]));Y<-sort(unique(d[[yy]]));id<-paste0("Figure_A_interpretation_",view);title<-"Contribution and interaction parameter combinations"
        sim_figure(file.path(figures,paste0(id,".png")),title,2,2,cfg,function(){
          for(k in seq_along(measures)) {
            vals<-d[[measures[k]]];lim<-sim_range(vals);pal<-hcl.colors(48,"YlOrRd",rev=TRUE)
            plot(c(.5,length(X)+.5),c(.5,length(Y)+.5),type="n",axes=FALSE,xlab=xx,ylab=yy,main=labels[k])
            axis(1,seq_along(X),format(X,trim=TRUE));axis(2,seq_along(Y),format(Y,trim=TRUE),las=1)
            for(i in which(is.finite(vals))) {
              a<-match(d[[xx]][i],X);b<-match(d[[yy]][i],Y);ci<-max(1,min(48,1+floor(47*(vals[i]-lim[1])/diff(lim))))
              rect(a-.5,b-.5,a+.5,b+.5,col=pal[ci],border="white");text(a,b,formatC(vals[i],digits=3,format="fg"),cex=.95)
            }
          }
        });caption(id,title,"Cell values report contribution-share error and corrected interaction detection. Each panel uses its own color scale; numerical values define the comparisons.")
      }
    }
    return(invisible(x))
  }
  plot_nulls<-if(is.null(cfg$report$C_null_counts))sort(unique(result$plan$tasks$null_count)) else cfg$report$C_null_counts
  if(module=="C")for(n in names(x))if(nrow(x[[n]])&&"null_count"%in%names(x[[n]]))
    x[[n]]<-x[[n]][x[[n]]$null_count%in%plot_nulls,,drop=FALSE]
  # Contribution bars: one truth bar followed by the fitted methods for each setting.
  shares<-display(x$contributions)
  settings<-if(module=="B")unique(shares$scenario_gf) else sort(unique(shares$null_count))
  if(nrow(shares)&&length(settings)) {
    id<-paste0("Figure_",module,"_contribution_shares");title<-"Exposure contribution shares"
    palette<-c("#0072B2","#56B4E9","#D55E00","#E69F00","#009E73","#CC79A7","#882255","#44AA99")
    sim_figure(file.path(figures,paste0(id,".png")),title,ceiling(length(settings)/2),min(2,length(settings)),cfg,function(){
      for(level in settings) {
        z<-if(module=="B")shares[shares$scenario_gf==level,] else shares[shares$null_count==level,]
        methods<-c("KDR","SplineCLR","LinearCLR");methods<-methods[methods%in%z$method]
        mat<-matrix(0,8,length(methods)+1L);rownames(mat)<-exposure_names
        put<-function(a,col){for(j in 1:8)mat[j,col]<<-sum(a$value[a$variable==exposure_names[j]])}
        a<-z[z$method==methods[1],];a$value<-a$true_share;put(a,1)
        for(k in seq_along(methods)){a<-z[z$method==methods[k],];a$value<-a$estimated_share;put(a,k+1)}
        keep<-rowSums(mat)>0;barplot(mat[keep,,drop=FALSE],col=palette[keep],border=NA,
          names.arg=c("Truth",sim_label(methods)),ylab="Contribution share (%)",ylim=c(0,150),yaxt="n",
          main=if(module=="B")level else paste(level,"null exposures"),cex.names=.88)
        axis(2,seq(0,100,20))
        legend("top",legend=sim_label(exposure_names)[keep],fill=palette[keep],bty="n",ncol=3,cex=.74)
      }
    });caption(id,title,note_contribution)
  }
  # Quantitative OR forests provide selected contrasts without duplicating curves.
  or<-display(x$OR)
  levels<-if(module=="B")unique(or$scenario_gf) else sort(unique(or$null_count))
  for(level in levels) {
    d<-if(module=="B")or[or$scenario_gf==level,] else or[or$null_count==level,]
    if(!nrow(d))next
    id<-paste0("Figure_",module,"_quantitative_OR_",if(module=="B")level else paste0("q",level))
    title<-paste("Quantitative odds-ratio recovery:",if(module=="B")level else paste(level,"null exposures"))
    cf<-cfg;cf$report$width_inches<-max(16,cf$report$width_inches)
    sim_figure(file.path(figures,paste0(id,".png")),title,1,2,cf,function(){
      layout(matrix(c(1,1,2,3),2,2,byrow=TRUE),heights=c(.48,4.3))
      par(cex=cf$report$base_cex,mar=c(0,0,0,0));plot.new()
      all_methods<-c("KDR","SplineCLR","LinearCLR");all_methods<-all_methods[all_methods%in%d$method]
      colors<-c(KDR="#0072B2",SplineCLR="#D55E00",LinearCLR="#009E73")
      legend("center",c("Truth",sim_label(all_methods)),col=c("black",colors[all_methods]),
        pch=c(4,rep(16,length(all_methods))),horiz=TRUE,bty="n",cex=.9)
      for(target in c("joint","single")) {
        z<-d[d$target==target,];group<-paste(z$exposure,z$contrast,sep=" | ");labs<-unique(group);J<-length(labs)
        vals<-c(z$true_OR,z$mean_OR,z$mean_CI_lower,z$mean_CI_upper);lim<-range(vals[is.finite(vals)&vals>0])
        par(mar=c(5.8,11,3.2,1.2));plot(lim,c(.4,J+.6),log="x",type="n",yaxt="n",xlab="Odds ratio (log scale)",ylab="",main=sim_label(target))
        labels<-vapply(labs,function(l){a<-z[match(l,group),];paste(if(a$exposure=="all")"Joint" else sim_label(a$exposure),a$contrast)},character(1))
        axis(2,seq_len(J),labels,las=1,cex.axis=.82);abline(v=1,lty=3,col="grey70")
        methods<-c("KDR","SplineCLR","LinearCLR");methods<-methods[methods%in%z$method];cols<-c(KDR="#0072B2",SplineCLR="#D55E00",LinearCLR="#009E73")
        for(k in seq_along(methods)) {
          a<-z[z$method==methods[k],];y<-match(paste(a$exposure,a$contrast,sep=" | "),labs)+(k-(length(methods)+1)/2)*.15
          good<-is.finite(a$mean_CI_lower)&is.finite(a$mean_CI_upper)&a$mean_CI_lower>0
          segments(a$mean_CI_lower[good],y[good],a$mean_CI_upper[good],y[good],col=adjustcolor(cols[methods[k]],.65))
          points(a$mean_OR,y,col=cols[methods[k]],pch=16,cex=.7)
        }
        a<-z[!duplicated(group),];points(a$true_OR,seq_len(J),pch=4,lwd=1.7)
      }
    });caption(id,title,"Crosses denote true ORs; colored points denote repetition-mean ORs. Horizontal segments connect the separately averaged pointwise 95% interval endpoints. Full contrast definitions, errors and interval performance are tabulated.")
  }
  # Two detection objectives, with one row per specified interaction.
  pair<-display(x$pairs);direction<-display(x$directions)
  if(nrow(pair)&&nrow(direction)) {
    pairs<-unique(pair$component);id<-paste0("Figure_",module,"_interaction_detection");title<-"Interaction significance and direction recovery"
    sim_figure(file.path(figures,paste0(id,".png")),title,length(pairs),2,cfg,function(){
      for(component in pairs)for(mode in c("overall","direction")) {
        z<-if(mode=="overall")pair[pair$component==component,] else direction[direction$component==component,]
        measure<-if(mode=="overall")"corrected_detection_percent" else "corrected_correct_percent"
        levels<-if(module=="B")paste0("G",1:4) else sort(unique(z$null_count))
        group<-if(module=="B")z$scenario_gf else z$null_count
        methods<-c("KDR","SplineCLR","LinearCLR");methods<-methods[methods%in%z$method];cols<-c(KDR="#0072B2",SplineCLR="#D55E00",LinearCLR="#009E73")
        plot(c(.7,length(levels)+.3),c(0,110),type="n",xaxt="n",xlab=if(module=="B")"Nuisance scenario" else "Null exposures",
          ylab=if(mode=="overall")"Overall detection (%)" else "Correct-direction detection (%)",
          main=paste(paste(sim_label(strsplit(component,"__",fixed=TRUE)[[1]]),collapse=" × "),if(mode=="overall")"| significance" else "| direction"))
        axis(1,seq_along(levels),levels);abline(h=c(0,50,100),lty=3,col="grey90")
        for(k in seq_along(methods)) {
          take<-which(z$method==methods[k]);ix<-match(group[take],levels);o<-order(ix);take<-take[o];ix<-ix[o]
          lines(ix,z[[measure]][take],col=cols[methods[k]],lwd=2,type="b",pch=16)
        }
        legend("bottomleft",sim_label(methods),col=cols[methods],lty=1,pch=16,bty="n",cex=.75)
      }
    });caption(id,title,"Left: Holm-corrected pairwise interaction detection. Right: correct significant directions after Holm correction of the full interaction grid. Rates are calculated per repetition; tables give uncertainty and wrong-direction rates. Additive CLR models have exactly zero four-point interaction contrasts.")
  }
  if(module=="C")for(q in intersect(sort(unique(result$plan$tasks$null_count)),plot_nulls)) {
    pp<-result$plan;pp$tasks<-pp$tasks[pp$tasks$null_count==q,,drop=FALSE]
    cf<-cfg;cf$report$response_prefix<-paste0("C_q",q)
    out<-sim_B_response_collect(run,cf,pp,sim_B_response_settings(cf$report,cf))
    for(name in c("metrics","status"))if(nrow(out[[name]]))sim_export_table(out[[name]],
      paste0("Table_C_q",q,"_response_",name),paste("Fixed-background responses:",q,"null exposures"),
      "The eight active-exposure curves and four interactions are evaluated with all other candidate exposures at their reference. Values use the G1 KDR-CC fit at this dimension; availability and independent-repetition errors are reported.",tables,captions)
    if(nrow(out$points))kdr_csv(out$points,file.path(tables,paste0("C_q",q,"_response_points.csv")))
    sim_plot_B_responses(out,cf,figures,captions)
  }
  invisible(x)
}
