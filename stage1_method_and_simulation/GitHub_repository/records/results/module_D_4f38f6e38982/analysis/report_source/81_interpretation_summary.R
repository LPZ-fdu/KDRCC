# Repetition-level summaries preserve the simulation unit.
sim_interpretation_collect <- function(records,plan) {
  keys<-c("config_id","method","scenario_gf","null_count","p","replicate",sim_parameter_names)
  out<-setNames(rep(list(list()),5),c("contributions","OR","pairs","directions","scores"))
  for(i in seq_along(records)) {
    a<-records[[i]];if(!isTRUE(a$interpretation_ok))next
    e<-a$interpretation;j<-plan[i,keys,drop=FALSE]
    attach<-function(z)if(!is.null(z)&&nrow(z))cbind(j[rep(1,nrow(z)),,drop=FALSE],z) else NULL
    out$contributions[[i]]<-attach(e$contributions);out$OR[[i]]<-attach(e$OR)
    pairs<-data.frame();directions<-data.frame()
    if(isTRUE(e$detection$ok)) {
      pairs<-e$detection$value$pairs;pts<-e$detection$value$points
      directions<-kdr_bind(lapply(sim_groups(pts,"component"),function(ix){z<-pts[ix,]
        data.frame(component=z$component[1],nonzero_truth_points=sum(z$true_signal),
          valid_direction_points=sum(z$true_signal & z$point_ok),
          nominal_correct=sim_mean(z$nominal_correct_direction),nominal_wrong=sim_mean(z$nominal_wrong_direction),
          corrected_correct=sim_mean(z$holm_grid_correct_direction),corrected_wrong=sim_mean(z$holm_grid_wrong_direction))}))
      out$pairs[[i]]<-attach(pairs);out$directions[[i]]<-attach(directions)
    }
    share_error<-e$contributions$estimated_share-e$contributions$true_share
    valid_interaction<-nrow(pairs)>0&&all(pairs$test_ok)&&nrow(directions)>0&&
      all(directions$valid_direction_points==directions$nonzero_truth_points)
    out$scores[[i]]<-attach(data.frame(contribution_share_error=sqrt(mean(share_error^2))/100,
      null_contribution_fraction=sum(e$contributions$estimated_share[!e$contributions$true_active])/100,
      interaction_miss=if(valid_interaction)1-mean(pairs$holm_pairs_detected[pairs$true_interaction]) else NA_real_,
      direction_miss=if(valid_interaction)1-mean(directions$corrected_correct) else NA_real_,
      direction_wrong=if(valid_interaction)mean(directions$corrected_wrong) else NA_real_,
      contribution_identity_error=e$contribution_identity_error,interaction_available=valid_interaction))
  }
  lapply(out,kdr_bind)
}
sim_extended_summary <- function(x) {
  cols<-c("config_id","method","scenario_gf","null_count","p")
  result<-list()
  result$contributions<-kdr_bind(lapply(sim_groups(x$contributions,c(cols,"variable","true_active")),function(ix){
    z<-x$contributions[ix,];good<-is.finite(z$estimated_share);z<-z[good,,drop=FALSE];if(!nrow(z))return(NULL)
    err<-z$estimated_share-z$true_share
    cbind(z[1,c(cols,"variable","true_active"),drop=FALSE],data.frame(repetitions=nrow(z),
      true_share=mean(z$true_share),estimated_share=mean(z$estimated_share),
      share_bias=mean(err),share_RMSE=sqrt(mean(err^2)),share_SD=sim_sd(z$estimated_share),share_MCSE=sim_mcse(z$estimated_share)))
  }))
  result$pairs<-kdr_bind(lapply(sim_groups(x$pairs,c(cols,"component")),function(ix){
    z<-x$pairs[ix,];valid<-z$test_ok;v<-as.numeric(z$holm_pairs_detected[valid]);n<-length(v);rate<-sim_mean(v)
    crit<-qnorm(.975);den<-1+crit^2/n;center<-(rate+crit^2/(2*n))/den
    half<-crit*sqrt(rate*(1-rate)/n+crit^2/(4*n^2))/den
    cbind(z[1,c(cols,"component","exposure_1","exposure_2"),drop=FALSE],data.frame(
      records=nrow(z),usable_tests=n,raw_detection_percent=100*sim_mean(as.numeric(z$raw_detected[valid])),
      corrected_detection_percent=100*rate,lower_MC=100*max(0,center-half),upper_MC=100*min(1,center+half)))
  }))
  result$directions<-kdr_bind(lapply(sim_groups(x$directions,c(cols,"component")),function(ix){
    z<-x$directions[ix,];v<-z$corrected_correct
    cbind(z[1,c(cols,"component"),drop=FALSE],data.frame(repetitions=nrow(z),
      nominal_correct_percent=100*sim_mean(z$nominal_correct),nominal_wrong_percent=100*sim_mean(z$nominal_wrong),
      corrected_correct_percent=100*sim_mean(v),corrected_wrong_percent=100*sim_mean(z$corrected_wrong),
      MCSE_correct_percent=100*sim_mcse(v)))
  }))
  result$OR<-kdr_bind(lapply(sim_groups(x$OR,c(cols,"target","exposure","contrast")),function(ix){
    z<-x$OR[ix,];err<-z$estimate-z$truth
    cbind(z[1,c(cols,"target","exposure","contrast"),drop=FALSE],data.frame(repetitions=nrow(z),
      true_OR=z$true_OR[1],mean_OR=mean(z$OR),median_OR=median(z$OR),
      mean_CI_lower=sim_mean(z$OR_lower),mean_CI_upper=sim_mean(z$OR_upper),
      RMS_Bias=abs(mean(err)),RMSE=sqrt(mean(err^2)),empirical_SD=sim_sd(z$estimate),
      model_SE_RMS=sqrt(sim_mean(z$SE^2)),coverage_percent=100*sim_mean(z$covered),
      correct_significant_percent=100*sim_mean(z$correct_significant),wrong_significant_percent=100*sim_mean(z$wrong_significant),
      mean_absolute_OR_error=mean(z$OR_absolute_relative_error),median_absolute_OR_error=median(z$OR_absolute_relative_error)))
  }))
  result$scores<-kdr_bind(lapply(sim_groups(x$scores,cols),function(ix){
    z<-x$scores[ix,];names<-c("contribution_share_error","null_contribution_fraction","interaction_miss","direction_miss","direction_wrong")
    cbind(z[1,c(cols,sim_parameter_names),drop=FALSE],data.frame(interpretation_repetitions=nrow(z),
      interaction_repetitions=sum(z$interaction_available)),as.data.frame(lapply(z[names],sim_mean)))
  }))
  result
}
