sim_effect_summary <- function(records,planned) {
  records<-Filter(function(a)isTRUE(a$ok)&&is.data.frame(a$outputs)&&nrow(a$outputs)>0,records)
  if(!length(records))return(list(metrics=data.frame(),points=data.frame(),repetitions=data.frame(),or_errors=list()))
  metrics<-points<-repetitions<-or_errors<-list()
  for(target in c("joint","individual")) {
    parts<-lapply(records,function(a)a$outputs[a$outputs$target==target & a$outputs$eligible,,drop=FALSE])
    keep<-vapply(parts,function(x)nrow(x)>0&&all(is.finite(x$estimate))&&all(is.finite(x$truth)),logical(1))
    parts<-parts[keep];rs<-records[keep];if(!length(parts))next
    first<-parts[[1]];J<-nrow(first);R<-length(parts)
    interval_ids<-if("interval_eligible"%in%names(first))which(first$interval_eligible) else seq_len(J)
    stopifnot(all(vapply(parts,function(x)identical(x$index,first$index)&&max(abs(x$truth-first$truth))<1e-12,logical(1))))
    estimate<-do.call(rbind,lapply(parts,`[[`,"estimate"));se<-do.call(rbind,lapply(parts,`[[`,"se"))
    truth<-first$truth;E<-sweep(estimate,2,truth);bias<-colMeans(E)
    mse<-rowMeans(E^2);rmse<-sqrt(mean(mse));rms_bias<-sqrt(mean(bias^2))
    norm<-sqrt(mean(truth^2));sdj<-apply(estimate,2,sim_sd)
    model_ms<-colMeans(se^2,na.rm=TRUE);model_ms[!is.finite(model_ms)]<-NA_real_
    emp_sd<-sqrt(sim_mean(sdj[interval_ids]^2));model_se<-sqrt(sim_mean(model_ms[interval_ids]))
    relative<-100*abs(expm1(E));cv<-do.call(rbind,lapply(parts,`[[`,"covered"))
    det<-do.call(rbind,lapply(parts,`[[`,"correct_significant"));wrong<-do.call(rbind,lapply(parts,`[[`,"wrong_significant"))
    width<-do.call(rbind,lapply(parts,function(x)x$upper-x$lower))
    rr<-data.frame(replicate=vapply(rs,function(x)x$job$replicate,integer(1)),target=target,
      RMSE=sqrt(mse),mean_OR_error=apply(relative,1,sim_mean),median_OR_error=apply(relative,1,median),
      correct_significant=apply(det,1,sim_mean),wrong_significant=apply(wrong,1,sim_mean),
      coverage=apply(cv,1,sim_mean),mean_width=apply(width,1,sim_mean))
    jack_bias<-if(R>1L)vapply(seq_len(R),function(r)sqrt(mean(((colSums(E)-E[r,])/(R-1L))^2)),numeric(1)) else NA_real_
    jack_rmse<-if(R>1L)sqrt((sum(mse)-mse)/(R-1L)) else NA_real_
    jkse<-function(x)if(length(x)>1)sqrt((length(x)-1)/length(x)*sum((x-mean(x))^2)) else NA_real_
    metrics[[target]]<-data.frame(target=target,planned_repetitions=planned,
      available_effect_repetitions=R,available_interval_repetitions=sum(rowSums(is.finite(se[,interval_ids,drop=FALSE]))==length(interval_ids)),evaluation_points=J,
      interval_evaluation_points=length(interval_ids),nonfinite_OR_errors=sum(!is.finite(relative)),
      RMS_Bias=rms_bias,RMSE=rmse,normalized_RMS_Bias=rms_bias/max(norm,1e-30),normalized_RMSE=rmse/max(norm,1e-30),
      MCSE_RMS_Bias=jkse(jack_bias),MCSE_RMSE=jkse(jack_rmse),RMS_empirical_SD=emp_sd,RMS_model_SE=model_se,
      SE_ratio=if(is.finite(emp_sd)&&emp_sd>0)model_se/emp_sd else NA_real_,
      mean_OR_absolute_relative_error=mean(relative),median_OR_absolute_relative_error=median(relative),
      OR_error_Q25=unname(quantile(relative,.25)),OR_error_Q75=unname(quantile(relative,.75)),OR_error_P95=unname(quantile(relative,.95)),
      MCSE_mean_OR_error=sim_mcse(rr$mean_OR_error),correct_significant_percent=100*sim_mean(det),
      MCSE_correct_significant_percent=100*sim_mcse(rr$correct_significant),wrong_significant_percent=100*sim_mean(wrong),
      coverage_percent=100*sim_mean(cv),MCSE_coverage_percent=100*sim_mcse(rr$coverage),mean_interval_width=sim_mean(width))
    points[[target]]<-data.frame(target=target,index=first$index,q=first$q,truth=truth,mean_estimate=colMeans(estimate),bias=bias,
      empirical_SD=sdj,RMS_model_SE=sqrt(model_ms),coverage_percent=100*apply(cv,2,sim_mean),
      correct_significant_percent=100*apply(det,2,sim_mean),
      lower_repetition_quantile=apply(estimate,2,quantile,.025),upper_repetition_quantile=apply(estimate,2,quantile,.975))
    repetitions[[target]]<-rr;or_errors[[target]]<-as.vector(relative)
  }
  list(metrics=kdr_bind(metrics),points=kdr_bind(points),repetitions=kdr_bind(repetitions),or_errors=or_errors)
}

sim_collect <- function(run) {
  cfg<-readRDS(file.path(run,"config.rds"));pp<-readRDS(file.path(run,"plan.rds"));plan<-pp$tasks
  records<-setNames(lapply(plan$key,function(k)sim_read(file.path(run,"fits",paste0(k,".rds")))),plan$key)
  status<-kdr_bind(lapply(seq_len(nrow(plan)),function(i){
    a<-records[[i]];j<-plan[i,,drop=FALSE]
    data.frame(j,completed=!is.null(a),point_ok=isTRUE(a$ok),interval_ok=isTRUE(a$interval_ok),
      selection_ok=isTRUE(a$selection$valid),interpretation_ok=isTRUE(a$interpretation_ok),
      interaction_ok=isTRUE(a$interpretation$detection$ok)&&all(a$interpretation$detection$value$pairs$test_ok),reused=isTRUE(a$reused),
      error=if(is.null(a))"Pending" else if(is.null(a$error))"" else a$error,
      interval_error=if(is.null(a$interval_error))"" else a$interval_error,
      selection_error=if(is.null(a$selection$error))"" else a$selection$error,
      elapsed_seconds=if(is.null(a$elapsed_seconds))NA_real_ else a$elapsed_seconds,
      stringsAsFactors=FALSE)
  }))
  selection<-kdr_bind(lapply(seq_along(records),function(i){
    a<-records[[i]];if(is.null(a$selection$table)||!nrow(a$selection$table))return(NULL)
    z<-a$selection$table
    cbind(plan[rep(i,nrow(z)),c("config_id","null_count","p","replicate","method","scenario_gf"),drop=FALSE],z)
  }))
  design_cols<-c("module","config_id","method","scenario_gf","null_count","p",sim_parameter_names)
  summaries<-points<-reps<-or_boxes<-list()
  for(ix in sim_groups(plan,design_cols)) {
    key<-plan[ix[1],design_cols,drop=FALSE]
    result<-sim_effect_summary(records[ix],length(ix))
    append_key<-function(x)if(nrow(x))cbind(key[rep(1L,nrow(x)),,drop=FALSE],x) else NULL
    summaries[[length(summaries)+1L]]<-append_key(result$metrics)
    points[[length(points)+1L]]<-append_key(result$points)
    reps[[length(reps)+1L]]<-append_key(result$repetitions)
    for(t in names(result$or_errors)) {
      box_key<-paste(key$config_id,key$method,key$scenario_gf,t,sep="|")
      if(key$module=="C")box_key<-paste(box_key,key$null_count,sep="|")
      or_boxes[[box_key]]<-result$or_errors[[t]]
    }
  }
  resources<-sim_resource_tables(run);peaks<-resources$peaks
  peak_for<-function(kind,key,column) {
    if(!nrow(peaks)||!column%in%names(peaks))return(NA_real_)
    a<-peaks[peaks$kind==kind & peaks$key==key,column];if(length(a))a[1] else NA_real_
  }
  cost<-kdr_bind(lapply(seq_along(records),function(i){
    a<-records[[i]];if(is.null(a)||!isTRUE(a$ok))return(NULL)
    j<-plan[i,,drop=FALSE]
    gm<-if(j$scenario_gf%in%c("G1","G2"))"correct" else "wrong";fw<-j$scenario_gf%in%c("G2","G4")
    auxkey<-paste0(j$aux_key,"_",gm,"_f",as.integer(fw))
    rss<-c(peak_for("fit",j$key,"sampled_working_set_gb"),peak_for("kernel",j$kernel_key,"sampled_working_set_gb"),
      peak_for("auxiliary",auxkey,"sampled_working_set_gb"))
    disk<-sim_read(file.path(run,"resource_checkpoints",paste0("fit_",j$key,".rds")))
    val<-function(n)if(is.null(a[[n]]))NA_real_ else a[[n]]
    data.frame(j,learning_seconds=val("learning_seconds"),auxiliary_seconds=val("auxiliary_seconds"),
      estimation_seconds=val("computation_seconds"),
      component_total_seconds=val("learning_seconds")+val("auxiliary_seconds")+val("computation_seconds"),
      sampled_peak_memory_GB=if(any(is.finite(rss)))max(rss,na.rm=TRUE) else NA_real_,
      cache_GB=if(is.null(disk))NA_real_ else disk$cache_bytes/1024^3,
      reused=isTRUE(a$reused))
  }))
  interpretation<-sim_interpretation_collect(records,plan)
  list(interpretation=interpretation,extended=sim_extended_summary(interpretation),config=cfg,plan=pp,status=status,selection=selection,effects=kdr_bind(summaries),
    points=kdr_bind(points),repetitions=kdr_bind(reps),or_errors=or_boxes,cost=cost,resources=resources)
}

sim_selection_summary <- function(x,status) {
  if(!nrow(x))return(list(overall=data.frame(),variables=data.frame()))
  overall<-variables<-list()
  for(ix in sim_groups(x,c("config_id","null_count","p"))) {
    z<-x[ix,];valid_reps<-names(Filter(function(v)length(v)==z$p[1] && all(!is.na(v)),split(z$selected,z$replicate)))
    z<-z[as.character(z$replicate)%in%valid_reps,,drop=FALSE]
    if(!nrow(z))next
    a<-kdr_bind(lapply(split(z,z$replicate),function(v)data.frame(replicate=v$replicate[1],
      TPR=mean(v$selected[v$true_active]),FPR=if(any(!v$true_active))mean(v$selected[!v$true_active]) else NA_real_)))
    counts<-status[status$config_id==z$config_id[1]&status$null_count==z$null_count[1]&status$method=="KDR"&status$scenario_gf=="G1",]
    ci<-function(v){n<-sum(is.finite(v));m<-sim_mean(v);h<-if(n>1)qt(.975,n-1)*sim_mcse(v) else NA_real_;c(m,m-h,m+h)}
    t<-ci(a$TPR);f<-ci(a$FPR)
    overall[[length(overall)+1L]]<-data.frame(config_id=z$config_id[1],null_count=z$null_count[1],p=z$p[1],
      planned_repetitions=nrow(counts),usable_repetitions=nrow(a),
      mean_TPR=t[1],TPR_MCSE=sim_mcse(a$TPR),TPR_lower=max(0,t[2]),TPR_upper=min(1,t[3]),
      mean_FPR=f[1],FPR_MCSE=sim_mcse(a$FPR),FPR_lower=max(0,f[2]),FPR_upper=min(1,f[3]))
    for(v in exposure_names) {
      q<-z[z$variable==v,];n<-nrow(q);p<-mean(q$selected);crit<-qnorm(.975)
      den<-1+crit^2/n;c0<-(p+crit^2/(2*n))/den;h<-crit*sqrt(p*(1-p)/n+crit^2/(4*n^2))/den
      variables[[length(variables)+1L]]<-data.frame(config_id=z$config_id[1],null_count=z$null_count[1],variable=v,
        usable_repetitions=n,selected=sum(q$selected),rate=p,lower=max(0,c0-h),upper=min(1,c0+h))
    }
  }
  list(overall=kdr_bind(overall),variables=kdr_bind(variables))
}

sim_recommend <- function(result,run) {
  x<-result$effects;if(!nrow(x))return(NULL)
  cfg<-result$config;out<-list()
  for(id in unique(x$config_id)) {
    z<-x[x$config_id==id,];j<-z[z$target=="joint",];a<-z[z$target=="individual",]
    if(nrow(j)!=1||nrow(a)!=1)next
    success<-min(j$available_effect_repetitions/j$planned_repetitions,a$available_effect_repetitions/a$planned_repetitions)
    out[[length(out)+1L]]<-data.frame(config_id=id,N=j$N,L=j$L,M=j$M,kappa=j$kappa,lambda=j$lambda,
      success_fraction=success,available_repetitions=min(j$available_effect_repetitions,a$available_effect_repetitions),
      joint_nrmse=j$normalized_RMSE,individual_nrmse=a$normalized_RMSE,
      joint_nrms_bias=j$normalized_RMS_Bias,individual_nrms_bias=a$normalized_RMS_Bias,
      joint_or_error=j$mean_OR_absolute_relative_error/100,individual_or_error=a$mean_OR_absolute_relative_error/100,
      joint_coverage=j$coverage_percent,individual_coverage=a$coverage_percent)
  }
  d<-kdr_bind(out);if(!nrow(d))return(NULL)
  e<-result$extended$scores
  if(nrow(e)) {
    ec<-c("config_id","interpretation_repetitions","interaction_repetitions",
      "contribution_share_error","interaction_miss","direction_miss","direction_wrong")
    d<-merge(d,e[ec],by="config_id",all.x=TRUE,sort=FALSE)
    d$available_repetitions<-pmin(d$available_repetitions,d$interpretation_repetitions,d$interaction_repetitions)
    d$success_fraction<-pmin(d$success_fraction,d$available_repetitions/cfg$A$repetitions)
  } else for(n in c("contribution_share_error","interaction_miss","direction_miss","direction_wrong"))d[[n]]<-NA_real_
  w<-cfg$A$recommendation$weights
  stopifnot(all(names(w)%in%names(d)),all(w>=0),sum(w)>0)
  d$score<-drop(as.matrix(d[names(w)])%*%(w/sum(w)))
  d$eligible<-d$success_fraction>=cfg$A$recommendation$min_success_fraction &
    d$available_repetitions>=cfg$A$recommendation$min_repetitions & is.finite(d$score)
  d$eligible[is.na(d$eligible)]<-FALSE
  d<-d[order(!d$eligible,d$score),];d$rank<-seq_len(nrow(d))
  eligible<-d[d$eligible,,drop=FALSE]
  recommendation<-if(nrow(eligible))list(config_id=eligible$config_id[1],
    parameters=as.list(eligible[1,sim_parameter_names,drop=FALSE]),
    source_run=basename(run),score=eligible$score[1],rule=cfg$A$recommendation) else NULL
  kdr_csv(d,file.path(run,"analysis","tables","configuration_ranking.csv"))
  kdr_save(recommendation,file.path(run,"analysis","recommended_configuration.rds"))
  if(!is.null(recommendation))dput(recommendation,file=file.path(run,"analysis","recommended_configuration.R"))
  list(ranking=d,recommendation=recommendation)
}
