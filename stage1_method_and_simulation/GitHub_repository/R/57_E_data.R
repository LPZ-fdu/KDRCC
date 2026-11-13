# Training-only nuisance dictionary and grouped validation data.
sim_E_nuisance_spec <- function(d,s,cfg) {
  thin<-s[c("kind","center","scale","centers","rho","M","B")]
  thin$B<-thin$B[,seq_len(min(ncol(thin$B),cfg$E$nuisance_kernel_rank)),drop=FALSE]
  list(splines=lapply(seq_len(ncol(d$A)),function(j)list(
    knots=unique(as.numeric(quantile(d$A[,j],cfg$E$nuisance_knots))),
    boundary=as.numeric(quantile(d$A[,j],cfg$E$nuisance_boundary)))),kernel=thin)
}
sim_E_nuisance_raw <- function(A,W,spec) {
  B<-do.call(cbind,lapply(seq_len(ncol(A)),function(j)splines::ns(A[,j],
    knots=spec$splines[[j]]$knots,Boundary.knots=spec$splines[[j]]$boundary,intercept=FALSE)))
  cbind(W=W,W2=W^2,B,kdr_basis(A,spec$kernel))
}
sim_E_take_groups <- function(d,ids) {
  ids<-sort(as.integer(ids));take<-d$group%in%ids
  out<-d
  for(k in c("A","W","date","event"))out[[k]]<-if(is.matrix(d[[k]]))d[[k]][take,,drop=FALSE] else d[[k]][take]
  out$truth<-out$probability<-NULL
  out$group<-match(d$group[take],ids);out$fold<-d$fold[ids];out$n<-length(ids)
  out$pairs<-if(length(ids))kdr_pairs(out) else NULL;out
}
sim_E_raw <- function(cfg,cache,N,p,r,role="tuning",fold=NULL,outer=NULL,validation=FALSE) {
  raw<-sim_raw(cfg,cache,N,p,r,role)
  function(FUN,progress=FALSE)raw(function(d){
    if(!is.null(fold)) {
      labels<-outer[d$offset+seq_len(d$n)]
      ids<-which(if(validation)labels==fold else labels!=fold)
      if(!length(ids))return(invisible(NULL))
      d<-sim_E_take_groups(d,ids)
    }
    d$truth<-d$probability<-NULL;FUN(d)
  },progress)
}
sim_E_score <- function(fit,s,raw) {
  loss<-0;n<-0L
  raw(function(d){
    a<-kdr_cl_block(c(fit$beta,fit$gamma),cbind(kdr_basis(d$A,s),d$W,d$W^2),d$event,d$group,FALSE)
    loss<<-loss+a$loss;n<<-n+d$n
  })
  if(n==0L||!is.finite(loss))stop("Invalid validation likelihood.")
  list(loss_sum=loss,n=n,mean_loss=loss/n)
}
sim_E_prepare_data <- function(cfg,cache,job) {
  sizes<-c(learning=job$L,tuning=cfg$E$tuning_sets,estimation=job$N)
  total<-sum(ceiling(sizes/cfg$execution$block_groups));done<-0L
  sim_major(1,"Independent learning, configuration and final-estimation samples")
  for(role in names(sizes))for(b in seq_len(ceiling(sizes[[role]]/cfg$execution$block_groups))) {
    path<-sim_pool_file(cache,role,job$replicate,b)
    if(is.null(sim_read(path))) {
      d<-kdr_generate(cfg$execution$block_groups,8L,cfg,kdr_seed(cfg,"data_pool",role,job$replicate,b))
      d$truth<-d$probability<-NULL;kdr_save(d,path)
    }
    done<-done+1L;sim_status(paste("Preparing",role,"block",b),done,total,TRUE)
  }
}
sim_E_prepare_evaluations <- function(cfg,run,scenarios) {
  base<-cfg;base$baseline$dgp$copula_eta<-NULL;base$batch<-"application_evaluation"
  ref<-sim_evaluation(base,8L,file.path(run,"evaluation","common_reference.rds"))
  for(i in seq_len(nrow(scenarios))) {
    a<-as.list(scenarios[i,,drop=FALSE]);path<-file.path(run,"evaluation",paste0(a$scenario_id,".rds"))
    if(!is.null(sim_read(path)))next
    cf<-base;if(a$dependence!="Original")cf<-sim_D_level_config(cf,list(correlation_id=a$dependence))
    cf$batch<-paste0("application_evaluation_",a$scenario_id);p<-8L+a$null_count
    d<-kdr_generate(cf$evaluation$groups,p,cf,kdr_seed(cf,"evaluation"),events=FALSE)
    out<-ref
    if(p>8L) {
      out$reference<-c(ref$reference,apply(d$A[,-(1:8),drop=FALSE],2,median))
      extra<-t(vapply(out$quantile,function(q)apply(d$A[,-(1:8),drop=FALSE],2,quantile,q),numeric(p-8L)))
      out$joint<-cbind(ref$joint,extra)
    }
    rows<-match(seq_len(cf$evaluation$function_points),d$group);out$individual<-d$A[rows,,drop=FALSE]
    truth0<-kdr_truth(matrix(out$reference,1,p),cf)
    out$true_joint<-kdr_truth(out$joint,cf)-truth0;out$true_individual<-kdr_truth(out$individual,cf)-truth0
    out$overall_correlation<-cor(d$A);out$interpretation<-sim_interpretation_context(out,cf)
    # Geometric support diagnostic on the eight active exposure coordinates.
    mu<-colMeans(d$A[,1:8]);S<-cov(d$A[,1:8]);inv<-kdr_inverse(S)
    dist2<-function(A){Z<-sweep(A[,1:8,drop=FALSE],2,mu);rowSums((Z%*%inv)*Z)}
    cutoff<-as.numeric(quantile(dist2(d$A),.95));out$support<-data.frame(index=seq_len(nrow(out$joint)),
      quantile=out$quantile,distance_squared=dist2(out$joint),reference_P95=cutoff,
      outside_reference_P95=dist2(out$joint)>cutoff)
    kdr_save(out,path)
  }
}
