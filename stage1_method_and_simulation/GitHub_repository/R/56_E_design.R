# Independent, outcome-based configuration selection.
sim_E_defaults <- function() list(
  repetitions=100L, null_counts=c(8L,16L,40L,80L),
  dependence=c("Original","DM","DH"), include_baseline=TRUE,
  tuning_sets=30000L, validation_folds=3L,
  kernel_candidates=data.frame(M=c(96L,144L,144L,144L,192L,288L),
    kappa=c(1e-4,3e-5,1e-4,3e-4,1e-4,1e-4)),
  lambda=c(.3,1,3), manual_candidates=NULL,
  nuisance_knots=c(.25,.5,.75), nuisance_boundary=c(.005,.995),
  nuisance_kernel_rank=24L, fallback_to_fixed=TRUE)

sim_E_candidates <- function(cfg,resolved) {
  p<-resolved$parameters;e<-cfg$E
  if(is.null(e$manual_candidates)) {
    d<-kdr_bind(lapply(e$lambda,function(v){x<-e$kernel_candidates;x$lambda<-v;x}))
  } else d<-e$manual_candidates
  if(!all(c("M","kappa","lambda")%in%names(d)))stop("E candidates require M, kappa and lambda.")
  d<-d[,c("M","kappa","lambda"),drop=FALSE]
  d<-rbind(d,as.data.frame(p[c("M","kappa","lambda")]))
  if(!is.numeric(d$M)||any(!is.finite(d$M)|d$M!=floor(d$M)))stop("E center counts must be finite integers.")
  d<-unique(d);d<-d[order(d$M,d$kappa,d$lambda),,drop=FALSE];rownames(d)<-NULL
  d$M<-as.integer(d$M);d$N<-p$N;d$L<-p$L
  d$candidate_id<-sprintf("E%02d",seq_len(nrow(d)))
  d$is_fixed<-d$M==p$M & d$kappa==p$kappa & d$lambda==p$lambda
  d$kernel_id<-vapply(seq_len(nrow(d)),function(i)substr(object_id(as.list(d[i,c("L","M","kappa")])),1,12),character(1))
  d
}
sim_E_scenarios <- function(cfg) {
  if(!is.numeric(cfg$E$null_counts)||any(!is.finite(cfg$E$null_counts)|cfg$E$null_counts<0|
      cfg$E$null_counts!=floor(cfg$E$null_counts)))stop("E null counts must be nonnegative integers.")
  d<-list()
  if(isTRUE(cfg$E$include_baseline))d[[length(d)+1L]]<-data.frame(
    scenario_id="E0",scenario_label="Eight exposures",null_count=0L,dependence="Original")
  for(q in cfg$E$null_counts)d[[length(d)+1L]]<-data.frame(
    scenario_id=paste0("EQ",q),scenario_label=paste(q,"null exposures"),null_count=as.integer(q),dependence="Original")
  for(id in setdiff(cfg$E$dependence,"Original")) {
    x<-sim_D_levels(cfg);i<-match(id,x$correlation_id)
    if(is.na(i))stop("Unknown E dependence ID: ",id)
    d[[length(d)+1L]]<-data.frame(scenario_id=paste0("E",id),
      scenario_label=paste(x$correlation_label[i],"correlation"),null_count=0L,dependence=id)
  }
  d<-kdr_bind(d)
  if(!nrow(d)||anyDuplicated(d$scenario_id)||any(d$null_count<0|d$null_count!=floor(d$null_count)))
    stop("E scenarios must be nonempty with unique IDs and integer null counts.")
  d
}
sim_plan_E <- function(cfg,resolved) {
  e<-cfg$E
  stopifnot(e$repetitions>=1,e$repetitions==floor(e$repetitions),e$validation_folds>=2,
    e$validation_folds==floor(e$validation_folds),e$tuning_sets>=300,
    e$tuning_sets==floor(e$tuning_sets),e$nuisance_kernel_rank>=1,
    all(e$nuisance_knots>0 & e$nuisance_knots<1),length(e$nuisance_boundary)==2L,
    e$nuisance_boundary[1]<min(e$nuisance_knots),e$nuisance_boundary[2]>max(e$nuisance_knots))
  cand<-sim_E_candidates(cfg,resolved);sc<-sim_E_scenarios(cfg);rows<-list()
  for(i in seq_len(nrow(sc)))for(r in seq_len(e$repetitions))for(strategy in c("Fixed","Adaptive")) {
    a<-as.data.frame(resolved$parameters);a<-cbind(a,sc[i,,drop=FALSE]);a$p<-8L+a$null_count
    a$replicate<-as.integer(r);a$strategy<-strategy;a$config_id<-paste(a$scenario_id,strategy,sep="_")
    a$method<-"KDR";a$scenario_gf<-"G1";a$module<-"E";a$batch<-paste0("application_selection_",a$scenario_id)
    a$key<-paste0(a$config_id,"_r",sprintf("%03d",r));a$kernel_key<-paste0(a$scenario_id,"_r",sprintf("%03d",r))
    a$aux_key<-a$kernel_key;rows[[length(rows)+1L]]<-a
  }
  tasks<-kdr_bind(rows);test<-tasks[rep(1,nrow(cand)),,drop=FALSE]
  test[sim_parameter_names]<-cand[sim_parameter_names];sim_validate(cfg,test)
  list(tasks=tasks,design=list(configurations=cand,membership=data.frame(),scenarios=sc),
    cv_fits_per_repetition=nrow(cand)*e$validation_folds)
}
sim_E_config <- function(cfg,scenario,parameters) {
  if(scenario$dependence!="Original")cfg<-sim_D_level_config(cfg,list(correlation_id=scenario$dependence))
  cfg$batch<-paste0("application_selection_",scenario$scenario_id)
  for(k in c("N","L","M"))cfg[[k]]<-as.integer(parameters[[k]])
  cfg$kernel$kappa<-parameters$kappa;cfg$penalty<-list(a=parameters$lambda,reference_N=parameters$N)
  cfg$nuisance$data_derived<-TRUE;cfg
}
sim_E_choose <- function(scores,candidates,folds,fallback=TRUE) {
  d<-kdr_bind(lapply(candidates$candidate_id,function(id) {
    z<-scores[scores$candidate_id==id,,drop=FALSE]
    good<-nrow(z)==folds && !anyDuplicated(z$fold) && all(z$ok & is.finite(z$loss_sum) & z$n>0)
    data.frame(candidate_id=id,eligible=good,valid_folds=sum(z$ok & is.finite(z$loss_sum)),
      score=if(good)sum(z$loss_sum)/sum(z$n) else NA_real_)
  }))
  ix<-which(d$eligible)
  if(!length(ix)) {
    if(!fallback)stop("No E candidate completed every validation fold.")
    id<-candidates$candidate_id[candidates$is_fixed];used_fallback<-TRUE
  } else {id<-d$candidate_id[ix[which.min(d$score[ix])]];used_fallback<-FALSE}
  list(candidate=candidates[match(id,candidates$candidate_id),,drop=FALSE],
    ranking=merge(candidates,d,by="candidate_id",sort=FALSE),fallback=used_fallback,
    selection_ok=!used_fallback)
}
