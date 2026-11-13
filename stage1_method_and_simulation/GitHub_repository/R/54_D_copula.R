# Calibrated dependence of final lag-defined exposures, using fixed marginal maps.
sim_D_common_score <- function(group,date,phi,seed) with_seed(stage_seed(seed,"D_common_score"),{
  u<-rnorm(length(group));day<-as.numeric(date)
  for(i in seq_along(u)[-1L])if(group[i]==group[i-1L]) {
    a<-phi^(day[i]-day[i-1L]);u[i]<-a*u[i-1L]+sqrt(1-a*a)*u[i]
  }
  u
})
sim_D_gaussian_scores <- function(A,maps) {
  z<-vapply(seq_len(8L),function(j)qnorm(approx(maps$values[[j]],maps$probability,
    A[,j],rule=2,ties="ordered")$y),numeric(nrow(A)))
  colnames(z)<-exposure_names;z
}
sim_D_inverse_scores <- function(z,maps) {
  a<-vapply(seq_len(8L),function(j)approx(maps$probability,maps$values[[j]],
    pnorm(z[,j]),rule=2,ties="ordered")$y,numeric(nrow(z)))
  colnames(a)<-exposure_names;a
}
sim_D_transform <- function(z,u,maps,eta,direction) {
  stopifnot(is.finite(eta),eta>=0,eta<1,length(direction)==8L,all(abs(direction)==1))
  sim_D_inverse_scores(sqrt(1-eta)*z+sqrt(eta)*outer(u,direction),maps)
}
sim_D_apply_copula <- function(A,group,date,cfg,seed) {
  eta<-cfg$baseline$dgp$copula_eta
  if(is.null(eta)||eta==0)return(A)
  maps<-cfg$D$copula$maps
  if(is.null(maps))stop("D marginal maps have not been calibrated.")
  z<-sim_D_gaussian_scores(A,maps)
  u<-sim_D_common_score(group,date,cfg$baseline$dgp$phi,seed)
  sim_D_transform(z,u,maps,eta,cfg$baseline$dgp$direction)
}
sim_D_cor_summary <- function(A) {
  R<-cor(A);v<-abs(R[upper.tri(R)])
  c(mean_absolute=mean(v),minimum_absolute=min(v),maximum_absolute=max(v))
}
sim_D_check_target <- function(summary,level) {
  all(is.finite(summary)) && abs(summary["mean_absolute"]-level$target_mean)<=level$mean_tolerance &&
    summary["minimum_absolute"]>=level$min_abs
}
sim_D_prepare_generator <- function(cfg,run) {
  path<-file.path(run,"calibration","copula.rds")
  old<-sim_read(path)
  identity<-object_id(list(baseline=cfg$baseline,seed=cfg$seed,schema=cfg$schema,
    levels=sim_D_levels(cfg)[,c("correlation_id","correlation_label","target_mean","min_abs","mean_tolerance")],
    settings=cfg$D$calibration))
  if(!is.null(old)) {
    if(!identical(identity,old$identity))stop("D calibration identity mismatch.")
    cfg$D$copula<-old;return(cfg)
  }
  cf<-cfg;cf$baseline$dgp$copula_eta<-NULL;cf$baseline$dgp$correlation_mix<-NULL
  cf$batch<-"correlation";settings<-cfg$D$calibration
  getdata<-function(role,n) {
    seed<-kdr_seed(cf,"copula_calibration",role)
    kdr_log("D calibration:",role,"| matched sets",n)
    d<-kdr_generate(n,8L,cf,seed,events=FALSE)
    list(d=d,u=sim_D_common_score(d$group,d$date,cf$baseline$dgp$phi,seed))
  }
  ref<-getdata("marginals",settings$reference_groups)$d
  maps<-list(values=lapply(1:8,function(j)sort(ref$A[,j])),
    probability=(seq_len(nrow(ref$A))-.5)/nrow(ref$A))
  base<-getdata("training",settings$calibration_groups)
  valid<-getdata("validation",settings$validation_groups)
  Z<-sim_D_gaussian_scores(base$d$A,maps);V<-sim_D_gaussian_scores(valid$d$A,maps)
  levels<-sim_D_levels(cf);rows<-marginals<-list()
  for(i in seq_len(nrow(levels))) {
    a<-levels[i,,drop=FALSE]
    make<-function(eta,z=Z,u=base$u)sim_D_transform(z,u,maps,eta,cf$baseline$dgp$direction)
    f<-function(eta)sim_D_cor_summary(make(eta))["mean_absolute"]-a$target_mean
    bracket<-c(0,settings$maximum_mix)
    if(f(bracket[1])*f(bracket[2])>0)stop("D target mean is outside the attainable range: ",a$correlation_id)
    eta<-uniroot(f,bracket,tol=settings$root_tolerance)$root
    train<-make(eta);check<-make(eta,V,valid$u)
    st<-sim_D_cor_summary(train);sv<-sim_D_cor_summary(check)
    if(!sim_D_check_target(st,a)||!sim_D_check_target(sv,a))
      stop("D correlation target/floor failed independent calibration: ",a$correlation_id,
        ". Review the calibration settings before fitting.")
    a$correlation_mix<-eta
    rows[[i]]<-cbind(a,
      as.data.frame(as.list(setNames(st,paste0("training_",names(st))))),
      as.data.frame(as.list(setNames(sv,paste0("validation_",names(sv))))))
    marginals[[i]]<-data.frame(correlation_id=a$correlation_id,variable=exposure_names,
      original_mean=colMeans(valid$d$A),generated_mean=colMeans(check),
      original_SD=apply(valid$d$A,2,sd),generated_SD=apply(check,2,sd))
    kdr_log(a$correlation_label,"| achieved mean absolute r",format(sv[1],digits=4),
      "| minimum absolute r",format(sv[2],digits=4))
  }
  out<-list(identity=identity,maps=maps,levels=kdr_bind(rows),marginals=kdr_bind(marginals),
    reference_groups=settings$reference_groups,settings=settings)
  kdr_save(out,path);kdr_csv(out$levels,file.path(run,"calibration","achieved_correlations.csv"))
  kdr_csv(out$marginals,file.path(run,"calibration","marginal_distribution_checks.csv"))
  cfg$D$copula<-out;cfg
}
