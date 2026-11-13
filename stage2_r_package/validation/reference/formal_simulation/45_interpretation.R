# Effect contrasts, exact contribution allocation and interaction detection.
sim_contribution_rule <- function() "active_reference_shapley_v1"

sim_active_contributions <- function(fit,s,ev,cfg) {
  active<-seq_along(exposure_names);b<-ev$reference;p<-length(b)
  if(p<length(active))stop("The contribution target requires eight active exposures.")
  ids<-ev$interpretation$ids
  if(is.null(ids))ids<-seq_len(min(nrow(ev$individual),cfg$interpretation$contribution_points))
  A<-ev$individual[ids,,drop=FALSE]
  if(ncol(A)!=p||!nrow(A)||any(!is.finite(A))||any(!is.finite(b)))
    stop("Invalid contribution evaluation combinations or reference.")
  inactive<-setdiff(seq_len(p),active)
  if(length(inactive))A[,inactive]<-matrix(b[inactive],nrow(A),length(inactive),byrow=TRUE)
  if(s$kind=="kernel") {
    # Restrict the fitted Gaussian function to the reference background.
    # Constant inactive-coordinate factors multiply the original center coefficients.
    sa<-s
    if(length(inactive)) {
      factor<-rep(1,nrow(s$centers))
      for(j in inactive)factor<-factor*exp(-s$rho[j]*((b[j]-s$center[j])/s$scale[j]-s$centers[,j])^2)
      sa$B<-sweep(s$B,1,factor,"*")
      sa$centers<-s$centers[,active,drop=FALSE]
      for(n in c("rho","center","scale"))sa[[n]]<-s[[n]][active]
    }
    phi<-sim_shapley_kernel(A[,active,drop=FALSE],b[active],b[active],sa,fit$beta,
      cfg$interpretation$chunk_points)$total
  } else {
    ref<-matrix(b,nrow(A),p,byrow=TRUE);base<-drop(kdr_basis(ref,s)%*%fit$beta)
    phi<-matrix(0,nrow(A),length(active))
    for(j in active){X<-ref;X[,j]<-A[,j];phi[,j]<-drop(kdr_basis(X,s)%*%fit$beta)-base}
  }
  delta<-drop((kdr_basis(A,s)-kdr_basis(matrix(b,nrow(A),p,byrow=TRUE),s))%*%fit$beta)
  identity_error<-max(abs(rowSums(phi)-delta))
  if(any(!is.finite(phi))||identity_error>1e-8*max(1,max(abs(delta))))
    stop("Active-exposure contribution identity failed.")
  truth<-sim_true_contributions(A,b,cfg)[,active,drop=FALSE]
  true_delta<-kdr_truth(A,cfg)-kdr_truth(matrix(b,1,p),cfg)
  if(max(abs(rowSums(truth)-true_delta))>1e-8*max(1,max(abs(true_delta))))
    stop("True active-exposure contribution identity failed.")
  importance<-colMeans(abs(phi));true_importance<-colMeans(abs(truth))
  share<-if(sum(importance)>0)100*importance/sum(importance) else rep(NA_real_,length(active))
  true_share<-if(sum(true_importance)>0)100*true_importance/sum(true_importance) else rep(NA_real_,length(active))
  list(contributions=data.frame(variable=exposure_names,true_active=TRUE,
      true_share=true_share,estimated_share=share,absolute_contribution=importance),
    contribution_identity_error=identity_error,contribution_scope=sim_contribution_rule(),
    contribution_reference=b,contribution_evaluation_ids=ids)
}

sim_true_contributions <- function(A,b,cfg) {
  p<-ncol(A);z<-matrix(0,nrow(A),p);t<-cfg$baseline$truth
  for(j in 1:6)z[,j]<-t$pollutant_coefficients[j]*(kdr_anchor(A[,j],t$anchors[[j]])-kdr_anchor(b[j],t$anchors[[j]]))
  temp<-function(v){u<-(v-t$temperature$optimum)/t$temperature$distance
    t$temperature$cold_log_coefficient*pmin(u,0)^2+t$temperature$hot_log_coefficient*pmax(u,0)^2}
  z[,7]<-temp(A[,7])-temp(b[7])
  X<-sweep(sweep(A[,1:8,drop=FALSE],2,cfg$baseline$reference$center),2,cfg$baseline$reference$scale,"/")
  B<-(b[1:8]-cfg$baseline$reference$center)/cfg$baseline$reference$scale
  for(i in seq_len(nrow(t$interactions))) {
    a<-t$interactions[i,];j<-a$q;k<-a$r;c<-a$coefficient
    z[,j]<-z[,j]+c*(X[,j]-B[j])*(X[,k]+B[k])/2
    z[,k]<-z[,k]+c*(X[,k]-B[k])*(X[,j]+B[j])/2
  }
  z
}
sim_OR_grid <- function(ev,cfg) {
  p<-length(ev$reference);b<-ev$reference;A<-B<-list();tab<-list()
  at<-function(q)vapply(seq_len(p),function(j)approx(ev$quantile,ev$joint[,j],q,rule=2)$y,numeric(1))
  add<-function(x,y,id,target,exposure){i<-length(A)+1L;A[[i]]<<-x;B[[i]]<<-y
    tab[[i]]<<-data.frame(contrast=id,target=target,exposure=exposure)}
  for(q in cfg$interpretation$joint_OR_quantiles)add(at(q),b,paste0("P",100*q," vs P50"),"joint","all")
  add(at(.75),at(.25),"P75 vs P25","joint","all")
  for(j in 1:8) {
    if(j==7L)for(q in c(.1,.9)) {
      a<-r<-b;a[j]<-at(q)[j];r[j]<-cfg$baseline$truth$temperature$optimum
      add(a,r,paste0("P",100*q," vs optimum"),"single",exposure_names[j])
    } else {a<-r<-b;a[j]<-at(.75)[j];r[j]<-at(.25)[j];add(a,r,"P75 vs P25","single",exposure_names[j])}
  }
  A<-do.call(rbind,A);B<-do.call(rbind,B);t<-kdr_bind(tab)
  t$truth<-kdr_truth(A,cfg)-kdr_truth(B,cfg);t$true_OR<-exp(t$truth)
  list(A=A,B=B,points=t)
}
sim_interpretation_context <- function(ev,cfg) {
  st<-sim_B_response_settings(list(B_single_points=cfg$interpretation$single_points,
    B_interaction_points=cfg$interpretation$interaction_points),cfg)
  grid<-sim_B_response_grid(ev,cfg,st);ids<-seq_len(min(nrow(ev$individual),cfg$interpretation$contribution_points))
  A<-ev$individual[ids,,drop=FALSE];active<-seq_along(exposure_names)
  inactive<-setdiff(seq_len(ncol(A)),active)
  if(length(inactive))A[,inactive]<-matrix(ev$reference[inactive],nrow(A),length(inactive),byrow=TRUE)
  phi<-sim_true_contributions(A,ev$reference,cfg)[,active,drop=FALSE]
  d<-kdr_truth(A,cfg)-kdr_truth(matrix(ev$reference,1,ncol(A)),cfg)
  stopifnot(max(abs(rowSums(phi)-d))<1e-8*max(1,max(abs(d))))
  v<-colMeans(abs(phi));list(response=grid,OR=sim_OR_grid(ev,cfg),ids=ids,
    true_contribution_share=100*v/sum(v))
}
sim_contrast_prediction <- function(D,b,V,truth,cfg) {
  estimate<-drop(D%*%b);se<-rep(NA_real_,nrow(D))
  if(is.matrix(V)&&all(dim(V)==ncol(D))&&all(is.finite(V))) {
    vv<-rowSums((D%*%V)*D);sc<-max(1,max(abs(vv)))
    if(any(vv < -1e-10*sc))stop("Materially negative contrast variance.")
    se<-sqrt(pmax(vv,0))
  }
  q<-qnorm((1+cfg$evaluation$confidence)/2);lo<-estimate-q*se;hi<-estimate+q*se
  sig<-lo>0|hi<0;signal<-abs(truth)>cfg$evaluation$zero_tolerance
  data.frame(estimate=estimate,SE=se,lower=lo,upper=hi,
    correct_significant=ifelse(signal,sig&sign(estimate)==sign(truth),NA),
    wrong_significant=ifelse(signal,sig&sign(estimate)!=sign(truth),NA),
    covered=lo<=truth&hi>=truth)
}
sim_extended_fit <- function(fit,s,ev,cfg,V) {
  ctx<-ev$interpretation;g<-ctx$response;D<-sim_B_response_operator(g,s)
  response<-cbind(g$points,sim_contrast_prediction(D,fit$beta,V,g$points$truth,cfg))
  i<-g$points$target=="interaction"&g$points$eligible
  sub<-list(points=g$points[i,,drop=FALSE],zero_tolerance=cfg$evaluation$zero_tolerance)
  if(s$kind=="kernel")detection<-kdr_attempt(sim_interaction_infer(D[i,,drop=FALSE],fit$beta,V,sub,cfg$interpretation)) else {
    # An additive response has an exactly zero four-point interaction contrast.
    if(max(abs(D[i,,drop=FALSE]))>1e-9)stop("Additive comparator interaction identity failed.")
    pairs<-kdr_bind(lapply(split(sub$points,sub$points$component),function(z)data.frame(
      component=z$component[1],exposure_1=z$exposure_1[1],exposure_2=z$exposure_2[1],
      evaluation_points=nrow(z),true_interaction=any(abs(z$truth)>sub$zero_tolerance),
      statistic=0,p_value=1,p_holm_pairs=1,raw_detected=FALSE,holm_pairs_detected=FALSE,
      test_ok=TRUE,numerical_accuracy=0,algorithm="Additive model",ifault=0L,effective_covariance_rank=0)))
    z<-sub$points;z$estimate<-0;z$SE<-0;z$p_value<-1;z$p_holm_grid<-1
    z$true_signal<-abs(z$truth)>sub$zero_tolerance;z$point_ok<-TRUE
    z$nominal_correct_direction<-z$holm_grid_correct_direction<-ifelse(z$true_signal,FALSE,NA)
    z$nominal_wrong_direction<-z$holm_grid_wrong_direction<-ifelse(z$true_signal,FALSE,NA)
    detection<-list(ok=TRUE,error="",value=list(pairs=pairs,points=z))
  }
  contribution<-sim_active_contributions(fit,s,ev,cfg)
  og<-ctx$OR;OD<-kdr_basis(og$A,s)-kdr_basis(og$B,s)
  ort<-cbind(og$points,sim_contrast_prediction(OD,fit$beta,V,og$points$truth,cfg))
  ort$OR<-exp(ort$estimate);ort$OR_lower<-exp(ort$lower);ort$OR_upper<-exp(ort$upper)
  ort$OR_absolute_relative_error<-100*abs(expm1(ort$estimate-ort$truth))
  c(list(response=response,OR=ort,detection=detection),contribution)
}
