# Quadratic-form tests of fixed-background interaction contrasts.
sim_interaction_tail <- function(score,weights,s) {
  w<-weights[weights>0]
  if(!length(w))return(list(p=NA_real_,accuracy=NA_real_,algorithm="Degenerate covariance",ifault=NA_integer_))
  if(!is.finite(score)||score<0)stop("Invalid nonnegative quadratic statistic.")
  if(max(w)-min(w)<=1e-12*max(w))return(list(p=pchisq(score/mean(w),length(w),lower.tail=FALSE),
    accuracy=10*.Machine$double.eps,algorithm="Scaled chi-square",ifault=0L))
  norm<-sum(w);q<-score/norm;w<-w/norm
  if(q==0)return(list(p=1,accuracy=0,algorithm="Zero statistic",ifault=0L))
  d<-suppressWarnings(CompQuadForm::davies(q,w,lim=s$tail_max_terms,acc=s$tail_accuracy))
  if(d$ifault!=0L||!is.finite(d$Qq)||d$Qq < -s$tail_accuracy || d$Qq > 1+s$tail_accuracy)
    d<-suppressWarnings(CompQuadForm::davies(q,w,lim=s$tail_retry_terms,acc=s$tail_accuracy))
  if(d$ifault!=0L||!is.finite(d$Qq)||d$Qq < -s$tail_accuracy || d$Qq > 1+s$tail_accuracy)
    return(list(p=NA_real_,accuracy=NA_real_,algorithm="Davies failed",ifault=d$ifault))
  list(p=pmin(1,pmax(0,d$Qq)),accuracy=s$tail_accuracy,algorithm="Davies",ifault=d$ifault)
}

sim_interaction_infer <- function(D,beta,V,grid,s) {
  J<-nrow(D);K<-ncol(D)
  if(length(beta)!=K||any(!is.finite(beta))||any(!is.finite(D)))stop("Invalid saved coefficients or interaction operator.")
  if(!is.matrix(V)||any(dim(V)!=K)||any(!is.finite(V)))stop("Missing or invalid saved coefficient covariance.")
  V<-kdr_sym(V);e<-eigen(V,symmetric=TRUE);scale<-max(abs(e$values))
  if(scale<=0||min(e$values)< -1e-8*scale)stop("Coefficient covariance is zero or materially indefinite.")
  L<-sweep(e$vectors,2,sqrt(pmax(e$values,0)),"*")
  estimate<-drop(D%*%beta);Z<-D%*%L;variance<-rowSums(Z^2);se<-sqrt(variance)
  if(any(!is.finite(estimate))||any(!is.finite(se)))stop("Nonfinite interaction estimates or variances.")
  point_ok<-is.finite(se)&se>0
  pv<-rep(NA_real_,J);pv[point_ok]<-2*pnorm(abs(estimate[point_ok]/se[point_ok]),lower.tail=FALSE)
  adj<-p.adjust(pv,method="holm",n=J)
  q<-qnorm(1-s$alpha/2);points<-grid$points
  points$estimate<-estimate;points$SE<-se;points$lower<-estimate-q*se;points$upper<-estimate+q*se
  points$lower[!point_ok]<-points$upper[!point_ok]<-NA_real_
  points$p_value<-pv;points$p_holm_grid<-adj
  points$nominal_detected<-pv<s$alpha;points$holm_grid_detected<-adj<s$alpha
  signal<-abs(points$truth)>grid$zero_tolerance;correct<-sign(estimate)==sign(points$truth)
  for(mode in c("nominal","holm_grid")) {
    selected<-points[[paste0(mode,"_detected")]]
    points[[paste0(mode,"_correct_direction")]]<-ifelse(signal,selected&correct,NA)
    points[[paste0(mode,"_wrong_direction")]]<-ifelse(signal,selected&!correct,NA)
    points[[paste0(mode,"_false_positive")]]<-ifelse(!signal,selected,NA)
  }
  points$true_signal<-signal;points$point_ok<-point_ok
  pairs<-list();weights<-list()
  for(pair in unique(points$component)) {
    ix<-which(points$component==pair);n<-length(ix);A<-Z[ix,,drop=FALSE]/sqrt(n)
    # The smaller Gram matrix has the same nonzero mixture weights.
    gram<-if(n<=K)tcrossprod(A) else crossprod(A)
    ee<-eigen(kdr_sym(gram),symmetric=TRUE,only.values=TRUE)$values
    if(min(ee)< -1e-8*max(abs(ee)))stop("Invalid interaction covariance spectrum.")
    w<-pmax(ee,0);score<-mean(estimate[ix]^2);tail<-sim_interaction_tail(score,w,s);weights[[pair]]<-w
    pairs[[pair]]<-data.frame(component=pair,exposure_1=points$exposure_1[ix[1]],exposure_2=points$exposure_2[ix[1]],
      evaluation_points=n,true_interaction=any(signal[ix]),statistic=score,
      p_value=tail$p,numerical_accuracy=tail$accuracy,algorithm=tail$algorithm,ifault=tail$ifault,
      effective_covariance_rank=sum(ee>max(ee)*1e-10),stringsAsFactors=FALSE)
  }
  pairs<-kdr_bind(pairs);P<-nrow(pairs)
  # Check numerical integration uncertainty against both decision thresholds.
  decisions<-function(p,a) {
    lo<-pmax(0,p-a);hi<-pmin(1,p+a)
    list(raw=(lo<s$alpha)!=(hi<s$alpha),
      holm=(p.adjust(lo,"holm",n=P)<s$alpha)!=(p.adjust(hi,"holm",n=P)<s$alpha))
  }
  for(attempt in 1:2) {
    amb<-decisions(pairs$p_value,pairs$numerical_accuracy)
    bad<-which((amb$raw|amb$holm)&pairs$algorithm=="Davies")
    if(!length(bad))break
    fine<-s;fine$tail_accuracy<-s$tail_accuracy/100^attempt;fine$tail_max_terms<-s$tail_retry_terms
    for(i in which(pairs$algorithm=="Davies")) {
      a<-sim_interaction_tail(pairs$statistic[i],weights[[pairs$component[i]]],fine)
      pairs$p_value[i]<-a$p;pairs$numerical_accuracy[i]<-a$accuracy;pairs$algorithm[i]<-a$algorithm;pairs$ifault[i]<-a$ifault
    }
  }
  amb<-decisions(pairs$p_value,pairs$numerical_accuracy)
  pairs$p_holm_pairs<-p.adjust(pairs$p_value,"holm",n=P)
  pairs$raw_detected<-ifelse(amb$raw,NA,pairs$p_value<s$alpha)
  pairs$holm_pairs_detected<-ifelse(amb$holm,NA,pairs$p_holm_pairs<s$alpha)
  pairs$test_ok<-is.finite(pairs$p_value)&!is.na(pairs$raw_detected)&!is.na(pairs$holm_pairs_detected)
  gate<-pairs$holm_pairs_detected[match(points$component,pairs$component)]
  points$pair_detected<-gate
  points$pair_gated_nominal_correct<-ifelse(signal,gate&points$nominal_correct_direction,NA)
  points$pair_gated_nominal_wrong<-ifelse(signal,gate&points$nominal_wrong_direction,NA)
  list(pairs=pairs,points=points,covariance_negative_eigenvalue=min(e$values),
    inference_scope="Centered-Gaussian quadratic-form inference conditional on the independently learned kernel, evaluated for its operating characteristics under regularization.")
}

