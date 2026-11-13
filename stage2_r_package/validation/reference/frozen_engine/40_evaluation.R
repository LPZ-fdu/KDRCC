kdr_outputs <- function(fit,s,ev,cfg,V=NULL) {
  reference<-drop(kdr_basis(matrix(ev$reference,1,length(ev$reference)),s))
  operators<-list(joint=sweep(kdr_basis(ev$joint,s),2,reference),
    individual=sweep(kdr_basis(ev$individual,s),2,reference))
  truths<-list(joint=ev$true_joint,individual=ev$true_individual)
  crit<-qnorm(1-(1-cfg$evaluation$confidence)/2);out<-list()
  for(name in names(operators)) {
    D<-operators[[name]];estimate<-drop(D%*%fit$beta);truth<-truths[[name]]
    eligible<-if(name=="joint")abs(ev$quantile-.5)>1e-12 else rep(TRUE,length(truth))
    interval_eligible<-if(name=="joint")eligible else seq_along(truth)%in%ev$coverage_ids
    se<-rep(NA_real_,length(truth))
    if(!is.null(V))se[interval_eligible]<-sqrt(pmax(rowSums((D[interval_eligible,,drop=FALSE]%*%V)*D[interval_eligible,,drop=FALSE]),0))
    low<-estimate-crit*se;up<-estimate+crit*se
    signal<-abs(truth)>cfg$evaluation$zero_tolerance
    correct<-ifelse(!signal,NA,(truth>0 & low>0)|(truth<0 & up<0))
    wrong<-ifelse(!signal,NA,(truth>0 & up<0)|(truth<0 & low>0))
    out[[name]]<-data.frame(target=name,index=seq_along(truth),
      q=if(name=="joint")ev$quantile else NA_real_,truth=truth,estimate=estimate,se=se,
      lower=low,upper=up,eligible=eligible,interval_eligible=interval_eligible,
      covered=ifelse(interval_eligible & is.finite(se),low<=truth & up>=truth,NA),
      correct_significant=correct,wrong_significant=wrong,
      or_absolute_relative_error=100*abs(expm1(estimate-truth)))
  }
  kdr_bind(out)
}

kdr_quadratic_tail <- function(score,weights,cfg,seed) {
  w<-weights[weights>0]
  if(!length(w))return(list(p=if(score<=0)1 else 0,method="degenerate",accuracy=0,ifault=0L,B=0L))
  if(max(w)-min(w)<=1e-13*max(w))return(list(
    p=pchisq(score/mean(w),df=length(w),lower.tail=FALSE),method="Scaled chi-square",
    accuracy=10*.Machine$double.eps,ifault=0L,B=0L))
  norm<-sum(w);w<-w/norm;q<-score/norm
  if(q<=0)return(list(p=1,method="boundary",accuracy=0,ifault=0L,B=0L))
  d<-CompQuadForm::davies(q,w,lim=cfg$selection$max_terms,acc=cfg$selection$accuracy)
  if(d$ifault!=0L || !is.finite(d$Qq) || d$Qq<0 || d$Qq>1)
    d<-CompQuadForm::davies(q,w,lim=cfg$selection$retry_terms,acc=cfg$selection$accuracy)
  if(d$ifault==0L && is.finite(d$Qq) && d$Qq>=0 && d$Qq<=1)
    return(list(p=d$Qq,method="Davies",accuracy=cfg$selection$accuracy,ifault=d$ifault,B=0L))
  if(!cfg$selection$fallback)return(list(p=NA_real_,method="failed",accuracy=NA_real_,ifault=d$ifault,B=0L))
  B<-max(cfg$selection$fallback_min_draws,
    ceiling(cfg$selection$expected_tail_exceedances*cfg$total_dimension/cfg$selection$alpha)-1L)
  count<-with_seed(seed,{
    done<-0;exceed<-0
    while(done<B) {
      n<-min(cfg$selection$fallback_batch,B-done)
      values<-drop(matrix(rnorm(n*length(w))^2,nrow=n)%*%w)
      exceed<-exceed+sum(values>=q);done<-done+n
      kdr_heartbeat(paste("Quadratic-tail fallback draws",done,"of",B))
    };exceed
  })
  list(p=(1+count)/(B+1),method="Monte Carlo fallback",accuracy=NA_real_,ifault=d$ifault,B=B)
}

kdr_selection <- function(fit,s,ev,cfg,V,job) {
  p<-length(ev$reference);cfg$total_dimension<-p
  if(is.null(V))return(list(table=data.frame(),valid=FALSE,error="Covariance unavailable"))
  A<-ev$individual[seq_len(min(cfg$selection$evaluation_points,nrow(ev$individual))),,drop=FALSE]
  L<-kdr_vroot(V);out<-weights<-list()
  for(j in seq_len(p)) {
    D<-kdr_derivative(A,s,j);Q<-crossprod(D)/nrow(D)
    score<-max(0,drop(crossprod(fit$beta,Q%*%fit$beta)))
    e<-eigen(kdr_sym(crossprod(L,Q%*%L)),symmetric=TRUE,only.values=TRUE)$values
    if(min(e) < -1e-9*max(abs(e)))stop("Indefinite quadratic-form weights")
    weights[[j]]<-pmax(e,0)
    a<-kdr_quadratic_tail(score,weights[[j]],cfg,kdr_seed(cfg,"selection",job$key,j))
    out[[j]]<-data.frame(variable=kdr_names(p)[j],true_active=j<=8L,score=score,
      p_value=a$p,algorithm=a$method,absolute_accuracy=a$accuracy,ifault=a$ifault,draws=a$B)
    kdr_bar("Variable tests",j,p)
  }
  tab<-kdr_bind(out)
  # Check that the numerical tolerance cannot change the Holm decision.
  for(attempt in 1:2) {
    acc<-ifelse(is.finite(tab$absolute_accuracy),tab$absolute_accuracy,0)
    lo<-p.adjust(pmax(tab$p_value-acc,0),cfg$selection$adjust)<cfg$selection$alpha
    hi<-p.adjust(pmin(tab$p_value+acc,1),cfg$selection$adjust)<cfg$selection$alpha
    bad<-which(!is.na(lo)&!is.na(hi)&lo!=hi & tab$algorithm=="Davies")
    if(!length(bad))break
    # Recompute the full family because Holm is a step-down procedure.
    fine<-cfg;fine$selection$accuracy<-cfg$selection$accuracy/100^attempt
    fine$selection$max_terms<-cfg$selection$retry_terms
    for(j in which(tab$algorithm=="Davies")) {
      a<-kdr_quadratic_tail(tab$score[j],weights[[j]],fine,kdr_seed(cfg,"selection",job$key,j))
      tab$p_value[j]<-a$p;tab$absolute_accuracy[j]<-a$accuracy;tab$algorithm[j]<-a$method
      tab$ifault[j]<-a$ifault;tab$draws[j]<-a$B
    }
  }
  tab$p_adjusted<-p.adjust(tab$p_value,cfg$selection$adjust)
  tab$selected<-tab$p_adjusted<cfg$selection$alpha
  acc<-ifelse(is.finite(tab$absolute_accuracy),tab$absolute_accuracy,0)
  lo<-p.adjust(pmax(tab$p_value-acc,0),cfg$selection$adjust)<cfg$selection$alpha
  hi<-p.adjust(pmin(tab$p_value+acc,1),cfg$selection$adjust)<cfg$selection$alpha
  ambiguous<-which(!is.na(lo)&!is.na(hi)&lo!=hi & tab$algorithm=="Davies")
  if(length(ambiguous))tab$selected[ambiguous]<-NA
  list(table=tab,valid=all(!is.na(tab$selected)),error=if(all(!is.na(tab$selected)))"" else "Unresolved variable tests")
}
