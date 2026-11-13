.kdr_tail <- function(score,weights,control,family=1L,seed=1L,accuracy=control$tail_accuracy,alpha=.05) {
  w<-weights[weights>0]
  if(!length(w))return(list(p=if(score<=0)1 else 0,accuracy=0,method="Point mass",status=0L,draws=0L))
  if(max(w)-min(w)<=1e-13*max(w))return(list(p=pchisq(score/mean(w),length(w),lower.tail=FALSE),
    accuracy=10*.Machine$double.eps,method="Scaled chi-square",status=0L,draws=0L))
  norm<-sum(w);q<-score/norm;w<-w/norm
  if(q<=0)return(list(p=1,accuracy=0,method="Zero statistic",status=0L,draws=0L))
  d<-suppressWarnings(CompQuadForm::davies(q,w,lim=control$tail_terms,acc=accuracy))
  if(d$ifault!=0L||!is.finite(d$Qq)||d$Qq<0||d$Qq>1)
    d<-suppressWarnings(CompQuadForm::davies(q,w,lim=control$tail_retry_terms,acc=accuracy))
  if(d$ifault==0L&&is.finite(d$Qq)&&d$Qq>=0&&d$Qq<=1)return(list(p=d$Qq,accuracy=accuracy,method="Davies",status=d$ifault,draws=0L))
  B<-max(control$tail_draws,ceiling(200*family/alpha)-1L)
  count<-with_seed(seed,{
    done<-0;hit<-0
    while(done<B) {n<-min(10000L,B-done);v<-drop(matrix(rnorm(n*length(w))^2,nrow=n)%*%w);hit<-hit+sum(v>=q);done<-done+n}
    hit
  })
  p<-(count+1)/(B+1)
  list(p=p,accuracy=sqrt(p*(1-p)/B),method="Monte Carlo",status=d$ifault,draws=B)
}
.kdr_quadratic_family <- function(scores,weights,control,alpha,adjust) {
  n<-length(scores);a<-lapply(seq_len(n),function(j).kdr_tail(scores[j],weights[[j]],control,n,control$seed+10L+j,alpha=alpha))
  make<-function()data.frame(p_value=vapply(a,`[[`,numeric(1),"p"),accuracy=vapply(a,`[[`,numeric(1),"accuracy"),
    algorithm=vapply(a,`[[`,character(1),"method"),draws=vapply(a,`[[`,numeric(1),"draws"),stringsAsFactors=FALSE)
  z<-make()
  for(it in 1:2) {
    lo<-p.adjust(pmax(0,z$p_value-z$accuracy),adjust)<alpha
    hi<-p.adjust(pmin(1,z$p_value+z$accuracy),adjust)<alpha
    if(!any(lo!=hi))break
    for(j in which(z$algorithm=="Davies"))a[[j]]<-.kdr_tail(scores[j],weights[[j]],control,n,control$seed+10L+j,control$tail_accuracy/100^it,alpha=alpha)
    z<-make()
  }
  z$p_adjusted<-p.adjust(z$p_value,adjust);z$selected<-z$p_adjusted<alpha
  lo<-p.adjust(pmax(0,z$p_value-z$accuracy),adjust)<alpha
  hi<-p.adjust(pmin(1,z$p_value+z$accuracy),adjust)<alpha
  z$numerically_stable<-lo==hi;z$selected[!z$numerically_stable]<-NA
  z
}
kdrcc_variables <- function(object,data=object$evaluation,alpha=.05,adjust="holm") {
  if(alpha<=0||alpha>=1)stop("alpha must lie between zero and one.")
  A<-.kdr_matrix(object,data);V<-vcov(object)
  if(is.null(V))stop("Coefficient covariance is required for variable tests.")
  L<-kdr_vroot(V);p<-ncol(A);scores<-numeric(p);weights<-vector("list",p)
  for(j in seq_len(p)) {
    D<-kdr_derivative(A,object$kernel,j);Q<-crossprod(D)/nrow(D)
    scores[j]<-max(0,drop(crossprod(coef(object),Q%*%coef(object))))
    weights[[j]]<-pmax(eigen(kdr_sym(crossprod(L,Q%*%L)),symmetric=TRUE,only.values=TRUE)$values,0)
  }
  z<-cbind(data.frame(variable=colnames(A),statistic=scores,stringsAsFactors=FALSE),
    .kdr_quadratic_family(scores,weights,object$specification$control,alpha,adjust))
  structure(list(data=z,alpha=alpha,adjust=adjust,evaluation_points=nrow(A),
    target="Zero coordinate derivative of the regularized working response on the evaluation grid"),class="kdrcc_test")
}
kdrcc_interaction <- function(object,pairs,values=NULL,background=NULL,points=9L,alpha=.05,adjust="holm",level=.95) {
  pairs<-as.matrix(pairs)
  if(ncol(pairs)!=2L||nrow(pairs)<1L||any(!pairs%in%object$specification$exposures)||any(pairs[,1]==pairs[,2]))stop("pairs must be a two-column matrix of distinct fitted exposure names.")
  .kdr_scalar(points,"points",3,TRUE);if(alpha<=0||alpha>=1)stop("alpha must lie between zero and one.")
  b<-.kdr_reference(object,background,1L);E<-.kdr_matrix(object,object$evaluation)
  V<-vcov(object);if(is.null(V))stop("Coefficient covariance is required for interaction tests.")
  L<-kdr_vroot(V);output<-list();weights<-list();scores<-numeric(nrow(pairs))
  for(i in seq_len(nrow(pairs))) {
    j<-pairs[i,1];k<-pairs[i,2]
    getgrid<-function(v)sort(unique(c(if(is.null(values))seq(quantile(E[,v],.1),quantile(E[,v],.9),length.out=points) else values[[v]],b[1,v])))
    xy<-expand.grid(first=getgrid(j),second=getgrid(k));n<-nrow(xy)
    D0<-b[rep(1L,n),,drop=FALSE];A<-B<-C<-D0
    A[,j]<-B[,j]<-xy$first;A[,k]<-C[,k]<-xy$second
    D<-kdr_basis(A,object$kernel)-kdr_basis(B,object$kernel)-kdr_basis(C,object$kernel)+kdr_basis(D0,object$kernel)
    anchor<-abs(xy$first-b[1,j])<1e-12|abs(xy$second-b[1,k])<1e-12
    D[anchor,]<-0;pred<-.kdr_infer(object,D,level,TRUE)
    good<-!anchor & pred$se>0;pv<-rep(NA_real_,n);pv[good]<-2*pnorm(abs(pred$log_or[good]/pred$se[good]),lower.tail=FALSE)
    output[[i]]<-cbind(data.frame(pair=paste(j,k,sep="__"),exposure_1=j,exposure_2=k,
      concentration_1=xy$first,concentration_2=xy$second,eligible=!anchor,stringsAsFactors=FALSE),pred,p_value=pv)
    if(!any(!anchor))stop("Interaction grids require non-anchor concentrations.")
    ix<-which(!anchor);Z<-D[ix,,drop=FALSE]%*%L/sqrt(length(ix));gram<-if(nrow(Z)<=ncol(Z))tcrossprod(Z) else crossprod(Z)
    weights[[i]]<-pmax(eigen(kdr_sym(gram),symmetric=TRUE,only.values=TRUE)$values,0)
    scores[i]<-mean(pred$log_or[ix]^2)
  }
  all<-kdr_bind(output);all$p_adjusted<-NA_real_;ix<-which(all$eligible)
  all$p_adjusted[ix]<-p.adjust(all$p_value[ix],adjust,n=length(ix))
  all$significant<-all$p_adjusted<alpha
  all$direction<-ifelse(all$log_or>0,"positive",ifelse(all$log_or<0,"negative","zero"))
  overall<-cbind(data.frame(pair=vapply(output,function(x)x$pair[1],character(1)),statistic=scores,stringsAsFactors=FALSE),
    .kdr_quadratic_family(scores,weights,object$specification$control,alpha,adjust))
  structure(list(overall=overall,data=all,reference=b,alpha=alpha,level=level,adjust=adjust,
    target="Zero four-point interaction contrast of the regularized working response"),class="kdrcc_interaction")
}
print.kdrcc_test <- function(x,...) {cat("KDR-CC variable participation tests\n");print(x$data);invisible(x)}
print.kdrcc_interaction <- function(x,...) {cat("KDR-CC interaction tests\n");print(x$overall);invisible(x)}
