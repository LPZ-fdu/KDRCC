.kdr_matrix <- function(object,data) {
  vars<-object$specification$exposures
  if(is.matrix(data)) {
    if(is.null(colnames(data))) {if(ncol(data)!=length(vars))stop("The exposure matrix has an incorrect number of columns.");colnames(data)<-vars}
    data<-as.data.frame(data)
  }
  .kdr_columns(data,vars);A<-as.matrix(data[vars]);storage.mode(A)<-"double"
  if(any(!is.finite(A)))stop("Prediction exposures must be finite.")
  A
}
.kdr_reference <- function(object,reference,n) {
  p<-length(object$reference);vars<-object$specification$exposures
  if(is.null(reference))reference<-object$reference
  if(is.data.frame(reference)||is.matrix(reference)) {
    R<-.kdr_matrix(object,reference)
    if(nrow(R)==1L)R<-R[rep(1L,n),,drop=FALSE]
    if(nrow(R)!=n)stop("reference must have one row or one row per prediction.")
    return(R)
  }
  if(!is.null(names(reference)))reference<-reference[vars]
  if(length(reference)!=p||any(!is.finite(reference)))stop("reference must specify a finite complete exposure combination.")
  matrix(reference,n,p,byrow=TRUE,dimnames=list(NULL,vars))
}
.kdr_infer <- function(object,D,level=.95,interval=TRUE) {
  est<-drop(D%*%object$coefficients);se<-rep(NA_real_,nrow(D))
  if(interval) {
    if(is.null(object$covariance$matrix))stop("Coefficient covariance is unavailable: ",object$covariance$error)
    vv<-rowSums((D%*%object$covariance$matrix)*D)
    if(any(vv< -1e-9*max(1,max(abs(vv)))))stop("A contrast variance is negative.")
    se<-sqrt(pmax(vv,0))
  }
  if(length(level)!=1L||!is.finite(level)||level<=0||level>=1)stop("level must be between zero and one.")
  crit<-qnorm(1-(1-level)/2);lo<-est-crit*se;hi<-est+crit*se
  data.frame(log_or=est,se=se,lower=lo,upper=hi,or=exp(est),or_lower=exp(lo),or_upper=exp(hi))
}
predict.kdrcc <- function(object,newdata,reference=NULL,type=c("log_or","or"),
    interval=TRUE,level=object$specification$control$confidence,...) {
  type<-match.arg(type);A<-.kdr_matrix(object,newdata);R<-.kdr_reference(object,reference,nrow(A))
  D<-kdr_basis(A,object$kernel)-kdr_basis(R,object$kernel)
  out<-.kdr_infer(object,D,level,interval)
  out$estimate<-out[[type]];out$conf_low<-out[[if(type=="or")"or_lower" else "lower"]]
  out$conf_high<-out[[if(type=="or")"or_upper" else "upper"]]
  attr(out,"reference")<-R;attr(out,"level")<-level;attr(out,"scale")<-type
  out
}
.kdr_effect <- function(data,kind,reference,level,xvars=NULL) {
  structure(list(data=data,kind=kind,reference=reference,level=level,xvars=xvars),class="kdrcc_effect")
}
kdrcc_contrast <- function(object,newdata,reference=NULL,level=.95) {
  A<-.kdr_matrix(object,newdata);R<-.kdr_reference(object,reference,nrow(A))
  .kdr_effect(cbind(as.data.frame(A),predict(object,newdata=A,reference=R,level=level)),"contrast",R,level)
}
kdrcc_joint <- function(object,quantiles=seq(.1,.9,length.out=21L),variables=object$specification$exposures,
    grid=NULL,reference=NULL,level=.95) {
  if(any(!variables%in%object$specification$exposures)||anyDuplicated(variables))stop("variables must name unique fitted exposures.")
  if(is.null(grid)) {
    if(any(!is.finite(quantiles))||any(quantiles<0|quantiles>1))stop("quantiles must lie between zero and one.")
    A<-.kdr_reference(object,reference,length(quantiles));E<-.kdr_matrix(object,object$evaluation)
    for(v in variables)A[,v]<-as.numeric(quantile(E[,v],quantiles))
    x<-quantiles
  } else {A<-.kdr_matrix(object,grid);x<-seq_len(nrow(A))}
  out<-cbind(data.frame(position=x),as.data.frame(A),predict(object,A,reference=reference,level=level))
  .kdr_effect(out,"joint",.kdr_reference(object,reference,1L),level,"position")
}
kdrcc_single <- function(object,exposure=object$specification$exposures,values=NULL,
    background=NULL,points=41L,level=.95) {
  if(any(!exposure%in%object$specification$exposures)||anyDuplicated(exposure))stop("exposure must name unique fitted exposures.")
  .kdr_scalar(points,"points",3,TRUE);b<-.kdr_reference(object,background,1L)
  E<-.kdr_matrix(object,object$evaluation);out<-list()
  for(v in exposure) {
    x<-if(is.null(values))seq(quantile(E[,v],.1),quantile(E[,v],.9),length.out=points) else if(is.list(values))values[[v]] else values
    if(!length(x)||any(!is.finite(x)))stop("values must provide finite exposure concentrations.")
    A<-b[rep(1L,length(x)),,drop=FALSE];A[,v]<-x
    out[[v]]<-cbind(data.frame(exposure=v,concentration=x,stringsAsFactors=FALSE),predict(object,A,reference=b,level=level))
  }
  .kdr_effect(kdr_bind(out),"single",b,level,c("exposure","concentration"))
}
coef.kdrcc <- function(object,...) object$coefficients
vcov.kdrcc <- function(object,...) object$covariance$matrix
summary.kdrcc <- function(object,...) {
  structure(list(call=object$call,data=object$data_summary,
    parameters=c(M=object$kernel$dimension,kappa=object$kernel$kappa,lambda=object$specification$control$lambda),
    candidate_centers=nrow(object$kernel$centers),
    covariance=object$covariance$status,gradient=object$diagnostics$gradient,
    objective=object$diagnostics$objective,exposures=object$specification$exposures,
    inference_target=object$covariance$target),class="summary.kdrcc")
}
print.kdrcc <- function(x,...) {print(summary(x));invisible(x)}
print.summary.kdrcc <- function(x,...) {
  cat("Kernel-Based Doubly Robust Case-Crossover Regression\n")
  cat("Estimation matched sets:",x$data$N,"| Independent learning sets:",x$data$L,"\n")
  cat("Exposures:",paste(x$exposures,collapse=", "),"\n")
  cat("Configuration:",paste(names(x$parameters),format(x$parameters),sep="=",collapse=", "),"\n")
  if(!is.null(x$candidate_centers))
    cat("Kernel basis:",x$parameters[["M"]],"retained of",x$candidate_centers,"candidate centers\n")
  cat("Objective:",format(x$objective),"| Gradient:",format(x$gradient),"\n")
  cat("Covariance:",x$covariance,"\n");invisible(x)
}
print.kdrcc_effect <- function(x,...) {cat("KDR-CC",x$kind,"effects\n");print(utils::head(x$data,10L));invisible(x)}
kdrcc_diagnostics <- function(object) {
  list(data=object$data_summary,estimation=object$diagnostics$estimation,
    kernel_learning=object$diagnostics$kernel,
    kernel_basis=attr(object$kernel$B,"kernel_basis",exact=TRUE),
    covariance_status=object$covariance$status,covariance_error=object$covariance$error,
    distance_parameters=setNames(object$kernel$rho,object$specification$exposures),
    boundaries=c(lower=object$diagnostics$rho_lower,upper=object$diagnostics$rho_upper),
    minimum_curvature=object$covariance$minimum_curvature,curvature_condition=object$covariance$condition,
    auxiliary=lapply(object$nuisance$models,function(a)list(g_iterations=a$g$iterations,
      f_minimum_cell_mass=a$f$minimum_mass)),inference_target=object$covariance$target,
    elapsed_seconds=object$diagnostics$elapsed_seconds)
}
plot.kdrcc_effect <- function(x,scale=c("or","log_or"),...) {
  scale<-match.arg(scale);d<-x$data
  groups<-if(x$kind=="single")unique(d$exposure) else "Joint exposure"
  old<-par(mfrow=c(ceiling(length(groups)/min(3,length(groups))),min(3,length(groups))),mar=c(4.3,4.5,2.6,1))
  on.exit(par(old),add=TRUE)
  for(g in groups) {
    z<-if(x$kind=="single")d[d$exposure==g,,drop=FALSE] else d
    xx<-if(x$kind=="single")z$concentration else if("position"%in%names(z))z$position else seq_len(nrow(z))
    yy<-z[[scale]];lo<-z[[if(scale=="or")"or_lower" else "lower"]];hi<-z[[if(scale=="or")"or_upper" else "upper"]]
    plot(xx,yy,type="n",ylim=range(c(yy,lo,hi),finite=TRUE),xlab=if(x$kind=="single")g else "Exposure path position",
      ylab=if(scale=="or")"Odds ratio" else "Log odds ratio",main=g,...)
    if(all(is.finite(c(lo,hi))))polygon(c(xx,rev(xx)),c(lo,rev(hi)),col=adjustcolor("#0072B2",alpha.f=.18),border=NA)
    abline(h=if(scale=="or")1 else 0,lty=3,col="grey50");lines(xx,yy,col="#0072B2",lwd=2)
  }
  invisible(x)
}
plot.kdrcc <- function(x,...) {
  d<-x$diagnostics$estimation;old<-par(mfrow=c(1,2),mar=c(4,4.5,2.5,1));on.exit(par(old),add=TRUE)
  plot(d$iteration,d$objective,type="b",xlab="Iteration",ylab="Regularized criterion",main="Estimation",...)
  plot(d$iteration,pmax(d$gradient,.Machine$double.eps),type="b",log="y",xlab="Iteration",ylab="Maximum gradient",main="Numerical convergence",...)
  invisible(x)
}
