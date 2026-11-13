kdrcc_contribution <- function(object,data=object$evaluation,variables=object$specification$exposures,
    reference=NULL,chunk=128L) {
  vars<-object$specification$exposures
  if(!length(variables)||anyDuplicated(variables)||any(!variables%in%vars))stop("variables must name a nonempty unique subset of fitted exposures.")
  .kdr_scalar(chunk,"chunk",1,TRUE)
  A<-.kdr_matrix(object,data);b<-.kdr_reference(object,reference,1L)[1,]
  active<-match(variables,vars);inactive<-setdiff(seq_along(vars),active);s<-object$kernel
  if(length(inactive))A[,inactive]<-matrix(b[inactive],nrow(A),length(inactive),byrow=TRUE)
  factor<-rep(1,nrow(s$centers))
  if(length(inactive))for(j in inactive)factor<-factor*exp(-s$rho[j]*((b[j]-s$center[j])/s$scale[j]-s$centers[,j])^2)
  sa<-s;sa$B<-sweep(s$B,1,factor,"*");sa$centers<-s$centers[,active,drop=FALSE]
  for(n in c("rho","center","scale"))sa[[n]]<-s[[n]][active]
  phi<-sim_shapley_kernel(A[,active,drop=FALSE],b[active],b[active],sa,coef(object),chunk)$total
  colnames(phi)<-variables
  D<-kdr_basis(A,s)-kdr_basis(matrix(b,nrow(A),length(b),byrow=TRUE),s)
  delta<-drop(D%*%coef(object));err<-max(abs(rowSums(phi)-delta))
  if(err>1e-8*max(1,max(abs(delta))))stop("The contribution allocation identity failed.")
  imp<-colMeans(abs(phi));shares<-if(sum(imp)>0)100*imp/sum(imp) else rep(NA_real_,length(imp))
  structure(list(data=data.frame(variable=variables,mean_absolute_contribution=imp,
    mean_signed_contribution=colMeans(phi),share_percent=shares,row.names=NULL,stringsAsFactors=FALSE),
    allocations=phi,log_or=delta,reference=b,variables=variables,identity_error=err,
    target="Baseline Shapley allocation of the specified coordinate contrast; within-set mean absolute shares"),class="kdrcc_contribution")
}
print.kdrcc_contribution <- function(x,...) {cat("KDR-CC baseline contribution shares\n");print(x$data);invisible(x)}
plot.kdrcc_contribution <- function(x,...) {
  d<-x$data;barplot(d$share_percent,names.arg=d$variable,ylab="Mean absolute contribution (%)",col="#0072B2",las=2,...);invisible(x)
}
plot.kdrcc_interaction <- function(x,...) {
  ids<-unique(x$data$pair);old<-par(mfrow=c(ceiling(length(ids)/min(length(ids),2L)),min(length(ids),2L)),mar=c(4,4.5,3,1));on.exit(par(old),add=TRUE)
  for(id in ids){d<-x$data[x$data$pair==id,,drop=FALSE];u<-sort(unique(d$concentration_1));v<-sort(unique(d$concentration_2))
    z<-matrix(d$log_or,nrow=length(u),ncol=length(v));image(u,v,z,xlab=d$exposure_1[1],ylab=d$exposure_2[1],main=id,col=colorRampPalette(c("#2166AC","white","#B2182B"))(31),...)
    contour(u,v,z,add=TRUE,drawlabels=FALSE)}
  invisible(x)
}
