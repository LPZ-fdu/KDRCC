# Fixed-background response evaluation from saved module B fits.
sim_B_response_settings <- function(report,cfg) {
  defaults<-list(B_response_summaries=TRUE,B_single_points=41L,
    B_interaction_points=9L,B_interaction_pairs=NULL,B_repetition_ranges=TRUE)
  s<-sim_merge(defaults,report)
  for(n in c("B_response_summaries","B_repetition_ranges"))
    if(length(s[[n]])!=1L||is.na(s[[n]])||!is.logical(s[[n]]))stop(paste(n,"must be TRUE or FALSE."))
  for(n in c("B_single_points","B_interaction_points"))
    if(length(s[[n]])!=1L||!is.numeric(s[[n]])||!is.finite(s[[n]])||s[[n]]<3||s[[n]]!=floor(s[[n]]))
      stop(paste(n,"must be an integer of at least three."))
  pairs<-s$B_interaction_pairs
  if(is.null(pairs))pairs<-as.matrix(cfg$baseline$truth$interactions[c("q","r")])
  pairs<-as.matrix(pairs)
  if(ncol(pairs)!=2L)stop("B_interaction_pairs must have two columns of exposure names or indices.")
  if(is.character(pairs))pairs<-matrix(match(pairs,exposure_names),ncol=2L)
  if(any(!is.finite(pairs))||any(pairs!=floor(pairs))||any(pairs<1|pairs>8)||any(pairs[,1]==pairs[,2]))
    stop("Interaction pairs must contain two distinct active exposures from the eight-exposure model.")
  if(nrow(pairs))pairs<-t(apply(pairs,1,sort))
  pairs<-unique(matrix(as.integer(pairs),ncol=2L))
  colnames(pairs)<-c("first","second");s$B_interaction_pairs<-pairs;s
}

sim_B_response_grid <- function(ev,cfg,settings) {
  p<-length(ev$reference)
  if(p<8L)stop("Response evaluation requires the eight active exposures.")
  if(ncol(ev$joint)!=p||any(!is.finite(ev$reference))||any(!is.finite(ev$joint)))
    stop("Invalid saved evaluation coordinates.")
  grid<-function(j,n)sort(unique(c(seq(min(ev$joint[,j]),max(ev$joint[,j]),length.out=n),ev$reference[j])))
  at_reference<-function(x,j)abs(x-ev$reference[j])<=1e-12*max(1,abs(ev$reference[j]))
  matrices<-rep(list(list()),4L);points<-list();block<-0L
  add<-function(A,B,C,D,target,component,j,k=NA_integer_) {
    block<<-block+1L;n<-nrow(A)
    zero<-at_reference(A[,j],j)
    if(!is.na(k))zero<-zero|at_reference(A[,k],k)
    truth<-kdr_truth(A,cfg)-kdr_truth(B,cfg)
    if(target=="interaction")truth<-truth-kdr_truth(C,cfg)+kdr_truth(D,cfg)
    truth[zero]<-0
    points[[block]]<<-data.frame(target=target,component=component,
      exposure_1=exposure_names[j],exposure_2=if(is.na(k))"" else exposure_names[k],
      concentration_1=A[,j],concentration_2=if(is.na(k))NA_real_ else A[,k],
      truth=truth,eligible=!zero,stringsAsFactors=FALSE)
    for(i in 1:4)matrices[[i]][[block]]<<-list(A,B,C,D)[[i]]
  }
  for(j in seq_len(8L)) {
    x<-grid(j,settings$B_single_points);ref<-matrix(ev$reference,length(x),p,byrow=TRUE)
    A<-ref;A[,j]<-x
    add(A,ref,ref,ref,"single",exposure_names[j],j)
  }
  pairs<-settings$B_interaction_pairs
  for(i in seq_len(nrow(pairs))) {
    j<-pairs[i,1];k<-pairs[i,2]
    xy<-expand.grid(a=grid(j,settings$B_interaction_points),b=grid(k,settings$B_interaction_points))
    D<-matrix(ev$reference,nrow(xy),p,byrow=TRUE);A<-B<-C<-D
    A[,j]<-B[,j]<-xy$a;A[,k]<-C[,k]<-xy$b
    add(A,B,C,D,"interaction",paste(exposure_names[c(j,k)],collapse="__"),j,k)
  }
  tab<-kdr_bind(points);tab$index<-seq_len(nrow(tab))
  matrices<-lapply(matrices,function(a)do.call(rbind,a));names(matrices)<-c("both","first","second","reference")
  list(points=tab,coordinates=matrices,reference=ev$reference,
    concentration_ranges=data.frame(exposure=kdr_names(p),
      lower=apply(ev$joint,2,min),reference=ev$reference,upper=apply(ev$joint,2,max)),
    settings=settings)
}

sim_B_response_operator <- function(grid,s) {
  x<-grid$coordinates
  D<-kdr_basis(x$both,s)-kdr_basis(x$first,s)
  ix<-grid$points$target=="interaction"
  if(any(ix))D[ix,]<-D[ix,,drop=FALSE]-kdr_basis(x$second[ix,,drop=FALSE],s)+kdr_basis(x$reference[ix,,drop=FALSE],s)
  D[!grid$points$eligible,]<-0
  D
}

