# Baseline Shapley allocation for separable Gaussian kernels.
sim_shapley_gauss <- function(n) {
  if(n==1L)return(list(nodes=.5,weights=1))
  J<-matrix(0,n,n);i<-seq_len(n-1L);v<-i/sqrt(4*i^2-1)
  J[cbind(i,i+1L)]<-v;J[cbind(i+1L,i)]<-v;e<-eigen(J,symmetric=TRUE)
  o<-order(e$values);list(nodes=(e$values[o]+1)/2,weights=e$vectors[1,o]^2)
}
sim_shapley_enumerate <- function(A,reference,predict) {
  A<-as.matrix(A);p<-ncol(A);masks<-0:(2^p-1L)
  membership<-sapply(seq_len(p),function(j)bitwAnd(masks,bitwShiftL(1L,j-1L))!=0L)
  V<-matrix(0,nrow(A),length(masks));W<-matrix(0,length(masks),p)
  for(i in seq_along(masks)) {
    X<-A;off<-which(!membership[i,]);if(length(off))X[,off]<-matrix(reference[off],nrow(A),length(off),byrow=TRUE)
    V[,i]<-predict(X)
    s<-sum(membership[i,]);if(s<p)for(j in which(!membership[i,])) {
      w<-1/(p*choose(p-1L,s));plus<-i+2^(j-1L)
      W[plus,j]<-W[plus,j]+w;W[i,j]<-W[i,j]-w
    }
  }
  list(total=V%*%W,delta=V[,ncol(V)]-V[,1])
}
sim_shapley_main <- function(A,reference,anchor,predict) {
  p<-ncol(A);out<-matrix(0,nrow(A),p)
  for(j in seq_len(p)) {
    X<-matrix(anchor,nrow(A),p,byrow=TRUE);Y<-matrix(anchor,1,p)
    X[,j]<-A[,j];Y[,j]<-reference[j]
    out[,j]<-predict(X)-predict(Y)[1]
  }
  out
}
sim_shapley_kernel <- function(A,reference,anchor,s,beta,chunk=128L) {
  stopifnot(s$kind=="kernel",ncol(A)==length(s$rho),length(reference)==ncol(A),length(anchor)==ncol(A))
  p<-ncol(A);m<-nrow(s$centers);alpha<-drop(s$B%*%beta)
  if(length(alpha)!=m)stop("Kernel coefficient dimensions do not agree.")
  factor_at<-function(x,j)exp(-s$rho[j]*outer((x-s$center[j])/s$scale[j],s$centers[,j],"-")^2)
  f0<-lapply(seq_len(p),function(j)drop(factor_at(reference[j],j)))
  fc<-lapply(seq_len(p),function(j)drop(factor_at(anchor[j],j)))
  co<-lapply(seq_len(p),function(j)Reduce(`*`,fc[-j]))
  quadrature<-sim_shapley_gauss(ceiling(p/2))
  total<-main<-matrix(0,nrow(A),p)
  for(lo in seq.int(1L,nrow(A),by=as.integer(chunk))) {
    ids<-lo:min(nrow(A),lo+chunk-1L);n<-length(ids)
    f1<-lapply(seq_len(p),function(j)factor_at(A[ids,j],j))
    d<-lapply(seq_len(p),function(j)sweep(f1[[j]],2,f0[[j]]))
    for(j in seq_len(p))main[ids,j]<-drop(d[[j]]%*%(co[[j]]*alpha))
    for(q in seq_along(quadrature$nodes)) {
      t<-quadrature$nodes[q];z<-lapply(seq_len(p),function(j)sweep(t*f1[[j]],2,(1-t)*f0[[j]],"+"))
      pre<-vector("list",p);pre[[1]]<-matrix(1,n,m)
      if(p>1L)for(j in 2:p)pre[[j]]<-pre[[j-1L]]*z[[j-1L]]
      suf<-matrix(1,n,m)
      for(j in p:1L) {
        total[ids,j]<-total[ids,j]+quadrature$weights[q]*drop((d[[j]]*pre[[j]]*suf)%*%alpha)
        suf<-suf*z[[j]]
      }
    }
  }
  list(total=total,main=main,interaction=total-main)
}
