kdrcc_example <- function(n=800L,seed=13579L,discrete=FALSE) {
  .kdr_scalar(n,"n",20,TRUE)
  with_seed(seed,{
    size<-4L;g<-rep(seq_len(n),each=size)
    W<-if(discrete)sample(c(-1,0,1),length(g),TRUE) else rnorm(length(g))
    shared<-rep(rnorm(n),each=size)
    X1<-.3*shared+.3*W+rnorm(length(g))
    X2<-.2*shared-.2*W+rnorm(length(g))
    X3<-rnorm(length(g));X4<-rnorm(length(g))
    h<-.35*X1-.30*X2+.18*X1*X3
    eta<-h+.25*W+.2*W^2;y<-integer(length(g))
    for(i in seq_len(n)){ix<-which(g==i);pr<-exp(eta[ix]-max(eta[ix]));y[ix[sample.int(size,1,prob=pr)]]<-1L}
    data.frame(set_id=g,case=y,X1=X1,X2=X2,X3=X3,X4=X4,W=W)
  })
}
