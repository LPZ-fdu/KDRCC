# Algebraic verification of the scientific specification.
run_method_contract_checks <- function(root, output_root=tempdir()) {
  source(file.path(root,"00_parameters.R"),encoding="UTF-8")
  cfg<-simulation_config(root);cfg$batch<-"contract";checks<-list()
  check<-function(name,value){stopifnot(isTRUE(value));checks[[name]]<<-TRUE}
  Z<-c(-1,.2,1.4,-.7,.4,1.9);X<-rep(1:2,each=3);mass<-c(.1,.2,.2,.15,.25,.1)
  b<-.4;g0<-c(-.3,.7);pr<-plogis(b*Z+g0[X])
  f0<-vapply(1:2,function(c)sum(mass[X==c]*(1-pr[X==c])*Z[X==c])/sum(mass[X==c]*(1-pr[X==c])),numeric(1))
  expectation<-function(g,f)sum(mass*(Z-f[X])*(pr*exp(-b*Z-g[X])-(1-pr)))
  check("DR: correct g, misspecified f",abs(expectation(g0,c(2,-1)))<1e-12)
  check("DR: misspecified g, correct f",abs(expectation(c(.6,-.4),f0))<1e-12)
  g<-c(.6,-.4);f<-c(2,-1)
  factorized<-sum(vapply(1:2,function(c)sum(mass[X==c]*(1-pr[X==c]))*(exp(g0[c]-g[c])-1)*(f0[c]-f[c]),numeric(1)))
  check("Exact product-bias identity",abs(expectation(g,f)-factorized)<1e-12)
  cfg$M<-6L;d<-kdr_generate(120L,8L,cfg,9023L);s<-kdr_geometry(d,cfg,901L)
  s$B<-diag(6);s$T<-diag(6);s$penalty<-diag(6)
  z<-d;z$X<-kdr_basis(d$A,s);z$pairs<-kdr_pairs(d)
  z$pairs$Z<-z$X[z$pairs$plus,,drop=FALSE]-z$X[z$pairs$minus,,drop=FALSE]
  each<-function(FUN,progress=FALSE)FUN(z)
  aux<-list(g=lapply(1:3,function(k)list(gamma=c(.35,.8))),f=kdr_f_fit(each,cfg,d$n))
  fit<-list(beta=seq_len(6)*.01,lambda=.2);a<-kdr_curvature(fit,each,aux,s,d$n);eps<-1e-5
  grad<-function(b,ax=aux,it=each)kdr_objective(b,it,ax,s,fit$lambda,d$n)$gradient
  numH<-vapply(1:6,function(j){bp<-bm<-fit$beta;bp[j]<-bp[j]+eps;bm[j]<-bm[j]-eps;(grad(bp)-grad(bm))/(2*eps)},numeric(6))
  check("Exact regularized Hessian",max(abs(a$H-numH))<1e-7)
  for(k in 1:3) {
    ng<-vapply(1:2,function(j){ap<-am<-aux;ap$g[[k]]$gamma[j]<-ap$g[[k]]$gamma[j]+eps
      am$g[[k]]$gamma[j]<-am$g[[k]]$gamma[j]-eps;(grad(fit$beta,ap)-grad(fit$beta,am))/(2*eps)},numeric(6))
    check(paste("Nuisance g score sensitivity, fold",k),max(abs(a$Ag[[k]]-ng))<1e-7)
    for(c in c(1L,5L,9L)) {
      nf<-vapply(1:6,function(j){ap<-am<-aux;ap$f$means[[k]][c,j]<-ap$f$means[[k]][c,j]+eps
        am$f$means[[k]][c,j]<-am$f$means[[k]][c,j]-eps;(grad(fit$beta,ap)-grad(fit$beta,am))/(2*eps)},numeric(6))
      check(paste("Nuisance f score sensitivity, fold/cell",k,c),max(abs(a$Af[[k]][[c]]-nf))<1e-7)
    }
  }
  job<-list(N=20000L,L=20000L,M=144L,kappa=1e-4,lambda=1,batch="validation")
  cf<-sim_job_config(cfg,job)
  check("Configured lambda is the effective penalty",cf$penalty$a*cf$penalty$reference_N/cf$N==job$lambda)
  check("Fixed positive penalty has an explicit target shift",sum((.2*diag(6)%*%fit$beta)^2)>0)
  s2<-s;s2$penalty<-kdr_kernel_penalty(d,s,cfg,12L)
  check("Mixed-derivative penalty is PSD",min(eigen(s2$penalty,symmetric=TRUE,only.values=TRUE)$values)>-1e-10)
  check("Penalty trace normalization",abs(mean(diag(s2$penalty))-1)<1e-12)
  out<-data.frame(check=names(checks),passed=unlist(checks),row.names=NULL)
  dir.create(file.path(output_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  kdr_csv(out,file.path(output_root,"validation","method_contract_checks.csv"))
  kdr_log("Method specification checks passed:",nrow(out));invisible(out)
}
