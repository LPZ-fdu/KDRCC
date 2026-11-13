# Algebraic checks of the prespecified active-exposure allocation.
run_contribution_checks <- function(source_root,output_root=tempdir()) {
  e<-new.env(parent=.GlobalEnv);sys.source(file.path(source_root,"00_parameters.R"),e)
  cfg<-e$simulation_config(source_root);rows<-list()
  check<-function(name,ok){rows[[length(rows)+1L]]<<-data.frame(check=name,passed=isTRUE(ok));if(!isTRUE(ok))stop(name)}
  fixture<-with_seed(812731L,{
    b<-c(cfg$baseline$reference$center,.3,-.2);sc<-c(cfg$baseline$reference$scale,1,1)
    A<-matrix(b,5,10,byrow=TRUE)+sweep(matrix(rnorm(50),5),2,sc*.3,"*")
    s<-list(kind="kernel",center=b-.1*sc,scale=sc,centers=matrix(rnorm(50),5,10),
      rho=rep(.13,10),B=matrix(rnorm(25),5),M=5L,dimension=5L)
    list(ev=list(reference=b,individual=A,interpretation=list(ids=1:5)),s=s,fit=list(beta=rnorm(5)))
  })
  ev<-fixture$ev;s<-fixture$s;fit<-fixture$fit;A<-ev$individual;b<-ev$reference
  z<-sim_active_contributions(fit,s,ev,cfg)
  fixed<-A;fixed[,9:10]<-matrix(b[9:10],nrow(A),2,byrow=TRUE)
  direct<-sim_shapley_kernel(fixed,b,b,s,fit$beta)$total
  exact<-sim_shapley_enumerate(A[,1:8],b[1:8],function(X){
    full<-matrix(b,nrow(X),10,byrow=TRUE);full[,1:8]<-X;drop(kdr_basis(full,s)%*%fit$beta)
  })$total
  check("Active Gaussian restriction agrees with subset enumeration",max(abs(direct[,1:8]-exact))<1e-10)
  check("Inactive fixed-coordinate allocations are zero",max(abs(direct[,9:10]))<1e-12)
  v<-colMeans(abs(exact));check("Active shares use only eight coordinates",
    max(abs(z$contributions$estimated_share-100*v/sum(v)))<1e-10)
  check("Contribution identity",z$contribution_identity_error<1e-10)
  check("Shares sum to 100",abs(sum(z$contributions$estimated_share)-100)<1e-10)
  check("Eight prespecified exposure rows",identical(z$contributions$variable,exposure_names)&&all(z$contributions$true_active))
  ev2<-ev;ev2$individual[,9:10]<-ev2$individual[,9:10]+100
  z2<-sim_active_contributions(fit,s,ev2,cfg)
  check("Evaluation invariant to inactive observed coordinates",identical(z$contributions,z2$contributions))
  s8<-s;s8$centers<-s$centers[,1:8];for(n in c("center","scale","rho"))s8[[n]]<-s[[n]][1:8]
  ev8<-ev;ev8$reference<-b[1:8];ev8$individual<-A[,1:8]
  z8<-sim_active_contributions(fit,s8,ev8,cfg)
  all8<-colMeans(abs(sim_shapley_kernel(A[,1:8],b[1:8],b[1:8],s8,fit$beta)$total))
  check("Eight-dimensional allocation unchanged",max(abs(z8$contributions$estimated_share-100*all8/sum(all8)))<1e-10)
  sl<-list(kind="linear",center=b,scale=rep(1,10),B=diag(10))
  fl<-list(beta=c(rep(.01,8),20,-10));zl<-sim_active_contributions(fl,sl,ev,cfg)
  vl<-colMeans(abs(sweep(A[,1:8],2,b[1:8],"-")*.01))
  check("Linear comparator active allocation",max(abs(zl$contributions$estimated_share-100*vl/sum(vl)))<1e-10)
  check("Truth shares sum to 100",abs(sum(z$contributions$true_share)-100)<1e-10)
  dir.create(file.path(output_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  out<-kdr_bind(rows);kdr_csv(out,file.path(output_root,"validation","contribution_checks.csv"))
  message("Active contribution checks passed: ",nrow(out));invisible(out)
}
