run_module_D_checks <- function(source_root,output_root=tempdir()) {
  source(file.path(source_root,"00_parameters.R"),local=TRUE,encoding="UTF-8")
  cfg<-simulation_config(source_root);checks<-list()
  add<-function(name,value){checks[[length(checks)+1L]]<<-data.frame(check=name,passed=isTRUE(value));if(!isTRUE(value))stop(name)}
  resolved<-list(parameters=list(N=20000L,L=20000L,M=144L,kappa=1e-4,lambda=1),configuration_id="test")
  pp<-sim_plan(cfg,"D",resolved);p<-pp$tasks
  add("D has 1600 new fits",nrow(p)==1600L)
  add("D has 200 new learning tasks",length(unique(p$kernel_key))==200L)
  add("D has eight fits per scenario/repetition",all(table(p$correlation_id,p$replicate)==8L))
  add("D excludes reused baseline from fitting",!any(p$correlation_id=="D0"))
  add("D unique keys",!anyDuplicated(p$key))
  add("D separate caches",length(unique(vapply(sim_D_data_jobs(p),sim_data_key,character(1))))==200L)
  bad<-cfg;bad$D$correlation_levels$min_abs[1]<-.9
  add("Invalid floor rejected",inherits(try(sim_D_levels(bad),silent=TRUE),"try-error"))
  cf<-cfg;cf$batch<-"correlation";cb<-cfg;cb$batch<-"validation"
  add("Independent D namespace",kdr_seed(cf,"data_pool",1L)!=kdr_seed(cb,"data_pool",1L))
  n<-12000L;prob<-(seq_len(n)-.5)/n
  maps<-list(values=rep(list(qnorm(prob)),8),probability=prob)
  A<-with_seed(9871L,matrix(rnorm(n*8L),n,8L));Z<-sim_D_gaussian_scores(A,maps)
  group<-rep(seq_len(n/4),each=4);date<-as.Date("2021-01-01")+rep(c(0,7,14,21),n/4)
  U<-sim_D_common_score(group,date,.6,112L)
  add("Common score reproducible",identical(U,sim_D_common_score(group,date,.6,112L)))
  for(eta in c(.6,.85)) {
    X<-sim_D_transform(Z,U,maps,eta,cfg$baseline$dgp$direction)
    R<-cor(X);v<-abs(R[upper.tri(R)])
    add(paste("Target correlation",eta),abs(mean(v)-eta)<.02)
    add(paste("Minimum correlation",eta),min(v)>eta-.04)
    add(paste("Positive-definite correlation",eta),min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)>0)
    add(paste("Marginal variance preservation",eta),max(abs(apply(X,2,var)-1))<.06)
    add(paste("Signed dependence",eta),all(sign(R[upper.tri(R)])==sign(tcrossprod(cfg$baseline$dgp$direction))[upper.tri(R)]))
  }
  add("Identity without recoupling",identical(sim_D_apply_copula(A,group,date,cfg,1L),A))
  old<-kdr_generate(4L,8L,cfg,154L)
  cf<-cfg;cf$baseline$dgp$copula_eta<-.6;cf$D$copula<-list(maps=maps)
  altered<-kdr_generate(4L,8L,cf,154L)
  add("Matching and confounder retained",identical(old[c("W","group","date","fold")],altered[c("W","group","date","fold")]))
  add("One event per matched set",all(tapply(altered$event,altered$group,sum)==1))
  dir.create(file.path(output_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  out<-kdr_bind(checks);kdr_csv(out,file.path(output_root,"validation","module_D_checks.csv"))
  kdr_log("Module D checks passed:",nrow(out));invisible(out)
}
