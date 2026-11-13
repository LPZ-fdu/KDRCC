run_module_E_checks <- function(root,output_root=tempdir()) {
  source(file.path(root,"00_parameters.R"),encoding="UTF-8",local=TRUE)
  cfg<-simulation_config(root);cfg$batch<-"E_checks"
  resolved<-list(parameters=list(N=20000L,L=20000L,M=144L,kappa=1e-4,lambda=1))
  pp<-sim_plan_E(cfg,resolved);d<-pp$design$configurations
  stopifnot(nrow(d)==18L,nrow(pp$tasks)==1400L,nrow(pp$design$scenarios)==7L,
    pp$cv_fits_per_repetition==54L,sum(d$is_fixed)==1L,!anyDuplicated(pp$tasks$key))
  bad<-cfg;bad$E$null_counts<-2.5
  stopifnot(inherits(try(sim_E_scenarios(bad),silent=TRUE),"try-error"))
  bad<-cfg;bad$E$manual_candidates<-data.frame(M=10.5,kappa=1e-4,lambda=1)
  stopifnot(inherits(try(sim_E_candidates(bad,resolved),silent=TRUE),"try-error"))
  x<-expand.grid(candidate_id=d$candidate_id,fold=1:3,stringsAsFactors=FALSE)
  x$ok<-TRUE;x$n<-100L;x$loss_sum<-match(x$candidate_id,d$candidate_id)*10
  x$loss_sum[x$candidate_id==d$candidate_id[1]&x$fold==2]<-NA_real_
  choice<-sim_E_choose(x,d,3L)
  stopifnot(choice$candidate$candidate_id==d$candidate_id[2],!choice$fallback)
  x$ok[]<-FALSE;fallback<-sim_E_choose(x,d,3L)
  stopifnot(fallback$fallback,identical(fallback$candidate$candidate_id,d$candidate_id[d$is_fixed]))
  dat<-kdr_generate(60L,10L,cfg,903L);dat$pairs<-kdr_pairs(dat)
  groups<-seq(1L,dat$n,by=2L);part<-sim_E_take_groups(dat,groups)
  stopifnot(part$n==length(groups),all(drop(group_sum(part$event,part$group))==1),
    identical(part$A,dat$A[dat$group%in%groups,,drop=FALSE]),is.null(part$truth),is.null(part$probability))
  cfg$M<-6L;cfg$E$nuisance_kernel_rank<-3L
  s<-kdr_geometry(dat,cfg,333L);s$B<-diag(6L)
  cfg$nuisance$raw_spec<-sim_E_nuisance_spec(dat,s,cfg)
  a<-kdr_nuisance_raw(dat$A,dat$W,cfg);altered<-cfg;altered$baseline$truth<-NULL;altered$baseline$reference<-NULL
  stopifnot(identical(a,kdr_nuisance_raw(dat$A,dat$W,altered)),
    identical(as.numeric(a[,1]),as.numeric(dat$W)),identical(as.numeric(a[,2]),as.numeric(dat$W^2)))
  fake<-list(beta=rep(.01,6),gamma=c(.1,.2));raw<-function(FUN,progress=FALSE)FUN(dat)
  score<-sim_E_score(fake,s,raw)
  eta<-drop(kdr_basis(dat$A,s)%*%fake$beta)+.1*dat$W+.2*dat$W^2
  p<-exp(eta)/drop(group_sum(exp(eta),dat$group))[dat$group]
  stopifnot(abs(score$mean_loss+sum(log(p[dat$event==1]))/dat$n)<1e-12)
  ns<-c("learning","tuning","estimation","evaluation")
  stopifnot(!anyDuplicated(vapply(ns,function(n)kdr_seed(cfg,"data_pool",n,1L,1L),integer(1))))
  dir.create(file.path(output_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  writeLines(c("Module E checks passed:","18 candidates; seven settings; independent selection and final samples",
    "Complete-fold scoring, failure exclusion, deterministic fallback",
    "Matched-set preservation, truth-free nuisance dictionary, conditional validation score"),
    file.path(output_root,"validation","module_E_checks.txt"))
  kdr_log("Module E design, fold eligibility, truth isolation and scoring checks passed.")
  invisible(TRUE)
}
