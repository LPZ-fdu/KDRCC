run_unit_checks <- function(root,output_root=root) {
  source(file.path(root,"00_parameters.R"),encoding="UTF-8")
  cfg<-simulation_config(root);cfg$batch<-"unit-check";cfg$total_dimension<-108L
  checks<-list();check<-function(name,value){stopifnot(isTRUE(value));checks[[name]]<<-TRUE}
  check("Default design has 50 unique configurations",nrow(sim_design_A(cfg)$configurations)==50L)
  check("Seed is independent of process scheduling",identical(kdr_seed(cfg,"task",5),kdr_seed(cfg,"task",5)))
  alt<-cfg;alt$batch<-"validation"
  check("Independent module namespaces",kdr_seed(cfg,"task",5)!=kdr_seed(alt,"task",5))
  d<-kdr_generate(100L,12L,cfg,843L);e<-kdr_generate(100L,16L,cfg,843L)
  check("Null dimensions preserve existing columns",identical(d$A,e$A[,1:12]))
  check("One case per individual",all(drop(group_sum(d$event,d$group))==1L))
  strata<-paste(format(d$date,"%Y-%m"),as.POSIXlt(d$date)$wday)
  check("Calendar matching",all(vapply(split(strata,d$group),function(x)length(unique(x))==1L,logical(1))))
  pa<-kdr_pairs(d)
  check("Pair weights total one per set",max(abs(drop(group_sum(pa$w,pa$group))-1))<1e-12)
  half<-seq_len(length(pa$D)/2)
  check("Pair orientations are symmetric",identical(pa$plus[half],pa$minus[-half])&&identical(pa$minus[half],pa$plus[-half]))
  A<-matrix(rep(cfg$baseline$reference$center,3),nrow=3,byrow=TRUE)
  h<-kdr_truth(A,cfg);B<-A;B[,8]<-B[,8]+c(-10,0,10)
  check("Humidity has no standalone main effect at centered temperature",max(abs(kdr_truth(B,cfg)-h))<1e-12)
  C<-cbind(A,noise1=c(-100,0,100));check("Null exposures do not enter truth",identical(kdr_truth(C,cfg),h))
  fixture<-read.csv(file.path(root,"tests","fixtures","response_values.csv"),check.names=FALSE)
  check("Frozen response fixture",max(abs(kdr_truth(as.matrix(fixture[exposure_names]),cfg)-fixture$true_log_response))<1e-11)
  p1<-kdr_quadratic_tail(5,c(2,2),cfg,99L)
  check("Davies uses upper tail",abs(p1$p-pchisq(2.5,df=2,lower.tail=FALSE))<1e-7)
  cfg$selection$accuracy<-1e-10;p2<-kdr_quadratic_tail(5,c(2,2),cfg,99L)
  check("Tail accuracy stability",abs(p1$p-p2$p)<1e-7)
  # Analytic kernel derivatives are checked against centered finite differences.
  cfg$M<-8L;cfg$kernel$kappa<-.001;s<-kdr_geometry(d,cfg,543L);s$B<-diag(s$M)
  aa<-d$A[1:10,];j<-7;eps<-1e-5;ap<-am<-aa
  ap[,j]<-ap[,j]+eps*s$scale[j];am[,j]<-am[,j]-eps*s$scale[j]
  check("Kernel derivative in standardized coordinates",max(abs((kdr_basis(ap,s)-kdr_basis(am,s))/(2*eps)-kdr_derivative(aa,s,j)))<1e-7)
  Z<-sweep(sweep(aa,2,s$center),2,s$scale,"/")
  kk<-exp(-kdr_distance(sweep(Z,2,sqrt(s$rho),"*"),sweep(s$centers,2,sqrt(s$rho),"*")))
  kr<-exp(-drop(s$centers^2%*%s$rho));res<-seq_len(nrow(Z))*.03;bb<-seq_len(s$M)*.1
  old_grad<-vapply(seq_len(ncol(Z)),function(j){
    DD<-sweep(-s$rho[j]*outer(Z[,j],s$centers[,j],"-")^2*kk,2,-s$rho[j]*s$centers[,j]^2*kr)
    sum(res*drop(DD%*%bb))/3},numeric(1))
  new_grad<-kdr_rho_gradient(Z,s$centers,kk,s$rho,kr,res,bb,3)
  check("Vectorized distance gradient equals direct derivatives",max(abs(old_grad-new_grad))<1e-12)
  mix<-CompQuadForm::imhof(8,c(1,2,3),epsabs=1e-8,epsrel=1e-8)
  calc<-kdr_quadratic_tail(8,c(1,2,3),cfg,98L)
  check("Unequal-weight quadratic tail agrees with independent integration",abs(calc$p-mix$Qq)<1e-6)
  # Differentiation of the actual paired estimating function and penalized objective.
  z<-d;z$X<-kdr_basis(d$A,s);z$pairs<-pa
  z$pairs$Z<-z$X[pa$plus,,drop=FALSE]-z$X[pa$minus,,drop=FALSE]
  each<-function(FUN,progress=FALSE)FUN(z)
  aux<-list(g=lapply(1:3,function(i)list(gamma=c(.35,.8))),f=kdr_f_fit(each,cfg,d$n))
  beta<-seq_len(s$M)*.003;s$T<-diag(s$M);s$penalty<-diag(s$M)*.2
  obj<-kdr_objective(beta,each,aux,s,2,d$n)
  num<-vapply(seq_along(beta),function(i){bp<-bm<-beta;bp[i]<-bp[i]+eps;bm[i]<-bm[i]-eps
    (kdr_objective(bp,each,aux,s,2,d$n,FALSE)$value-kdr_objective(bm,each,aux,s,2,d$n,FALSE)$value)/(2*eps)},numeric(1))
  check("Penalized objective gradient",max(abs(obj$gradient-num))<1e-7)
  # Bias cancellation must not be mistaken for small total RMSE.
  mk<-function(r,estimate)list(ok=TRUE,job=list(replicate=as.integer(r)),outputs=data.frame(
    target="joint",index=1:2,q=c(.25,.75),truth=c(1,1),estimate=estimate,se=c(1,1),
    lower=estimate-1.96,upper=estimate+1.96,eligible=TRUE,covered=TRUE,
    correct_significant=FALSE,wrong_significant=FALSE))
  a<-sim_effect_summary(list(mk(1,c(0,2)),mk(2,c(2,0))),2)$metrics
  check("RMS bias separates cancellation from RMSE",a$RMS_Bias==0&&a$RMSE==1)
  check("Empirical variability uses independent repetitions",abs(a$RMS_empirical_SD-sqrt(2))<1e-12)
  base<-cfg;tmp<-tempfile();dir.create(tmp)
  base$configuration$fixed_file<-"absent.R";base$configuration$overrides$B<-list(N=1000L,L=1000L,M=8L,kappa=.001,lambda=3)
  rr<-sim_resolve_fixed(base,tmp,"B")
  check("Explicit manual parameter override",rr$parameters$N==1000L&&is.null(rr$base))
  x<-data.frame(config_id="x",null_count=0L,p=8L,replicate=rep(1:2,each=8),
    variable=rep(exposure_names,2),true_active=TRUE,selected=TRUE)
  st<-data.frame(config_id="x",null_count=0L,method="KDR",scenario_gf="G1",replicate=1:2)
  check("Zero-null FPR is undefined",is.na(sim_selection_summary(x,st)$overall$mean_FPR))
  pool<-file.path(output_root,"validation","unit_pool")
  cfg$execution$block_groups<-100L;kdr_save(sim_subset(d,100L),sim_pool_file(pool,"estimation",1,1))
  # Stored data pools have only the active columns; nulls are expanded by the reader.
  dd<-sim_subset(d,100L);dd$A<-dd$A[,1:8];kdr_save(dd,sim_pool_file(pool,"estimation",1,1))
  cached<-sim_raw(cfg,pool,100L,12L,1);nc<-cfg;nc$execution$raw_memory_MB<-0
  direct<-sim_raw(nc,pool,100L,12L,1);got<-list()
  cached(function(z)got[[1]]<<-z);cached(function(z)got[[2]]<<-z);direct(function(z)got[[3]]<<-z)
  check("Bounded RAM cache preserves repeated data reads",identical(got[[1]],got[[2]])&&identical(got[[1]],got[[3]]))
  # Independent subset enumeration validates exact kernel allocation and truth attribution.
  AA<-d$A[1:3,1:8,drop=FALSE];ref<-cfg$baseline$reference$center
  exact<-sim_shapley_enumerate(AA,ref,function(X)kdr_truth(X,cfg))$total
  check("Analytic true Shapley allocation",max(abs(exact-sim_true_contributions(AA,ref,cfg)))<1e-10)
  ss<-s;ss$rho<-s$rho[1:8];ss$center<-s$center[1:8];ss$scale<-s$scale[1:8];ss$centers<-s$centers[,1:8,drop=FALSE]
  b<-seq_len(s$M)/10
  kexact<-sim_shapley_enumerate(AA,ref,function(X)drop(kdr_basis(X,ss)%*%b))$total
  fast<-sim_shapley_kernel(AA,ref,ref,ss,b)$total
  check("Gaussian Shapley integral equals subset enumeration",max(abs(kexact-fast))<1e-10)
  hi<-ss;hi$rho<-rep(.03,88);hi$center<-rep(0,88);hi$scale<-rep(1,88)
  hi$centers<-matrix(sin(seq_len(s$M*88)),s$M,88);XX<-matrix(cos(seq_len(3*88)),3,88);ref88<-rep(0,88)
  hi$B<-diag(s$M);vv<-sim_shapley_kernel(XX,ref88,ref88,hi,b)$total
  dd<-drop((kdr_basis(XX,hi)-kdr_basis(matrix(ref88,3,88,byrow=TRUE),hi))%*%b)
  check("88-dimensional Shapley efficiency",max(abs(rowSums(vv)-dd))<1e-10)
  check("Fixed interaction coefficients",max(abs(cfg$baseline$truth$interactions$coefficient-
    c(.135743242415130250,.084016234171294737,.085641763170191634,.069229207135094992)))<1e-15)
  ev<-list(reference=ref,joint=rbind(ref-1,ref,ref+1),quantile=c(.1,.5,.9),individual=AA)
  cfg$interpretation$contribution_points<-3L;ctx<-sim_interpretation_context(ev,cfg)
  check("True contribution shares sum to 100",abs(sum(ctx$true_contribution_share)-100)<1e-10)
  D<-matrix(c(1,2,3,4),2);V<-diag(c(.1,.2));beta<-c(.3,.4)
  cp<-sim_contrast_prediction(D,beta,V,c(0,0),cfg)
  check("Contrast uncertainty includes covariance",max(abs(cp$SE-sqrt(diag(D%*%V%*%t(D)))))<1e-12)
  rec<-list(ok=TRUE,interval_ok=TRUE,interpretation_ok=TRUE,selection_required=FALSE,
    interpretation=list(detection=list(ok=TRUE,value=list(pairs=data.frame(test_ok=rep(TRUE,4))))))
  check("Completed inference is resumable",sim_complete(rec,TRUE))
  rec$interpretation$detection$value$pairs$test_ok[2]<-FALSE
  check("Unresolved interaction inference is retryable",!sim_complete(rec,TRUE)&&sim_complete(rec,FALSE))
  report<-data.frame(check=names(checks),passed=unlist(checks),row.names=NULL)
  dir.create(file.path(output_root,"validation"),recursive=TRUE,showWarnings=FALSE)
  kdr_csv(report,file.path(output_root,"validation","unit_check_results.csv"))
  kdr_log("Numerical and unit checks passed:",nrow(report));invisible(report)
}
