# Run from the Independent Validation directory after installing kdrcc.
options(warn=1)
library(kdrcc)
dir.create("results",showWarnings=FALSE)
checks<-list();started<-proc.time()[3]
record<-function(name,error=0,tolerance=1e-8) {
  ok<-length(error)==1L&&is.finite(error)&&error<=tolerance
  checks[[length(checks)+1L]]<<-data.frame(check=name,error=error,tolerance=tolerance,passed=ok,stringsAsFactors=FALSE)
  write.csv(do.call(rbind,checks),"results/checks.csv",row.names=FALSE)
  cat(format(Sys.time(),"[%H:%M:%S]"),if(ok)"PASS" else "FAIL",name,"|",format(error),"\n")
  if(!ok)stop("Validation failed: ",name)
}
ns<-asNamespace("kdrcc");fn<-function(n)get(n,ns)
old<-new.env(parent=globalenv())
for(f in sort(list.files("reference/frozen_engine",pattern="[.]R$",full.names=TRUE)))sys.source(f,old)
old$kdr_heartbeat<-function(...)invisible(NULL);old$kdr_bar<-function(...)invisible(NULL)
# Shared finite g dictionary for a direct prepared-design engine comparison.
dat<-kdrcc_example(480L,seed=1621L,discrete=TRUE)
sp<-kdrcc_split(dat,"set_id",learning=160L,seed=1622L)
ctrl<-kdrcc_control(centers=6L,kappa=1e-3,lambda=1,kernel_tolerance=1e-4)
nu<-kdrcc_nuisance(g_exposure=~X1+X2,f_method="cells")
fit<-kdrcc(sp$estimation,exposures=c("X1","X2","X3"),adjust=~W+I(W^2),learning_data=sp$learning,nuisance=nu,control=ctrl)
record("End-to-end discrete fit",as.numeric(!fit$diagnostics$converged))
d<-fn(".kdr_validate")(sp$estimation,"case","set_id",c("X1","X2","X3"),~W+I(W^2))
s<-fit$kernel;z<-fn(".kdr_prepare")(d,s);aux<-fit$nuisance;N<-d$n
old$kdr_nuisance_raw<-function(A,W,cfg)z$Traw
pa<-z$pa;pa$DW<-pa$DS;pa$cell<-max.col(z$RF,ties.method="first")
legacy<-list(A=d$A,W=d$data$W,event=d$event,group=d$group,n=N,fold=aux$folds,pairs=pa,X=z$X)
each<-function(FUN)FUN(legacy)
oa<-list(g=lapply(aux$models,function(a){g<-a$g;g$keep<-seq_len(ncol(z$Traw));g}),
 f=list(means=lapply(aux$models,function(a)t(a$f$Pi)),
 mass=do.call(rbind,lapply(aux$models,function(a)diag(solve(a$f$gram_inverse)))),
 ntrain=vapply(aux$models,function(a)a$f$ntrain,numeric(1)),wrong=FALSE))
cfg<-list(solver=list(maxit=ctrl$maxit,tolerance=ctrl$tolerance,step_tolerance=ctrl$step_tolerance,max_backtrack=20L),penalty=list(a=ctrl$lambda,reference_N=N))
b<-coef(fit)
st<-fn(".kdr_state")(b,z,aux);ot<-old$kdr_dr_state(b,each,oa,rep(1,9),N)
record("Doubly robust moment equivalence",max(abs(st$U-ot$U)),1e-10)
record("Moment Jacobian equivalence",max(abs(st$J-ot$J)),1e-10)
of<-old$kdr_fit(each,oa,s,cfg,N)
record("Regularized coefficient equivalence",max(abs(b-of$beta)),1e-8)
oc<-old$kdr_dr_covariance(of,each,oa,s,cfg,N)
record("Full auxiliary-adjusted covariance equivalence",max(abs(vcov(fit)-oc$V)),1e-8)
cur<-fn(".kdr_curvature")(b,z,aux,s,ctrl$lambda)
eps<-1e-5
fd<-sapply(seq_along(b),function(j){v<-rep(0,length(b));v[j]<-eps
 (fn(".kdr_objective")(b+v,z,aux,s,ctrl$lambda)$gradient-fn(".kdr_objective")(b-v,z,aux,s,ctrl$lambda)$gradient)/(2*eps)})
record("Exact Hessian central difference",max(abs(cur$H-fd)),2e-6)
# Refit a continuous conditional-mean basis and compare contracted and explicit influence formulae.
cd<-kdrcc_example(480L,seed=1631L)
cp<-kdrcc_split(cd,"set_id",learning=160L,seed=1632L)
cf<-kdrcc(cp$estimation,exposures=c("X1","X2","X3"),adjust=~W+I(W^2),learning_data=cp$learning,
 nuisance=kdrcc_nuisance(g_df=2L,kernel_rank=3L,f_df=2L),control=ctrl)
dz<-fn(".kdr_validate")(cp$estimation,"case","set_id",c("X1","X2","X3"),~W+I(W^2))
zz<-fn(".kdr_prepare")(dz,cf$kernel);aa<-cf$nuisance;nn<-dz$n;bb<-coef(cf)
cv<-fn(".kdr_curvature")(bb,zz,aa,cf$kernel,ctrl$lambda)
rows<-matrix(0,length(zz$pa$D),length(bb))
for(q in cv$state$parts)rows[q$ix,]<-(q$R*(q$w*q$t))%*%t(cv$JW)-q$Z*(q$w*q$e*drop(q$R%*%cv$r0))
SS<-fn("group_sum")(rows,zz$pa$group);ff<-aa$folds[zz$pa$group]
for(k in seq_along(aa$models)) {
 ax<-aa$models[[k]];g<-ax$g
 et<-drop(zz$Traw%*%g$raw_beta);mx<-as.numeric(tapply(et,zz$d$group,max));ee<-exp(et-mx[zz$d$group]);pr<-ee/drop(fn("group_sum")(ee,zz$d$group))[zz$d$group]
 IG<-fn("group_sum")((zz$Traw%*%g$influence_map)*(zz$d$event-pr),zz$d$group);IG[aa$folds==k,]<-0
 SS<-SS+(nn/g$ntrain)*IG%*%t(cv$sensitivity[[k]]$Ag)
 ix<-which(ff!=k&zz$pa$D==0);RF<-zz$RF[ix,,drop=FALSE];R<-zz$pa$Z[ix,,drop=FALSE]-RF%*%t(ax$f$Pi)
 for(j in seq_len(ncol(RF))) {
   Ab<--cv$sensitivity[[k]]$v[j]*cv$JW+outer(cv$sensitivity[[k]]$U[,j],cv$r0)
   IFr<-R*(zz$pa$w[ix]*drop(RF%*%ax$f$gram_inverse[,j]))
   ids<-sort(unique(zz$pa$group[ix]))
   SS[ids,]<-SS[ids,,drop=FALSE]+(nn/ax$f$ntrain)*fn("group_sum")(IFr,zz$pa$group[ix])%*%t(Ab)
 }
}
SSc<-sweep(SS,2,colMeans(SS));Hi<-solve(cv$H);VV<-Hi%*%crossprod(SSc)%*%t(Hi)/(nn*(nn-1))
record("Continuous f full influence contraction",max(abs(VV-vcov(cf))),1e-8)
# Protective and harmful responses, reference identity and serialization.
pred<-predict(cf,cp$estimation[1:10,]);saveRDS(cf,"results/example_model.rds",version=2)
record("Saved model prediction identity",max(abs(pred$log_or-predict(readRDS("results/example_model.rds"),cp$estimation[1:10,])$log_or)),1e-12)
ref<-as.data.frame(as.list(cf$reference));zero<-predict(cf,ref)
record("Reference contrast and variance",max(abs(c(zero$log_or,zero$se))),1e-12)
K<-cf$kernel;X<-as.matrix(cp$estimation[1:20,c("X1","X2","X3")]);j<-2L
h<-1e-5;xp<-xm<-X;xp[,j]<-xp[,j]+h*K$scale[j];xm[,j]<-xm[,j]-h*K$scale[j]
num<-(fn("kdr_basis")(xp,K)-fn("kdr_basis")(xm,K))/(2*h)
record("Standardized feature derivative",max(abs(num-fn("kdr_derivative")(X,K,j))),2e-7)
sh<-kdrcc_contribution(cf,data=as.data.frame(X));en<-old$sim_shapley_enumerate(X,cf$reference,function(A)drop(fn("kdr_basis")(A,K)%*%coef(cf)))
record("Exact Shapley versus subset enumeration",max(abs(sh$allocations-en$total)),1e-9)
record("Contribution allocation identity",sh$identity_error,1e-9)
vt<-kdrcc_variables(cf);record("Holm calculation",max(abs(vt$data$p_adjusted-p.adjust(vt$data$p_value,"holm"))),1e-12)
it<-kdrcc_interaction(cf,matrix(c("X1","X3"),1,2),points=3L)
record("Interaction reference axes",max(abs(it$data$log_or[!it$data$eligible])),1e-12)
# Independent selection and final estimation; repeat by serial and PSOCK execution.
td<-kdrcc_example(600L,seed=220L);tp<-kdrcc_split(td,"set_id",learning=180L,tuning=240L,seed=221L)
nu2<-kdrcc_nuisance(g_df=2L,kernel_rank=3L,f_df=2L);grid<-expand.grid(M=c(6L,8L),kappa=1e-3,lambda=c(.3,1))
tc<-ctrl;tc$seed<-27182818L
ta<-kdrcc_tune(tp$tuning,tp$learning,exposures=c("X1","X2","X3"),adjust=~W+I(W^2),nuisance=nu2,grid=grid,control=tc)
tc$workers<-2L
tb<-kdrcc_tune(tp$tuning,tp$learning,exposures=c("X1","X2","X3"),adjust=~W+I(W^2),nuisance=nu2,grid=grid,control=tc)
record("Serial and PSOCK selection scores",max(abs(ta$grid$mean_loss-tb$grid$mean_loss)),1e-10)
record("Serial and PSOCK selected configuration",max(abs(ta$selected-tb$selected)),1e-12)
tf<-kdrcc(tp$estimation,exposures=c("X1","X2","X3"),adjust=~W+I(W^2),kernel=ta$kernel,control=ta$control)
record("Selected-kernel independent final fit",as.numeric(is.null(tf$tuning)||!tf$diagnostics$converged))
write.csv(ta$grid,"results/tuning_candidates.csv",row.names=FALSE)
# Saved formal Module B fit: public package outputs against archived analysis outputs.
fixture<-readRDS("reference/formal_B_fixture.rds")
obj<-structure(list(coefficients=fixture$beta,covariance=list(matrix=fixture$V),kernel=fixture$kernel,
 reference=fixture$evaluation$reference,evaluation=as.data.frame(fixture$evaluation$individual),
 specification=list(exposures=colnames(fixture$evaluation$individual),control=kdrcc_control())),class="kdrcc")
for(target in c("joint","individual")) {
 A<-fixture$evaluation[[target]];p<-predict(obj,A);q<-fixture$outputs[fixture$outputs$target==target,]
 record(paste("Formal B",target,"prediction"),max(abs(p$log_or-q$estimate)),1e-10)
 ix<-which(is.finite(q$se));record(paste("Formal B",target,"standard error"),max(abs(p$se[ix]-q$se[ix])),1e-9)
}
fv<-kdrcc_variables(obj)
record("Formal B variable statistics",max(abs(fv$data$statistic-fixture$selection$score)),1e-9)
record("Formal B variable p values",max(abs(fv$data$p_value-fixture$selection$p_value)),2e-7)
record("Formal B Holm decisions",as.numeric(!identical(fv$data$selected,fixture$selection$selected)))
fc<-kdrcc_contribution(obj,data=as.data.frame(fixture$evaluation$individual[fixture$contribution_ids,,drop=FALSE]))
record("Formal B active contribution shares",max(abs(fc$data$share_percent-fixture$contributions$estimated_share)),1e-7)
og<-fixture$evaluation$interpretation$OR;po<-predict(obj,og$A,reference=og$B)
record("Formal B quantitative OR",max(abs(po$or-fixture$OR$OR)),1e-9)
saveRDS(list(continuous=cf,discrete=fit,variables=vt,interaction=it,contribution=sh,tuning=ta),"results/validation_objects.rds",version=2)
write.csv(pred,"results/example_predictions.csv",row.names=FALSE)
write.csv(vt$data,"results/example_variable_tests.csv",row.names=FALSE)
write.csv(sh$data,"results/example_contributions.csv",row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),paste("Validation elapsed seconds:",round(proc.time()[3]-started,3))),"results/session_info.txt")
cat("ALL",length(checks),"CHECKS PASSED\n")
