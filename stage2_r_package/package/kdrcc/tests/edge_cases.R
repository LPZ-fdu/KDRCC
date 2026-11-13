library(kdrcc)
# One exposure and a factor-valued adjustment design.
d<-kdrcc_example(480L,seed=1801L,discrete=TRUE)
d$category<-factor(d$W,levels=c(-1,0,1))
sp<-kdrcc_split(d,"set_id",learning=160L,seed=1802L)
cc<-kdrcc_control(centers=6L,kappa=1e-3,lambda=1,kernel_tolerance=1e-4)
fit<-kdrcc(sp$estimation,exposures="X1",adjust=~category,learning_data=sp$learning,
 nuisance=kdrcc_nuisance(g_exposure=~X1,f_method="cells"),control=cc)
stopifnot(fit$diagnostics$converged,fit$covariance$status=="available")
stopifnot(kdrcc_contribution(fit)$identity_error<1e-8)
row<-sp$estimation[1:10,]
stopifnot(max(abs(predict(fit,row,reference=row)$log_or))<1e-12)
permuted<-sp$estimation[rev(seq_len(nrow(sp$estimation))),]
# A supplied learned object preserves its auxiliary specification during reuse.
f2<-kdrcc(permuted,exposures="X1",adjust=~category,kernel=fit$kernel,control=cc)
stopifnot(f2$diagnostics$converged,identical(f2$specification$nuisance,fit$specification$nuisance))
