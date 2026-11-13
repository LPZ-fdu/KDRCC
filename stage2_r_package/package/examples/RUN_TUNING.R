# Matched-set validation with separate learning, selection and final data.
library(kdrcc)
d<-kdrcc_example(1200L,seed=8501L)
parts<-kdrcc_split(d,"set_id",learning=300L,tuning=450L,seed=8502L)
vars<-c("X1","X2","X3","X4")
nu<-kdrcc_nuisance(g_df=2L,kernel_rank=4L,f_df=2L)
grid<-expand.grid(M=c(8L,12L),kappa=c(1e-3,3e-3),lambda=c(.3,1))
cc<-kdrcc_control(centers=8L,kappa=1e-3,lambda=1,workers=2L)
selection<-kdrcc_tune(parts$tuning,learning_data=parts$learning,exposures=vars,
  adjust=~W+I(W^2),nuisance=nu,grid=grid,control=cc,validation_folds=3L,fallback=FALSE)
fit_tuned<-kdrcc(parts$estimation,exposures=vars,adjust=~W+I(W^2),
  kernel=selection$kernel,control=selection$control)
print(selection)
summary(fit_tuned)
