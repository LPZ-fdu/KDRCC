# Complete public-package example; outputs are saved beside this script.
run_kdrcc_example <- function() {
  script<-NULL
  for(fr in sys.frames())if(!is.null(fr$ofile))script<-fr$ofile
  args<-grep("^--file=",commandArgs(FALSE),value=TRUE)
  if(is.null(script)&&length(args))script<-sub("^--file=","",args[1])
  root<-if(is.null(script))getwd() else dirname(normalizePath(script,winslash="/",mustWork=TRUE))
  library(kdrcc)
  out<-file.path(root,"output");dir.create(out,showWarnings=FALSE)
  d<-read.csv(file.path(root,"example_matched_data.csv"),stringsAsFactors=FALSE)
  vars<-c("X1","X2","X3","X4")
  sizes<-table(d$set_id);events<-tapply(d$case,d$set_id,sum)
  stopifnot(all(sizes>=2L),all(events==1L),all(complete.cases(d)))
  parts<-kdrcc_split(d,"set_id",learning=250L,seed=703L)
  control<-kdrcc_control(centers=10L,kappa=1e-3,lambda=1,
    kernel_tolerance=1e-4,verbose=TRUE)
  nuisance<-kdrcc_nuisance(g_df=2L,kernel_rank=4L,f_df=2L)
  ref<-vapply(parts$learning[vars],median,numeric(1))
  fit<-kdrcc(parts$estimation,exposures=vars,adjust=~W+I(W^2),
    learning_data=parts$learning,nuisance=nuisance,control=control,reference=ref)
  # Supply the evaluation distribution used for quantiles and default summaries.
  fit$evaluation<-parts$estimation[seq_len(min(1000L,nrow(parts$estimation))),vars,drop=FALSE]
  joint<-kdrcc_joint(fit)
  single<-kdrcc_single(fit)
  predictions<-predict(fit,parts$estimation[1:20,],type="or")
  variables<-kdrcc_variables(fit)
  interaction<-kdrcc_interaction(fit,matrix(c("X1","X3"),1L,2L),points=7L)
  contribution<-kdrcc_contribution(fit,variables=c("X1","X2","X3"))
  a<-as.data.frame(as.list(ref));b<-a
  a$X1<-unname(quantile(parts$learning$X1,.75))
  b$X1<-unname(quantile(parts$learning$X1,.25))
  contrast<-kdrcc_contrast(fit,a,reference=b)
  diagnostics<-kdrcc_diagnostics(fit)
  objects<-list(joint=joint,single=single,predictions=predictions,variables=variables,
    interaction=interaction,contribution=contribution,contrast=contrast,diagnostics=diagnostics)
  write.csv(parts$assignment,file.path(out,"sample_assignment.csv"),row.names=FALSE)
  write.csv(joint$data,file.path(out,"joint_effect.csv"),row.names=FALSE)
  write.csv(single$data,file.path(out,"single_effects.csv"),row.names=FALSE)
  write.csv(predictions,file.path(out,"observed_combination_predictions.csv"),row.names=FALSE)
  write.csv(variables$data,file.path(out,"variable_participation.csv"),row.names=FALSE)
  write.csv(interaction$overall,file.path(out,"interaction_overall.csv"),row.names=FALSE)
  write.csv(interaction$data,file.path(out,"interaction_points.csv"),row.names=FALSE)
  write.csv(contribution$data,file.path(out,"contribution_shares.csv"),row.names=FALSE)
  write.csv(contrast$data,file.path(out,"specified_contrast.csv"),row.names=FALSE)
  saveRDS(fit,file.path(out,"fitted_model.rds"),version=2)
  saveRDS(objects,file.path(out,"analysis_objects.rds"),version=2)
  writeLines(c(capture.output(summary(fit)),capture.output(sessionInfo())),file.path(out,"session_info.txt"))
  figure<-function(name,width,height,draw) {
    png(file.path(out,name),width=width,height=height,res=200,pointsize=12)
    on.exit(dev.off());draw()
  }
  figure("example_joint.png",1600,1100,function()plot(joint))
  figure("example_single.png",2200,1600,function()plot(single))
  figure("example_contribution.png",1600,1100,function()plot(contribution))
  figure("example_interaction.png",1600,1100,function()plot(interaction))
  figure("example_diagnostics.png",2000,1000,function()plot(fit))
  cat("Example outputs:",out,"\n")
  invisible(list(fit=fit,results=objects,output=out))
}
example_analysis<-run_kdrcc_example()
