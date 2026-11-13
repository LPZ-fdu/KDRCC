.kdr_subset_data <- function(d,ids) {
  rows<-which(d$group%in%ids)
  out<-d;out$data<-d$data[rows,,drop=FALSE];out$A<-d$A[rows,,drop=FALSE]
  out$event<-d$event[rows];out$group<-match(d$group[rows],ids)
  ord<-order(out$group);out$data<-out$data[ord,,drop=FALSE];out$A<-out$A[ord,,drop=FALSE]
  out$event<-out$event[ord];out$group<-out$group[ord];out$n<-length(ids)
  out$counts<-tabulate(out$group);out$ids<-d$ids[ids];out
}
.kdr_cv_task <- function(job) {
  ctrl<-job$control;ctrl$centers<-job$grid$M[1];ctrl$kappa<-job$grid$kappa[1]
  s<-tryCatch(.kdr_learn(job$learning,job$adjust,job$nuisance,ctrl),error=function(e)e)
  rows<-list()
  for(k in seq_len(job$K)) {
    train<-.kdr_subset_data(job$selection,which(job$fold!=k))
    valid<-.kdr_subset_data(job$selection,which(job$fold==k))
    prep<-tryCatch({
      if(inherits(s,"error"))stop(s)
      z<-.kdr_prepare(train,s)
      f<-with_seed(ctrl$seed+2L,sample(rep(seq_len(ctrl$folds),length.out=train$n)))
      list(z=z,aux=.kdr_auxiliary(z,s,f,ctrl),valid=.kdr_prepare(valid,s))
    },error=function(e)e)
    for(i in seq_len(nrow(job$grid))) {
      cc<-ctrl;cc$lambda<-job$grid$lambda[i];start<-proc.time()[3]
      result<-tryCatch({
        if(inherits(prep,"error"))stop(prep)
        fit<-.kdr_solve(prep$z,prep$aux,s,cc);vd<-prep$valid
        X<-cbind(vd$X,vd$d$S)
        sc<-kdr_cl_block(c(fit$beta,fit$gamma),X,valid$event,valid$group,FALSE)
        if(!is.finite(sc$loss))stop("The held-out conditional negative log likelihood is nonfinite.")
        cv<-tryCatch(.kdr_covariance(fit$beta,prep$z,prep$aux,s,cc$lambda),error=identity)
        if(inherits(cv,"error")) {
          list(loss=sc$loss,error="",ok=TRUE,covariance_ok=FALSE,
            covariance_error=conditionMessage(cv),curvature_condition=NA_real_)
        } else {
          list(loss=sc$loss,error="",ok=TRUE,covariance_ok=TRUE,
            covariance_error="",curvature_condition=cv$condition)
        }
      },error=function(e)list(loss=NA_real_,error=conditionMessage(e),ok=FALSE,
        covariance_ok=NA,covariance_error="",curvature_condition=NA_real_))
      rows[[length(rows)+1L]]<-data.frame(candidate=job$grid$candidate[i],fold=k,
        n=valid$n,training_n=train$n,loss_sum=result$loss,ok=result$ok,
        covariance_ok=result$covariance_ok,covariance_error=result$covariance_error,
        curvature_condition=result$curvature_condition,
        error=result$error,elapsed_seconds=unname(proc.time()[3]-start),stringsAsFactors=FALSE)
    }
  }
  list(scores=kdr_bind(rows),kernel=if(inherits(s,"error"))NULL else s)
}
kdrcc_tune <- function(data,learning_data,case="case",strata="set_id",exposures,
    adjust=~1,nuisance=kdrcc_nuisance(),
    grid=expand.grid(M=c(96L,144L),kappa=c(3e-5,1e-4),lambda=c(.3,1,3)),
    control=kdrcc_control(),validation_folds=3L,fallback=TRUE) {
  .kdr_scalar(validation_folds,"validation_folds",2,TRUE)
  if(!is.data.frame(grid)||!all(c("M","kappa","lambda")%in%names(grid))||!nrow(grid))stop("grid requires M, kappa and lambda columns.")
  grid<-unique(grid[,c("M","kappa","lambda"),drop=FALSE]);grid<-grid[order(grid$M,grid$kappa,grid$lambda),,drop=FALSE]
  for(i in seq_len(nrow(grid))) {
    .kdr_scalar(grid$M[i],"M",2,TRUE);.kdr_scalar(grid$kappa[i],"kappa",0,strict=TRUE);.kdr_scalar(grid$lambda[i],"lambda")
  }
  grid$candidate<-seq_len(nrow(grid));grid$profile<-match(paste(grid$M,grid$kappa),unique(paste(grid$M,grid$kappa)))
  ld<-.kdr_validate(learning_data,case,strata,exposures,adjust);td<-.kdr_validate(data,case,strata,exposures,adjust)
  if(td$n<=validation_folds*control$folds)stop("The selection sample is too small for nested matched-set folds.")
  if(identical(ld$A,td$A)&&identical(ld$event,td$event))stop("Learning and selection data must use distinct matched sets.")
  fold<-with_seed(control$seed+3L,sample(rep(seq_len(validation_folds),length.out=td$n)))
  jobs<-lapply(unique(grid$profile),function(j)list(grid=grid[grid$profile==j,,drop=FALSE],learning=ld,
    selection=td,adjust=adjust,nuisance=nuisance,control=control,K=validation_folds,fold=fold))
  old<-options(kdrcc.trace=control$verbose);on.exit(options(old),add=TRUE)
  if(control$workers>1L && length(jobs)>1L) {
    cl<-parallel::makePSOCKcluster(min(control$workers,length(jobs)));on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(paths,locale){
      if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE",locale)
      .libPaths(paths);loadNamespace("kdrcc");NULL
    },.libPaths(),Sys.getlocale("LC_CTYPE"))
    result<-parallel::parLapplyLB(cl,jobs,.kdr_cv_task)
  } else result<-lapply(jobs,.kdr_cv_task)
  scores<-kdr_bind(lapply(result,`[[`,"scores"));grid$mean_loss<-Inf
  grid$score_complete<-grid$eligible<-FALSE
  grid$covariance_folds<-0L
  dimension<-vapply(result,function(x)if(is.null(x$kernel))NA_integer_ else x$kernel$dimension,integer(1))
  grid$kernel_dimension<-dimension[match(grid$profile,unique(grid$profile))]
  for(i in seq_len(nrow(grid))) {
    z<-scores[scores$candidate==grid$candidate[i],,drop=FALSE]
    grid$score_complete[i]<-nrow(z)==validation_folds&&all(z$ok)&&all(is.finite(z$loss_sum))
    grid$covariance_folds[i]<-sum(z$covariance_ok %in% TRUE)
    grid$eligible[i]<-grid$score_complete[i]&&grid$covariance_folds[i]==validation_folds
    if(grid$score_complete[i])grid$mean_loss[i]<-sum(z$loss_sum)/sum(z$n)
  }
  selected<-which.min(ifelse(grid$eligible,grid$mean_loss,Inf));used_fallback<-!any(grid$eligible)
  if(used_fallback) {
    if(!isTRUE(fallback))stop("No configuration passed held-out scoring and covariance checks in every fold.")
    warning("No configuration passed every fold. Using the supplied control configuration without validation.",call.=FALSE)
    best<-c(M=control$centers,kappa=control$kappa,lambda=control$lambda)
    s<-.kdr_learn(ld,adjust,nuisance,control)
  } else {
    best<-unlist(grid[selected,c("M","kappa","lambda")],use.names=TRUE)
    s<-result[[match(grid$profile[selected],unique(grid$profile))]]$kernel
  }
  control$centers<-unname(best["M"]);control$kappa<-unname(best["kappa"]);control$lambda<-unname(best["lambda"])
  record<-list(candidate_grid=grid,fold_scores=scores,selected=best,validation_folds=fold,
    learning_n=ld$n,selection_n=td$n,fallback=used_fallback,
    criterion="Matched-set-weighted held-out conditional negative log likelihood among candidates with covariance available in every training fold")
  s$tuning_record<-record
  structure(list(kernel=s,control=control,grid=grid,scores=scores,selected=best,
    fallback=used_fallback,record=record),class="kdrcc_tuning")
}
print.kdrcc_tuning <- function(x,...) {cat("KDR-CC configuration selection\n");print(x$grid);cat("Selected:\n");print(x$selected);invisible(x)}
