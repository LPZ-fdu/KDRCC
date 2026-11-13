# Refit archived formal replicates through the installed kdrcc package.
rp_log<-function(...) {cat(format(Sys.time(),"[%H:%M:%S]"),...,"\n");flush.console()}
rp_save<-function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile("record_",dirname(path));saveRDS(x,tmp,version=2)
  if(!file.copy(tmp,path,overwrite=TRUE))stop("Checkpoint write failed: ",path)
  unlink(tmp);invisible(path)
}
rp_engine<-function(root,report=FALSE) {
  e<-new.env(parent=globalenv())
  for(f in sort(list.files(file.path(root,"reference/formal_simulation"),pattern="[.]R$",full.names=TRUE)))sys.source(f,e)
  e$kdr_bar<-function(...)invisible(NULL);e$kdr_heartbeat<-function(...)invisible(NULL)
  e$sim_parameter_names<-c("N","L","M","kappa","lambda")
  if(report)for(f in sort(list.files(file.path(root,"reference/formal_report"),pattern="[.]R$",full.names=TRUE)))sys.source(f,e)
  e
}
rp_reference<-function(root,r)readRDS(file.path(root,"reference/stage_one_G1/fits",sprintf("r%03d.rds",r)))
rp_complete<-function(x)isTRUE(x$ok)&&isTRUE(x$interval_ok)&&isTRUE(x$interpretation_ok)&&
  isTRUE(x$selection$valid)&&isTRUE(x$interpretation$detection$ok)
rp_package_guard<-function() {
  if(!requireNamespace("kdrcc",quietly=TRUE))stop("Install kdrcc before running this validation.")
  if(as.character(utils::packageVersion("kdrcc"))!="0.1.0")stop("This adapter specifies kdrcc 0.1.0; verify the adapter for the installed package version.")
  for(n in c(".kdr_validate",".kdr_learn","kdr_kernel_penalty","kdr_derivative","kdr_vroot",".kdr_tail"))
    if(!exists(n,asNamespace("kdrcc"),inherits=FALSE))stop("Required installed-package component is missing: ",n)
}
rp_status<-function(run,r,stage,detail="")rp_save(list(replicate=r,stage=stage,detail=detail,time=Sys.time()),file.path(run,"status",sprintf("r%03d.rds",r)))
rp_maxdiff<-function(a,b) {
  if(length(a)!=length(b)||!identical(dim(a),dim(b)))return(Inf)
  if(any(is.finite(a)!=is.finite(b)))return(Inf)
  i<-is.finite(a)&is.finite(b);if(!any(i))return(if(all(is.na(a)&is.na(b)))0 else Inf)
  max(abs(a[i]-b[i]))
}
rp_callback<-function(cfg) {
  e<-new.env(parent=asNamespace("stats"))
  e$spec<-cfg$baseline;e$names<-c("PM25","PMcoarse","O3","NO2","SO2","CO","temperature","humidity")
  e$anchor<-function(x,a) {
    f<-splinefun(a$x,a$percent,method="monoH.FC");lo<-min(a$x);hi<-max(a$x);tau<-(hi-lo)/3
    y<-f(pmax(lo,pmin(hi,x)));l<-x<lo;u<-x>hi
    y[l]<-f(lo)+f(lo,deriv=1)*tau*(exp((x[l]-lo)/tau)-1)
    y[u]<-f(hi)+f(hi,deriv=1)*tau*(1-exp(-(x[u]-hi)/tau))
    if(any(y<=-100))stop("Invalid response-surface tail")
    log1p(y/100)
  }
  environment(e$anchor)<-e
  f<-function(data) {
    A<-as.matrix(data[names]);z<-sweep(sweep(A,2,spec$reference$center),2,spec$reference$scale,"/")
    main<-sapply(1:6,function(j)anchor(A[,j],spec$truth$anchors[[j]]))
    t<-(A[,7]-spec$truth$temperature$optimum)/10;h<-(A[,8]-60)/20
    weather<-cbind(pmin(t,0)^2,pmax(t,0)^2,h^2);pp<-combn(8L,2L)
    ix<-sapply(seq_len(ncol(pp)),function(i)z[,pp[1,i]]*z[,pp[2,i]])
    cbind(main,weather,z,ix)
  }
  environment(f)<-e;f
}
rp_frame<-function(d) {
  data.frame(set_id=d$group,case=d$event,as.data.frame(d$A),W=d$W,check.names=FALSE)
}

rp_generate_role<-function(e,cfg,n,role,r,run) {
  bs<-cfg$execution$block_groups;parts<-vector("list",ceiling(n/bs));done<-0L
  for(b in seq_along(parts)) {
    x<-e$kdr_generate(bs,8L,cfg,e$kdr_seed(cfg,"data_pool",role,r,b))
    parts[[b]]<-e$sim_subset(x,min(bs,n-(b-1L)*bs))
    done<-done+parts[[b]]$n;rp_status(run,r,"Data generation",paste(role,done,"/",n,"matched sets"))
  }
  d<-list()
  for(name in c("A","W","date","event","truth","probability","fold"))
    d[[name]]<-if(name=="A")do.call(rbind,lapply(parts,`[[`,name)) else do.call(c,lapply(parts,`[[`,name))
  offsets<-c(0L,head(cumsum(vapply(parts,function(x)x$n,integer(1))),-1L))
  d$group<-unlist(Map(function(x,o)x$group+o,parts,offsets),use.names=FALSE);d$n<-as.integer(n);d
}
rp_control<-function(cfg,job,seed)kdrcc::kdrcc_control(
  centers=job$M,kappa=job$kappa,lambda=job$lambda,folds=cfg$nuisance$folds,
  coverage=cfg$kernel$coverage_fraction,seed=seed,workers=1L,verbose=TRUE,
  maxit=cfg$solver$maxit,tolerance=cfg$solver$tolerance,step_tolerance=cfg$solver$step_tolerance,
  rho_bounds=cfg$kernel$rho_bounds,initial_scale=cfg$kernel$initial_scale,
  kernel_maxit=cfg$kernel$maxit,kernel_continuations=cfg$kernel$continuations,
  kernel_tolerance=cfg$kernel$projected_gradient_tol,penalty_rows=cfg$kernel$penalty_rows,
  evaluation_points=cfg$selection$evaluation_points,confidence=cfg$evaluation$confidence,
  tail_accuracy=cfg$selection$accuracy,tail_terms=cfg$selection$max_terms,
  tail_retry_terms=cfg$selection$retry_terms,tail_draws=cfg$selection$fallback_min_draws)
rp_selection<-function(fit,ev,cfg,job,e) {
  ans<-kdrcc::kdrcc_variables(fit,data=as.data.frame(ev$individual),alpha=cfg$selection$alpha,adjust=cfg$selection$adjust)
  z<-ans$data;z$default_stream_p_value<-z$p_value
  # A fallback uses the original replicate-and-variable named stream.
  for(j in which(z$algorithm=="Monte Carlo")) {
    A<-ev$individual;D<-getFromNamespace("kdr_derivative","kdrcc")(A,fit$kernel,j)
    Q<-crossprod(D)/nrow(D);L<-getFromNamespace("kdr_vroot","kdrcc")(vcov(fit))
    w<-pmax(eigen((crossprod(L,Q%*%L)+t(crossprod(L,Q%*%L)))/2,symmetric=TRUE,only.values=TRUE)$values,0)
    a<-getFromNamespace(".kdr_tail","kdrcc")(z$statistic[j],w,fit$specification$control,
      nrow(z),e$kdr_seed(cfg,"selection",job$key,j),alpha=cfg$selection$alpha)
    z$p_value[j]<-a$p;z$accuracy[j]<-a$accuracy;z$algorithm[j]<-a$method;z$draws[j]<-a$draws
  }
  z$p_adjusted<-p.adjust(z$p_value,cfg$selection$adjust);z$selected<-z$p_adjusted<cfg$selection$alpha
  lo<-p.adjust(pmax(0,z$p_value-z$accuracy),cfg$selection$adjust)<cfg$selection$alpha
  hi<-p.adjust(pmin(1,z$p_value+z$accuracy),cfg$selection$adjust)<cfg$selection$alpha
  z$numerically_stable<-lo==hi;z$selected[!z$numerically_stable]<-NA
  table<-data.frame(variable=z$variable,true_active=TRUE,score=z$statistic,p_value=z$p_value,
    algorithm=z$algorithm,absolute_accuracy=z$accuracy,draws=z$draws,p_adjusted=z$p_adjusted,
    selected=z$selected,numerically_stable=z$numerically_stable,stringsAsFactors=FALSE)
  list(table=table,valid=all(!is.na(z$selected)),error=if(all(!is.na(z$selected)))"" else "Numerically unresolved participation decisions",public_output=z)
}
rp_outputs<-function(fit,ev,cfg) {
  out<-list()
  for(target in c("joint","individual")) {
    a<-predict(fit,ev[[target]],reference=ev$reference,level=cfg$evaluation$confidence);truth<-ev[[paste0("true_",target)]]
    eligible<-if(target=="joint")abs(ev$quantile-.5)>1e-12 else rep(TRUE,length(truth))
    interval<-if(target=="joint")eligible else seq_along(truth)%in%ev$coverage_ids
    a$se[!interval]<-a$lower[!interval]<-a$upper[!interval]<-NA_real_
    signal<-abs(truth)>cfg$evaluation$zero_tolerance
    out[[target]]<-data.frame(target=target,index=seq_along(truth),q=if(target=="joint")ev$quantile else NA_real_,
      truth=truth,estimate=a$log_or,se=a$se,lower=a$lower,upper=a$upper,eligible=eligible,interval_eligible=interval,
      covered=ifelse(interval&is.finite(a$se),a$lower<=truth&a$upper>=truth,NA),
      correct_significant=ifelse(signal,(truth>0&a$lower>0)|(truth<0&a$upper<0),NA),
      wrong_significant=ifelse(signal,(truth>0&a$upper<0)|(truth<0&a$lower>0),NA),
      or_absolute_relative_error=100*abs(expm1(a$log_or-truth)))
  }
  do.call(rbind,out)
}
rp_interpretation<-function(fit,ev,cfg,e) {
  vars<-names(ev$reference);grid<-ev$interpretation$response;pt<-grid$points
  si<-which(pt$target=="single");ii<-which(pt$target=="interaction")
  values<-setNames(lapply(vars,function(v)pt$concentration_1[pt$target=="single"&pt$exposure_1==v]),vars)
  single<-kdrcc::kdrcc_single(fit,values=values,background=ev$reference,level=cfg$evaluation$confidence)$data
  stopifnot(nrow(single)==length(si),identical(single$exposure,pt$exposure_1[si]),rp_maxdiff(single$concentration,pt$concentration_1[si])<1e-10)
  pairs<-grid$settings$B_interaction_pairs;pairs<-matrix(vars[pairs],ncol=2L)
  vi<-setNames(lapply(vars,function(v)sort(unique(c(pt$concentration_1[ii][pt$exposure_1[ii]==v],pt$concentration_2[ii][pt$exposure_2[ii]==v])))),vars)
  interaction<-kdrcc::kdrcc_interaction(fit,pairs,values=vi,background=ev$reference,
    alpha=cfg$interpretation$alpha,level=cfg$evaluation$confidence)
  ip<-interaction$data
  stopifnot(nrow(ip)==length(ii),identical(ip$pair,pt$component[ii]),
    rp_maxdiff(ip$concentration_1,pt$concentration_1[ii])<1e-10,rp_maxdiff(ip$concentration_2,pt$concentration_2[ii])<1e-10)
  response<-pt
  for(n in c("estimate","SE","lower","upper"))response[[n]]<-NA_real_
  response$estimate[si]<-single$log_or;response$SE[si]<-single$se;response$lower[si]<-single$lower;response$upper[si]<-single$upper
  response$estimate[ii]<-ip$log_or;response$SE[ii]<-ip$se;response$lower[ii]<-ip$lower;response$upper[ii]<-ip$upper
  signal<-abs(pt$truth)>cfg$evaluation$zero_tolerance
  sig<-response$lower>0|response$upper<0
  response$correct_significant<-ifelse(signal,sig&sign(response$estimate)==sign(pt$truth),NA)
  response$wrong_significant<-ifelse(signal,sig&sign(response$estimate)!=sign(pt$truth),NA)
  response$covered<-response$lower<=pt$truth&response$upper>=pt$truth
  ix<-which(ip$eligible);points<-pt[ii[ix],,drop=FALSE];points$estimate<-ip$log_or[ix];points$SE<-ip$se[ix]
  points$lower<-ip$lower[ix];points$upper<-ip$upper[ix];points$p_value<-ip$p_value[ix];points$p_holm_grid<-ip$p_adjusted[ix]
  points$nominal_detected<-points$p_value<cfg$interpretation$alpha;points$holm_grid_detected<-points$p_holm_grid<cfg$interpretation$alpha
  points$true_signal<-abs(points$truth)>cfg$evaluation$zero_tolerance;points$point_ok<-is.finite(points$SE)&points$SE>0
  for(mode in c("nominal","holm_grid")) {
    detected<-points[[paste0(mode,"_detected")]];correct<-sign(points$estimate)==sign(points$truth)
    points[[paste0(mode,"_correct_direction")]]<-ifelse(points$true_signal,detected&correct,NA)
    points[[paste0(mode,"_wrong_direction")]]<-ifelse(points$true_signal,detected&!correct,NA)
  }
  overall<-interaction$overall;first<-match(overall$pair,points$component)
  pair_table<-data.frame(component=overall$pair,exposure_1=points$exposure_1[first],exposure_2=points$exposure_2[first],
    evaluation_points=as.integer(table(factor(points$component,levels=overall$pair))),true_interaction=TRUE,
    statistic=overall$statistic,p_value=overall$p_value,p_holm_pairs=overall$p_adjusted,
    raw_detected=overall$p_value<cfg$interpretation$alpha,holm_pairs_detected=overall$selected,
    numerical_accuracy=overall$accuracy,algorithm=overall$algorithm,
    test_ok=is.finite(overall$p_value)&!is.na(overall$selected),stringsAsFactors=FALSE)
  ids<-ev$interpretation$ids;con<-kdrcc::kdrcc_contribution(fit,data=as.data.frame(ev$individual[ids,,drop=FALSE]),variables=vars)
  contribution<-data.frame(variable=vars,true_active=TRUE,true_share=ev$interpretation$true_contribution_share,
    estimated_share=con$data$share_percent,absolute_contribution=con$data$mean_absolute_contribution)
  og<-ev$interpretation$OR;a<-predict(fit,og$A,reference=og$B,level=cfg$evaluation$confidence);or<-og$points
  or$estimate<-a$log_or;or$SE<-a$se;or$lower<-a$lower;or$upper<-a$upper
  signal<-abs(or$truth)>cfg$evaluation$zero_tolerance;sig<-a$lower>0|a$upper<0
  or$correct_significant<-ifelse(signal,sig&sign(a$log_or)==sign(or$truth),NA)
  or$wrong_significant<-ifelse(signal,sig&sign(a$log_or)!=sign(or$truth),NA)
  or$covered<-a$lower<=or$truth&a$upper>=or$truth
  or$OR<-a$or;or$OR_lower<-a$or_lower;or$OR_upper<-a$or_upper
  or$OR_absolute_relative_error<-100*abs(expm1(a$log_or-or$truth))
  list(response=response,contributions=contribution,OR=or,
    detection=list(ok=all(pair_table$test_ok),error="",value=list(pairs=pair_table,points=points)),
    contribution_identity_error=con$identity_error)
}
rp_compare_record<-function(reference,package,settings) {
  rows<-list();a<-settings$agreement
  put<-function(name,x,y,tolerance,relative=FALSE) {
    d<-rp_maxdiff(as.numeric(x),as.numeric(y));scale<-if(relative)1+max(abs(c(x,y)),na.rm=TRUE) else 1
    rows[[length(rows)+1L]]<<-data.frame(check=name,maximum_absolute_difference=d,tolerance=tolerance*scale,
      passed=is.finite(d)&&d<=tolerance*scale,stringsAsFactors=FALSE)
  }
  r<-reference$record;p<-package
  put("Joint and observed-combination log OR",r$outputs$estimate,p$outputs$estimate,a$log_or)
  put("Joint and observed-combination standard errors",r$outputs$se,p$outputs$se,a$standard_error)
  put("Joint and observed-combination confidence limits",c(r$outputs$lower,r$outputs$upper),c(p$outputs$lower,p$outputs$upper),a$log_or+2*a$standard_error)
  put("Single-exposure and four-point interaction log OR",r$interpretation$response$estimate,p$interpretation$response$estimate,a$log_or)
  put("Response standard errors",r$interpretation$response$SE,p$interpretation$response$SE,a$standard_error)
  put("Response confidence limits",c(r$interpretation$response$lower,r$interpretation$response$upper),c(p$interpretation$response$lower,p$interpretation$response$upper),a$log_or+2*a$standard_error)
  put("Regularized coefficients",r$fit$beta,p$fit$coefficients,a$coefficient_relative,TRUE)
  put("Coefficient covariance",r$inference$V,p$fit$covariance$matrix,a$covariance_relative,TRUE)
  put("Independently learned centers",reference$kernel$centers,p$fit$kernel$centers,a$log_or)
  put("Independently learned distance parameters",reference$kernel$rho,p$fit$kernel$rho,a$coefficient_relative,TRUE)
  put("Feature whitening",reference$kernel$B,p$fit$kernel$B,a$coefficient_relative,TRUE)
  put("Moment transformation",reference$kernel$T,p$fit$kernel$T,a$coefficient_relative,TRUE)
  put("Original-stream penalty matrix",reference$kernel$penalty,p$fit$kernel$penalty,a$penalty_relative,TRUE)
  put("Quantitative log-OR contrasts",r$interpretation$OR$estimate,p$interpretation$OR$estimate,a$log_or)
  put("Contribution shares (percentage points)",r$interpretation$contributions$estimated_share,p$interpretation$contributions$estimated_share,a$contribution_percentage_points)
  put("Variable quadratic statistics",r$selection$table$score,p$selection$table$score,a$coefficient_relative,TRUE)
  pv_tol<-a$p_value
  for(j in seq_len(nrow(r$selection$table))) {
    rr<-r$selection$table[j,];pp<-p$selection$table[j,]
    mc<-if(rr$draws>0)sqrt(rr$p_value*(1-rr$p_value)/rr$draws) else 0
    put(paste("Variable p value",rr$variable),rr$p_value,pp$p_value,max(pv_tol,4*mc))
  }
  put("Variable Holm decisions",as.numeric(r$selection$table$selected),as.numeric(p$selection$table$selected),0)
  put("Overall interaction p values",r$interpretation$detection$value$pairs$p_value,p$interpretation$detection$value$pairs$p_value,a$p_value)
  put("Overall interaction Holm decisions",as.numeric(r$interpretation$detection$value$pairs$holm_pairs_detected),as.numeric(p$interpretation$detection$value$pairs$holm_pairs_detected),0)
  put("Pointwise interaction Holm p values",r$interpretation$detection$value$points$p_holm_grid,p$interpretation$detection$value$points$p_holm_grid,a$p_value)
  put("Pointwise interaction directions",sign(r$interpretation$detection$value$points$estimate),sign(p$interpretation$detection$value$points$estimate),0)
  do.call(rbind,rows)
}
rp_one_task<-function(task) {
  rp_package_guard();root<-task$root;run<-task$run;r<-task$replicate;settings<-task$settings
  log<-file(file.path(run,"logs",sprintf("r%03d.log",r)),"at",encoding="UTF-8")
  sink(log);on.exit({sink();close(log)},add=TRUE)
  started<-proc.time()[3];warnings<-character();e<-rp_engine(root)
  reference<-rp_reference(root,r);job<-reference$record$job
  cfg<-readRDS(file.path(root,"reference/stage_one_G1/configuration.rds"))
  ev<-readRDS(file.path(root,"reference/stage_one_G1/evaluation.rds"))
  path<-file.path(run,"fits",sprintf("r%03d.rds",r));kernel_path<-file.path(run,"kernels",sprintf("r%03d.rds",r))
  center_seed<-e$kdr_seed(cfg,"centers",r);penalty_seed<-e$kdr_seed(cfg,"penalty",r)
  out<-withCallingHandlers(tryCatch({
    rp_status(run,r,"Data generation")
    ld<-rp_generate_role(e,cfg,job$L,"learning",r,run);d<-rp_generate_role(e,cfg,job$N,"estimation",r,run)
    stopifnot(!identical(ld$A,d$A))
    identities<-data.frame(role=c("learning","estimation"),matched_sets=c(ld$n,d$n),
      identity=c(e$object_id(ld[c("A","W","group","event","fold")]),e$object_id(d[c("A","W","group","event","fold")])))
    cfg$M<-job$M;cfg$kernel$kappa<-job$kappa
    control<-rp_control(cfg,job,center_seed);adj<-~W+I(W^2);environment(adj)<-baseenv()
    nu<-kdrcc::kdrcc_nuisance(g_exposure=rp_callback(cfg),f_method="cells")
    rp_status(run,r,"Independent package kernel learning")
    ld_valid<-getFromNamespace(".kdr_validate","kdrcc")(rp_frame(ld),"case","set_id",e$exposure_names,adj)
    if(settings$resume&&file.exists(kernel_path))kernel<-readRDS(kernel_path) else {
      kernel<-getFromNamespace(".kdr_learn","kdrcc")(ld_valid,adj,nu,control)
      kernel$penalty<-getFromNamespace("kdr_kernel_penalty","kdrcc")(ld_valid,kernel,cfg,penalty_seed)
      kernel$validation_streams<-list(centers=center_seed,penalty=penalty_seed)
      rp_save(kernel,kernel_path)
    }
    # All fitted coefficients and auxiliaries are estimated from the regenerated data.
    nu<-kernel$nuisance
    rp_status(run,r,"Package cross-fitting, regularized estimation and covariance")
    fit<-kdrcc::kdrcc(rp_frame(d),exposures=e$exposure_names,adjust=adj,kernel=kernel,nuisance=nu,
      control=control,reference=ev$reference,fold_id=d$fold)
    if(fit$covariance$status!="available")stop("Package covariance is unavailable: ",fit$covariance$error)
    rp_status(run,r,"Package predictions, contributions and tests")
    outputs<-rp_outputs(fit,ev,cfg);interpretation<-rp_interpretation(fit,ev,cfg,e)
    selection<-rp_selection(fit,ev,cfg,job,e)
    result<-list(job=job,ok=TRUE,fit=fit,outputs=outputs,interpretation=interpretation,
      interpretation_ok=TRUE,interval_ok=TRUE,selection=selection,selection_required=TRUE,
      input_identities=identities,random_streams=list(centers=center_seed,penalty=penalty_seed,
        selection=vapply(seq_len(8L),function(j)e$kdr_seed(cfg,"selection",job$key,j),integer(1))),
      elapsed_seconds=unname(proc.time()[3]-started),warnings=warnings,
      session=capture.output(sessionInfo()),R=as.character(getRversion()),package_version=as.character(utils::packageVersion("kdrcc")))
    result$agreement<-rp_compare_record(reference,result,settings)
    result
  },error=function(err)list(job=job,ok=FALSE,error=conditionMessage(err),warnings=warnings,
    elapsed_seconds=unname(proc.time()[3]-started),calls=capture.output(sys.calls()))),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  out$warnings<-unique(warnings);out$completed_at<-Sys.time();rp_save(out,path)
  rp_status(run,r,if(isTRUE(out$ok))"Complete" else "Failed",if(isTRUE(out$ok))paste("agreement checks",sum(out$agreement$passed),"/",nrow(out$agreement)) else out$error)
  list(replicate=r,ok=isTRUE(out$ok),error=if(is.null(out$error))"" else out$error,seconds=out$elapsed_seconds)
}
rp_find_run<-function(root,settings) {
  if(!is.null(settings$output_directory))return(normalizePath(settings$output_directory,winslash="/",mustWork=TRUE))
  path<-file.path(root,"results/LATEST_G1_REPRODUCTION.txt")
  if(!file.exists(path))stop("Run RUN_END_TO_END.R before preparing the reproduced results.")
  saved<-readLines(path,warn=FALSE)[1];run<-file.path(root,"results",basename(saved))
  normalizePath(run,winslash="/",mustWork=TRUE)
}
rp_run<-function(root,settings) {
  rp_package_guard();e<-rp_engine(root)
  reps<-sort(unique(as.integer(settings$replicates)))
  if(!length(reps)||anyNA(reps)||any(!reps%in%1:100))stop("replicates must select existing formal repetition IDs from 1 to 100.")
  if(length(settings$workers)!=1L||!is.finite(settings$workers)||settings$workers<1||settings$workers!=floor(settings$workers))stop("workers must be a positive integer.")
  paths<-file.path(root,"reference/stage_one_G1/fits",sprintf("r%03d.rds",reps))
  if(!all(file.exists(paths)))stop("The complete archived G1 reference cohort is required.")
  signature<-list(replicates=reps,reference_MD5=unname(tools::md5sum(paths)),
    configuration=unname(tools::md5sum(file.path(root,"reference/stage_one_G1/configuration.rds"))),
    adapter=unname(tools::md5sum(file.path(root,"scripts/04_reproduce_G1.R"))),agreement=settings$agreement,
    package_version=as.character(utils::packageVersion("kdrcc")),R=as.character(getRversion()))
  id<-substr(e$object_id(signature),1,12)
  run<-if(is.null(settings$output_directory))file.path(root,"results",paste0("G1_reproduction_",id)) else settings$output_directory
  dir.create(run,recursive=TRUE,showWarnings=FALSE);run<-normalizePath(run,winslash="/",mustWork=TRUE)
  for(n in c("fits","kernels","logs","status","analysis"))dir.create(file.path(run,n),showWarnings=FALSE)
  sigpath<-file.path(run,"signature.rds")
  if(settings$resume&&file.exists(sigpath)&&!identical(signature,readRDS(sigpath)))stop("Output identity differs; choose a new output directory or resume=FALSE.")
  if(settings$resume&&file.exists(file.path(run,"analysis/summary.rds"))) {
    saved<-readRDS(file.path(run,"analysis/summary.rds"))
    if(identical(as.integer(saved$replicates),reps)&&all(saved$status$complete)) {
      writeLines(basename(run),file.path(root,"results/LATEST_G1_REPRODUCTION.txt"))
      rp_log("Completed matched G1 summary reused:",run)
      return(invisible(run))
    }
  }
  rp_save(signature,sigpath);rp_save(settings,file.path(run,"execution_settings.rds"))
  rp_save(list(planned_replicates=reps,started=Sys.time(),settings=settings),file.path(run,"invocations",paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),".rds")))
  writeLines(basename(run),file.path(root,"results/LATEST_G1_REPRODUCTION.txt"))
  pending<-reps
  if(settings$resume)pending<-Filter(function(r){p<-file.path(run,"fits",sprintf("r%03d.rds",r));if(!file.exists(p))return(TRUE)
    x<-tryCatch(readRDS(p),error=function(e)NULL);is.null(x)||(!rp_complete(x)&&isTRUE(settings$retry_failed))},reps)
  pending<-as.integer(pending);n<-length(reps);cached<-n-length(pending);started<-Sys.time();completed<-cached
  rp_log("Formal G1 reproduction |",n,"matched original replicates | N=20000 L=20000 M=144 kappa=0.0001 lambda=1")
  rp_log("Reused checkpoints",cached,"| new tasks",length(pending),"| requested processes",settings$workers)
  print_progress<-function(detail="") {
    elapsed<-as.numeric(difftime(Sys.time(),started,units="secs"));fresh<-completed-cached
    eta<-if(fresh>0)elapsed/fresh*(n-completed) else NA_real_
    rp_log("Progress",paste0(completed,"/",n),"| elapsed",round(elapsed/60,1),"min | remaining",if(is.finite(eta))paste(round(eta/60,1),"min") else "estimating",detail)
  }
  if(length(pending)) {
    threads<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    original<-Sys.getenv(threads,unset=NA_character_)
    do.call(Sys.setenv,as.list(setNames(rep("1",length(threads)),threads)))
    on.exit({for(i in seq_along(threads))if(is.na(original[i]))Sys.unsetenv(threads[i]) else do.call(Sys.setenv,setNames(list(original[i]),threads[i]))},add=TRUE)
    workers<-min(as.integer(settings$workers),length(pending));jobs<-lapply(pending,function(r)list(root=root,run=run,replicate=r,settings=settings))
    if(workers==1L)for(job in jobs){z<-rp_one_task(job);completed<-completed+1L;print_progress(paste("replicate",z$replicate,if(z$ok)"saved" else paste("failed:",z$error)))} else {
      cl<-parallel::makePSOCKcluster(workers,outfile=file.path(run,"workers.log"));on.exit(parallel::stopCluster(cl),add=TRUE)
      parallel::clusterCall(cl,function(paths,locale,script){if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE",locale)
        if(getRversion()<"4.0.0")options(stringsAsFactors=FALSE)
        .libPaths(paths);sys.source(script,.GlobalEnv);NULL},.libPaths(),Sys.getlocale("LC_CTYPE"),file.path(root,"scripts/04_reproduce_G1.R"))
      active<-seq_len(workers);next_job<-workers+1L;slot_replicate<-pending[seq_len(workers)]
      for(i in active)parallel:::sendCall(cl[[i]],rp_one_task,list(jobs[[i]]))
      last<-Sys.time()
      while(length(active)) {
        ready<-socketSelect(lapply(cl[active],function(x)x$con),timeout=1)
        for(i in active[ready]) {
          z<-parallel:::recvResult(cl[[i]]);completed<-completed+1L
          if(inherits(z,"try-error")||!is.list(z))stop("A worker connection returned an invalid result; saved checkpoints can be resumed.")
          print_progress(paste("replicate",z$replicate,if(z$ok)"saved" else paste("failed:",z$error)))
          if(next_job<=length(jobs)){slot_replicate[i]<-jobs[[next_job]]$replicate;parallel:::sendCall(cl[[i]],rp_one_task,list(jobs[[next_job]]));next_job<-next_job+1L} else active<-setdiff(active,i)
        }
        if(as.numeric(difftime(Sys.time(),last,units="secs"))>=settings$progress_seconds) {
          print_progress()
          for(i in active){p<-file.path(run,"status",sprintf("r%03d.rds",slot_replicate[i]))
            st<-if(file.exists(p))suppressWarnings(tryCatch(readRDS(p),error=function(err)NULL)) else NULL
            if(!is.null(st))rp_log("  Replicate",st$replicate,"|",st$stage,"|",st$detail)
          }
          last<-Sys.time()
        }
      }
    }
  }
  print_progress("All scheduled tasks processed")
  rp_log("Reproduction records:",run);invisible(run)
}
