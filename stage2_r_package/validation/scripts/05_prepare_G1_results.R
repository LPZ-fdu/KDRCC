# Matched-cohort summaries and comparisons, independent of model fitting.
rp_bind<-function(x) {x<-Filter(function(z)!is.null(z)&&NROW(z)>0,x);if(!length(x))data.frame() else do.call(rbind,x)}
rp_label<-function(x,label)if(nrow(x))cbind(implementation=label,x,stringsAsFactors=FALSE) else data.frame()
rp_escape<-function(x) {x<-gsub("&","&amp;",as.character(x),fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)}
rp_table<-function(d,id,title,note,folder) {
  write.csv(d,file.path(folder,paste0(id,".csv")),row.names=FALSE,na="")
  if(!nrow(d))return(invisible(NULL))
  fmt<-function(x)if(is.numeric(x))ifelse(is.finite(x),format(signif(x,6),trim=TRUE),"") else ifelse(is.na(x),"",as.character(x))
  cells<-lapply(d,fmt)
  body<-vapply(seq_len(nrow(d)),function(i)paste0("<tr>",paste0("<td>",vapply(cells,function(x)rp_escape(x[i]),character(1)),"</td>",collapse=""),"</tr>"),character(1))
  html<-c("<!doctype html><html lang='en'><meta charset='utf-8'>",
    "<style>body{font-family:'Times New Roman',serif;font-size:12pt;margin:28px}table{border-collapse:collapse}th,td{padding:6px 9px;text-align:right}th{border-top:1pt solid black;border-bottom:.75pt solid black}tr:last-child td{border-bottom:1pt solid black}th:first-child,td:first-child{text-align:left}.note{max-width:1200px;line-height:1.5}</style>",
    paste0("<h2>",rp_escape(title),"</h2><table><thead><tr>",paste0("<th>",rp_escape(gsub("_"," ",names(d))),"</th>",collapse=""),"</tr></thead><tbody>"),
    body,paste0("</tbody></table><p class='note'><b>Note.</b> ",rp_escape(note),"</p></html>"))
  writeLines(html,file.path(folder,paste0(id,".html")),useBytes=TRUE)
}
rp_point_summary<-function(records,e) {
  out<-list()
  for(target in c("joint","individual")) {
    parts<-lapply(records,function(x)x$outputs[x$outputs$target==target,])
    z<-parts[[1]];Y<-do.call(rbind,lapply(parts,`[[`,"estimate"))
    L<-do.call(rbind,lapply(parts,`[[`,"lower"));U<-do.call(rbind,lapply(parts,`[[`,"upper"))
    S<-do.call(rbind,lapply(parts,`[[`,"se"));E<-sweep(Y,2,z$truth)
    tab<-z[,c("target","index","q","truth","eligible","interval_eligible")]
    tab$usable_repetitions<-nrow(Y);tab$mean_estimate<-colMeans(Y)
    tab$empirical_SD<-apply(Y,2,e$sim_sd);tab$bias<-colMeans(E);tab$pointwise_RMSE<-sqrt(colMeans(E^2))
    tab$lower_repetition_quantile<-apply(Y,2,quantile,.025);tab$upper_repetition_quantile<-apply(Y,2,quantile,.975)
    tab$mean_CI_lower<-apply(L,2,e$sim_mean);tab$mean_CI_upper<-apply(U,2,e$sim_mean)
    tab$RMS_model_SE<-sqrt(apply(S^2,2,e$sim_mean));tab$interval_repetitions<-colSums(is.finite(L)&is.finite(U))
    if(target=="joint") {
      i<-which(!z$eligible);tab$mean_CI_lower[i]<-tab$mean_CI_upper[i]<-tab$RMS_model_SE[i]<-0
    }
    out[[target]]<-tab
  }
  rp_bind(out)
}
rp_summarize<-function(records,plan,ev,e) {
  effects<-e$sim_effect_summary(records,nrow(plan));raw<-e$sim_interpretation_collect(records,plan)
  extended<-e$sim_extended_summary(raw)
  Y<-do.call(rbind,lapply(records,function(x)x$interpretation$response$estimate))
  response<-e$sim_B_response_summary(Y,ev$interpretation$response,plan$replicate,nrow(plan))
  L<-do.call(rbind,lapply(records,function(x)x$interpretation$response$lower));U<-do.call(rbind,lapply(records,function(x)x$interpretation$response$upper))
  response$points$mean_CI_lower<-apply(L,2,e$sim_mean)[response$points$index]
  response$points$mean_CI_upper<-apply(U,2,e$sim_mean)[response$points$index]
  response$points$interval_repetitions<-colSums(is.finite(L)&is.finite(U))[response$points$index]
  selections<-rp_bind(lapply(seq_along(records),function(i){z<-records[[i]]$selection$table
    cbind(plan[rep(i,nrow(z)),c("config_id","null_count","p","replicate","method","scenario_gf"),drop=FALSE],z)}))
  status<-plan;status$selection_ok<-vapply(records,function(x)isTRUE(x$selection$valid),logical(1))
  select<-e$sim_selection_summary(selections,status)
  list(effects=effects,points=rp_point_summary(records,e),interpretation=raw,extended=extended,
    response=response,selection=select,selection_raw=selections)
}
rp_diagnostic_table<-function(records,label) {
  getnum<-function(x)if(is.null(x)||!length(x))NA_real_ else as.numeric(x[1])
  rp_bind(lapply(records,function(a){p<-inherits(a$fit,"kdrcc");f<-a$fit
    data.frame(implementation=label,replicate=a$job$replicate,
      iterations=if(p)nrow(f$diagnostics$estimation) else getnum(f$iterations),
      objective=if(p)getnum(f$diagnostics$objective) else getnum(f$objective),
      maximum_gradient=if(p)getnum(f$diagnostics$gradient) else getnum(f$gradient),
      minimum_curvature=if(p)getnum(f$covariance$minimum_curvature) else getnum(a$inference$minimum_curvature),
      curvature_condition=if(p)getnum(f$covariance$condition) else getnum(a$inference$hessian_condition),
      elapsed_seconds=getnum(a$elapsed_seconds),stringsAsFactors=FALSE)}))
}
rp_summary_difference<-function(reference,package) {
  out<-list()
  compare<-function(label,a,b,keys,values) {
    if(!nrow(a)||!nrow(b))return(invisible(NULL))
    za<-a[,c(keys,values),drop=FALSE];zb<-b[,c(keys,values),drop=FALSE]
    m<-merge(za,zb,by=keys,suffixes=c("_archived","_package"),all=TRUE,sort=FALSE)
    for(n in values) {
      x<-m[[paste0(n,"_archived")]];y<-m[[paste0(n,"_package")]];good<-is.finite(x)&is.finite(y)
      out[[length(out)+1L]]<<-data.frame(result=label,quantity=n,comparison_rows=nrow(m),
        finite_pairs=sum(good),maximum_absolute_difference=if(any(good))max(abs(x[good]-y[good])) else NA_real_,
        mean_absolute_difference=if(any(good))mean(abs(x[good]-y[good])) else NA_real_)
    }
  }
  compare("Joint and observed combinations",reference$points,package$points,c("target","index"),c("mean_estimate","empirical_SD","mean_CI_lower","mean_CI_upper"))
  compare("Single and interaction curves",reference$response$points,package$response$points,c("target","component","index"),c("mean_estimate","pointwise_RMSE","mean_CI_lower","mean_CI_upper"))
  compare("Effect accuracy",reference$effects$metrics,package$effects$metrics,"target",c("RMS_Bias","RMSE","mean_OR_absolute_relative_error","coverage_percent"))
  compare("Contribution shares",reference$extended$contributions,package$extended$contributions,"variable",c("estimated_share","share_RMSE"))
  compare("Overall interactions",reference$extended$pairs,package$extended$pairs,"component","corrected_detection_percent")
  compare("Interaction directions",reference$extended$directions,package$extended$directions,"component",c("corrected_correct_percent","corrected_wrong_percent"))
  compare("Variable identification",reference$selection$variables,package$selection$variables,"variable","rate")
  rp_bind(out)
}
rp_slim<-function(x) {
  x$effects$repetitions<-NULL;x$interpretation<-NULL;x$selection_raw<-NULL
  x$response$repetitions<-NULL;x
}
rp_cleanup<-function(run) {
  base<-normalizePath(run,winslash="/",mustWork=TRUE)
  for(n in c("fits","kernels","logs","status","invocations","data")) {
    path<-file.path(base,n)
    if(!dir.exists(path))next
    resolved<-normalizePath(path,winslash="/",mustWork=TRUE)
    if(dirname(tolower(resolved))!=tolower(base)||!basename(resolved)%in%c("fits","kernels","logs","status","invocations","data"))
      stop("Unexpected checkpoint directory: ",resolved)
    unlink(resolved,recursive=TRUE)
    if(dir.exists(resolved))stop("Unable to remove completed checkpoint directory: ",resolved)
  }
  worker_log<-file.path(base,"workers.log")
  if(file.exists(worker_log))unlink(worker_log)
}
rp_render<-function(root,summary,analysis,settings,e) {
  sys.source(file.path(root,"scripts/06_plot_G1_results.R"),environment())
  rp_figures(summary,analysis,settings,e)
  rp_write_captions(analysis,length(summary$replicates))
}
rp_report<-function(root,run,settings) {
  rp_package_guard();e<-rp_engine(root,report=TRUE)
  final<-file.path(run,"analysis/summary.rds")
  if(file.exists(final)) {
    saved<-readRDS(final)
    if(all(saved$status$complete)&&!length(list.files(file.path(run,"fits"),pattern="^r[0-9]+[.]rds$"))) {
      rp_render(root,saved,dirname(final),settings,e)
      rp_cleanup(run)
      rp_log("Final G1 figures regenerated from summary data:",dirname(final))
      return(invisible(dirname(final)))
    }
  }
  signature<-readRDS(file.path(run,"signature.rds"));reps<-signature$replicates
  plan<-readRDS(file.path(root,"reference/stage_one_G1/plan.rds"));plan<-plan[match(reps,plan$replicate),,drop=FALSE]
  cfg<-readRDS(file.path(root,"reference/stage_one_G1/configuration.rds"));ev<-readRDS(file.path(root,"reference/stage_one_G1/evaluation.rds"))
  refs<-lapply(reps,function(r)rp_reference(root,r)$record)
  packages<-lapply(reps,function(r){p<-file.path(run,"fits",sprintf("r%03d.rds",r));if(!file.exists(p))return(NULL)
    tryCatch(readRDS(p),error=function(err)NULL)})
  analysis<-file.path(run,"analysis");tables<-file.path(analysis,"tables");figures<-file.path(analysis,"figures")
  dir.create(tables,recursive=TRUE,showWarnings=FALSE);dir.create(figures,recursive=TRUE,showWarnings=FALSE)
  status<-rp_bind(lapply(seq_along(reps),function(i){x<-packages[[i]]
    data.frame(replicate=reps[i],completed=!is.null(x),fit_ok=isTRUE(x$ok),
      covariance_ok=isTRUE(x$interval_ok),interpretation_ok=isTRUE(x$interpretation_ok),
      selection_ok=isTRUE(x$selection$valid),interaction_ok=isTRUE(x$interpretation$detection$ok),
      error=if(is.null(x))"Pending" else if(is.null(x$error))"" else x$error,stringsAsFactors=FALSE)}))
  status$complete<-with(status,fit_ok&covariance_ok&interpretation_ok&selection_ok&interaction_ok)
  availability<-data.frame(component=c("Planned repetitions","Completed fits","Coefficient covariance","Response and contribution outputs","Variable tests","Interaction tests","Complete matched repetitions"),
    count=c(length(reps),sum(status$fit_ok),sum(status$covariance_ok),sum(status$interpretation_ok),sum(status$selection_ok),sum(status$interaction_ok),sum(status$complete)))
  rp_table(availability,"Table_G1_availability","Availability of matched G1 reproductions",
    "Counts refer to the selected original repetition cohort. Completed matched repetitions determine all comparisons.",tables)
  checks<-rp_bind(lapply(seq_along(packages),function(i){x<-packages[[i]]
    if(is.null(x$agreement))return(NULL);cbind(replicate=reps[i],x$agreement)}))
  if(nrow(checks)) {
    ag<-rp_bind(lapply(split(checks,checks$check),function(x)data.frame(check=x$check[1],repetitions=nrow(x),passed=sum(x$passed),
      maximum_absolute_difference=max(x$maximum_absolute_difference),smallest_tolerance=min(x$tolerance),largest_tolerance=max(x$tolerance))))
    rp_table(ag,"agreement_summary","Agreement across matched repetitions","Counts retain every matched repetition. Maximum differences and fixed tolerance ranges summarize numerical agreement; decisions use zero tolerance.",tables)
  }
  complete<-with(status,fit_ok&covariance_ok&interpretation_ok&selection_ok&interaction_ok)
  if(settings$require_complete_cohort&&!all(complete)) {
    write.csv(status,file.path(run,"pending_status.csv"),row.names=FALSE)
    stop("The selected cohort requires completed fits and inference for every repetition. Inspect Table_G1_availability.csv, resume failed or pending tasks, and rerun PREPARE_RESULTS.R.")
  }
  use<-which(complete);if(!length(use))stop("Completed matched reproductions are required for summary figures.")
  original<-rp_summarize(refs[use],plan[use,,drop=FALSE],ev,e)
  reproduced<-rp_summarize(packages[use],plan[use,,drop=FALSE],ev,e)
  messages<-unlist(lapply(packages[use],function(x)unique(x$warnings)),use.names=FALSE)
  warnings<-if(length(messages))as.data.frame(table(messages),stringsAsFactors=FALSE) else data.frame(messages=character(),Freq=integer())
  names(warnings)<-c("warning","repetitions")
  rp_table(warnings,"Table_G1_execution_warnings","Execution warning frequencies",
    "A warning is counted at most once per repetition. Frequencies summarize the completed package refits.",tables)
  both<-function(a,b)rp_bind(list(rp_label(a,"Archived G1"),rp_label(b,"R package")))
  note<-paste("The same",length(use),"original repetition IDs enter both implementations. Effects use log OR; OR error is 100*abs(exp(estimate-truth)-1).")
  rp_table(both(original$effects$metrics,reproduced$effects$metrics),"Table_G1_effect_accuracy","Joint and observed-combination effect accuracy",note,tables)
  rp_table(both(original$points,reproduced$points),"Table_G1_pointwise_effects","Pointwise joint and observed-combination estimates",
    paste(note,"Mean confidence limits average fitted pointwise interval endpoints separately; repetition limits are empirical 2.5th and 97.5th percentiles."),tables)
  rp_table(both(original$response$points,reproduced$response$points),"Table_G1_response_curves","Single-exposure and interaction response estimates",note,tables)
  rp_table(both(original$response$metrics,reproduced$response$metrics),"Table_G1_response_accuracy","Single-exposure and interaction response accuracy",note,tables)
  rp_table(both(original$extended$contributions,reproduced$extended$contributions),"Table_G1_contribution_shares","Active-exposure contribution shares",
    "Shares normalize evaluation-distribution mean absolute baseline Shapley allocations over the eight specified active exposures and sum to 100% within a fit.",tables)
  rp_table(both(original$extended$OR,reproduced$extended$OR),"Table_G1_quantitative_OR","Quantitative joint and single-exposure odds ratios",
    "OR estimates and mean interval endpoints use the same exposure contrasts and reference values. Error and coverage summaries use the generating truth.",tables)
  rp_table(both(original$extended$pairs,reproduced$extended$pairs),"Table_G1_interaction_detection","Overall interaction detection",
    "Holm correction covers the four prespecified pairs within each fit. Rates average over matched independent repetitions; Monte Carlo limits use Wilson intervals.",tables)
  rp_table(both(original$extended$directions,reproduced$extended$directions),"Table_G1_interaction_direction","Significant interaction directions",
    "Directions require a significant pointwise contrast and its true sign. Holm correction covers all eligible interaction grid contrasts; rates first average within a repetition.",tables)
  rp_table(both(original$selection$variables,reproduced$selection$variables),"Table_G1_variable_identification","Variable participation identification rates",
    "Identification uses the derivative quadratic statistic and Holm-adjusted p values. Rates are proportions of matched repetitions; limits are Wilson Monte Carlo intervals.",tables)
  difference<-rp_summary_difference(original,reproduced)
  rp_table(difference,"Table_G1_summary_agreement","Agreement of aggregated G1 results",
    "Absolute differences compare the same aggregation, reference combinations, concentration grids and repetition cohort.",tables)
  diagnostics<-rp_bind(list(rp_diagnostic_table(refs[use],"Archived G1"),rp_diagnostic_table(packages[use],"R package")))
  diagnostic_dir<-file.path(analysis,"model_diagnostics");dir.create(diagnostic_dir,showWarnings=FALSE)
  diagnostic_summary<-rp_bind(lapply(split(diagnostics,diagnostics$implementation),function(x)
    rp_bind(lapply(setdiff(names(x),c("implementation","replicate")),function(n) {
      v<-x[[n]][is.finite(x[[n]])];if(!length(v))return(NULL)
      data.frame(implementation=x$implementation[1],quantity=n,repetitions=length(v),
        mean=mean(v),median=median(v),SD=if(length(v)>1)sd(v) else NA_real_,minimum=min(v),maximum=max(v))
    }))))
  rp_table(diagnostic_summary,"Table_G1_model_diagnostics","Numerical estimation diagnostics",
    "Curvature and gradient refer to the final regularized criterion. Quantities are summarized across matched repetitions.",diagnostic_dir)
  summary<-list(reference=rp_slim(original),package=rp_slim(reproduced),configuration=cfg,evaluation=ev,
    replicates=reps[use],status=status[c("replicate","complete")],availability=availability,
    agreement=if(exists("ag"))ag else data.frame(),differences=difference,diagnostics=diagnostics,warnings=warnings,
    execution=list(finished=Sys.time(),session=capture.output(sessionInfo())))
  rp_save(summary,file.path(analysis,"summary.rds"))
  rp_render(root,summary,analysis,settings,e)
  writeLines(c("# Matched reproduction of formal KDR-CC G1 results","",
    paste("Selected original repetition IDs:",paste(reps,collapse=", ")),
    paste("Completed matched summaries:",length(use),"/",length(reps)),
    "Fixed configuration: N=20,000; L=20,000; M=144; kappa=0.0001; lambda=1.","",
    "Each installed-package model is learned from its regenerated original matched data. The center, penalty-row and auxiliary-fold streams follow the archived study. Original saved coefficients are used for comparison after fitting.","",
    "The figures show the matched-cohort mean effects, average pointwise confidence limits and between-repetition ranges. Tables retain aggregated differences, identification rates, contributions, interaction detection and quantitative ORs.","",
    paste("Per-repetition numerical checks passing:",sum(checks$passed),"/",nrow(checks)),
    paste("Generated:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")),"",capture.output(sessionInfo())),file.path(analysis,"REPORT.md"))
  html<-c("<!doctype html><html lang='en'><meta charset='utf-8'><title>KDR-CC G1 reproduction</title>",
    "<style>body{font-family:'Times New Roman',serif;max-width:1300px;margin:30px auto;padding:0 20px}img{max-width:100%}a{line-height:1.8}</style>",
    "<h1>Matched reproduction of formal KDR-CC G1 results</h1>",
    paste0("<p>",length(use)," matched independent repetitions. N=20,000; L=20,000; M=144; kappa=0.0001; lambda=1.</p>"),
    "<p><a href='REPORT.md'>Execution report</a> | <a href='Figure_and_table_captions.md'>Figure and table captions</a></p>")
  for(f in list.files(figures,pattern="[.]png$"))html<-c(html,paste0("<h2>",rp_escape(gsub("_"," ",sub("[.]png$","",f))),"</h2><a href='figures/",f,"'><img src='figures/",f,"'></a>"))
  html<-c(html,"<h2>Tables</h2>",paste0("<a href='tables/",list.files(tables,pattern="[.]html$"),"'>",list.files(tables,pattern="[.]html$"),"</a><br>"),"</html>")
  writeLines(html,file.path(analysis,"index.html"),useBytes=TRUE)
  if(all(status$complete)) {
    rp_cleanup(run)
    if(file.exists(file.path(run,"pending_status.csv")))unlink(file.path(run,"pending_status.csv"))
  }
  rp_log("Final matched G1 figures and tables:",analysis)
  invisible(analysis)
}
