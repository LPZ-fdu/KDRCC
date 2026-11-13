# Numerical fitting diagnostics, reported separately from scientific performance.
sim_report_model_diagnostics <- function(run,result,captions,module) {
  base<-file.path(run,"analysis","model_diagnostics");fig<-file.path(base,"figures");tab<-file.path(base,"tables")
  dir.create(fig,recursive=TRUE,showWarnings=FALSE);dir.create(tab,showWarnings=FALSE)
  plan<-result$plan$tasks;records<-traces<-list()
  for(i in seq_len(nrow(plan))) {
    a<-sim_read(file.path(run,"fits",paste0(plan$key[i],".rds")));j<-plan[i,c("key","config_id","method","scenario_gf","null_count","replicate"),drop=FALSE]
    get<-function(name)if(is.null(a$fit[[name]]))NA_real_ else a$fit[[name]]
    V<-a$inference$V;if(is.null(V))V<-a$fit$V
    eig<-if(is.matrix(V)&&all(is.finite(V)))eigen(kdr_sym(V),symmetric=TRUE,only.values=TRUE)$values else NA_real_
    records[[i]]<-cbind(j,data.frame(completed=!is.null(a),point_ok=isTRUE(a$ok),interval_ok=isTRUE(a$interval_ok),
      interpretation_ok=isTRUE(a$interpretation_ok),iterations=get("iterations"),objective=get("objective"),
      gradient=get("gradient"),equation_residual=get("residual"),jacobian_condition=get("jacobian_condition"),
      covariance_min_eigenvalue=if(any(is.finite(eig)))min(eig) else NA_real_,
      error=if(is.null(a))"Pending" else a$error,
      interpretation_error=if(is.null(a$interpretation_error))"" else a$interpretation_error,
      interaction_error=if(is.null(a$interpretation$detection$error))"" else a$interpretation$detection$error))
    h<-a$fit$history
    if(!is.null(h)&&nrow(h))traces[[i]]<-cbind(j[rep(1,nrow(h)),,drop=FALSE],h)
  }
  records<-kdr_bind(records);traces<-kdr_bind(traces);rho<-learning<-list()
  for(k in unique(plan$kernel_key[plan$method=="KDR"])) {
    s<-sim_read(file.path(run,"kernels",paste0(k,".rds")));if(is.null(s))next
    j<-plan[match(k,plan$kernel_key),c("replicate","p","L","M","kappa"),drop=FALSE]
    p<-length(s$rho);b<-result$config$kernel$rho_bounds
    lower_bound<-if(module%in%c("D","E"))s$base_rho*b[1] else b[1]
    upper_bound<-if(module%in%c("D","E"))s$base_rho*b[2] else b[2]
    rho[[length(rho)+1L]]<-cbind(j[rep(1,p),,drop=FALSE],data.frame(kernel_key=k,variable=kdr_names(p),rho=s$rho,
      lower_boundary=abs(log(s$rho/lower_bound))<=1e-6,upper_boundary=abs(log(s$rho/upper_bound))<=1e-6))
    d<-s$learning$diagnostics
    if(nrow(d))learning[[length(learning)+1L]]<-cbind(j[rep(1,nrow(d)),,drop=FALSE],kernel_key=k,d)
  }
  rho<-kdr_bind(rho);learning<-kdr_bind(learning)
  local<-new.env(parent=emptyenv())
  emit<-function(d,name,title,note)sim_export_table(d,paste0("Table_",module,"_model_",name),title,note,tab,local)
  emit(records,"status","Numerical fitting and covariance diagnostics",
    "One row per planned fit. Point convergence, interval availability, interpretation availability and recorded errors are separate fields. Penalized convergence is assessed by the objective gradient and scaled step; the unpenalized equation residual is reported alongside it.")
  emit(learning,"kernel_learning","Outcome-supervised kernel optimization",
    "Each independent kernel initialization reports its objective, projected gradient, optimizer status and bound hits. Shared kernels appear once rather than once per nuisance scenario.")
  emit(rho,"distance_parameters","Estimated exposure distance parameters",
    "Positive distance parameters control Gaussian kernel smoothness. Bound indicators use a 1e-6 tolerance on the log scale. Null variables have their own distance parameters.")
  if(nrow(traces))kdr_csv(traces,file.path(tab,"iteration_history.csv"))
  cfg<-result$config
  if(nrow(traces)) {
    groups<-unique(paste(traces$method,traces$scenario_gf,sep=" | "));id<-paste0("Figure_",module,"_model_convergence");title<-"Optimization convergence"
    sim_figure(file.path(fig,paste0(id,".png")),title,ceiling(length(groups)/2),min(2,length(groups)),cfg,function(){
      for(g in groups) {
        d<-traces[paste(traces$method,traces$scenario_gf,sep=" | ")==g,];its<-sort(unique(d$iteration))
        q<-vapply(its,function(i){z<-d$gradient[d$iteration==i];z<-z[is.finite(z)];if(length(z))quantile(pmax(z,1e-16),c(.1,.5,.9)) else rep(NA_real_,3)},numeric(3))
        good<-colSums(is.finite(q))==3;xx<-its[good];q<-q[,good,drop=FALSE]
        if(!length(xx)){plot.new();next}
        plot(range(c(xx,xx+1e-6)),range(q),type="n",log="y",xlab="Iteration",ylab="Maximum gradient magnitude",main=g)
        polygon(c(xx,rev(xx)),c(q[1,],rev(q[3,])),col=adjustcolor("#0072B2",.15),border=NA)
        lines(xx,q[2,],col="#0072B2",lwd=2);abline(h=cfg$solver$tolerance,lty=3,col="grey40")
      }
    })
    local[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note="The line is the median gradient magnitude at each recorded iteration; shading is the 10th to 90th percentile range among fits reaching that iteration. The horizontal line is the stopping gradient tolerance. Different objective functions are shown in separate panels.")
  }
  if(nrow(rho)) {
    rho$label<-ifelse(rho$variable%in%exposure_names,rho$variable,"null")
    z<-aggregate(rho[c("lower_boundary","upper_boundary")],list(variable=rho$label),mean)
    order<-match(c(exposure_names,"null"),z$variable);z<-z[order[!is.na(order)],]
    id<-paste0("Figure_",module,"_model_kernel_boundaries");title<-"Kernel distance-parameter boundary rates"
    sim_figure(file.path(fig,paste0(id,".png")),title,1,1,cfg,function(){
      par(mar=c(9,6.6,3.2,1.4));barplot(100*t(as.matrix(z[c("lower_boundary","upper_boundary")])),beside=TRUE,
        names.arg=ifelse(z$variable=="null","Null exposures",sim_label(z$variable)),las=2,
        col=c("#0072B2","#D55E00"),ylab="Boundary rate (%)",ylim=c(0,110),border=NA)
      legend("topright",c("Lower bound","Upper bound"),fill=c("#0072B2","#D55E00"),bty="n")
    });local[[id]]<-list(type="Figure",file=paste0("figures/",id,".png"),title=title,
      note="Rates summarize unique fitted kernels. The null-exposure bar pools all null-variable parameters. These numerical diagnostics describe fitted smoothness and are separate from the variable significance tests.")
  }
  cap<-kdr_bind(lapply(sort(ls(local)),function(id){a<-local[[id]]
    captions[[id]]<-list(type=a$type,file=paste0("model_diagnostics/",a$file),title=a$title,note=a$note)
    data.frame(id=id,type=a$type,file=a$file,title=a$title,note=a$note)}))
  kdr_csv(cap,file.path(base,"Figure_and_table_captions.csv"))
  if(nrow(cap))writeLines(c("<!doctype html><meta charset='utf-8'><h1>Model diagnostics</h1>",
    unlist(lapply(seq_len(nrow(cap)),function(i)c(paste0("<h2>",sim_escape(cap$title[i]),"</h2>"),
      if(cap$type[i]=="Figure")paste0("<img style='max-width:100%' src='",cap$file[i],"'>") else paste0("<a href='",cap$file[i],"'>Table</a>"),
      paste0("<p>",sim_escape(cap$note[i]),"</p>"))))),file.path(base,"index.html"),useBytes=TRUE)
  invisible(records)
}
