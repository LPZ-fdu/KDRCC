kdrcc_control <- function(centers=144L, kappa=1e-4, lambda=1, folds=3L,
    coverage=.25, seed=27182818L, workers=1L, verbose=FALSE,
    maxit=70L, tolerance=1e-7, step_tolerance=1e-6,
    rho_bounds=c(.01,100), initial_scale=.0625, kernel_maxit=80L,
    kernel_continuations=1L, kernel_tolerance=1e-5, penalty_rows=256L,
    evaluation_points=1000L, confidence=.95, tail_accuracy=1e-8,
    tail_terms=100000L, tail_retry_terms=1000000L, tail_draws=199999L,
    keep_data=FALSE) {
  .kdr_scalar(centers,"centers",2,TRUE);.kdr_scalar(kappa,"kappa",0,strict=TRUE)
  .kdr_scalar(lambda,"lambda");.kdr_scalar(folds,"folds",2,TRUE)
  .kdr_scalar(coverage,"coverage");if(coverage>1)stop("coverage must be at most one.")
  .kdr_scalar(seed,"seed",1,TRUE);if(seed>=2147483000)stop("seed exceeds its supported range.")
  .kdr_scalar(workers,"workers",1,TRUE);.kdr_scalar(maxit,"maxit",1,TRUE)
  for(n in c("tolerance","step_tolerance","initial_scale","kernel_tolerance","tail_accuracy"))
    .kdr_scalar(get(n),n,0,strict=TRUE)
  if(length(rho_bounds)!=2L||any(!is.finite(rho_bounds))||rho_bounds[1]<=0||rho_bounds[2]<=rho_bounds[1]||initial_scale<rho_bounds[1]||initial_scale>rho_bounds[2])stop("rho_bounds must contain initial_scale and two increasing positive values.")
  for(n in c("kernel_maxit","penalty_rows","evaluation_points","tail_terms","tail_retry_terms","tail_draws"))
    .kdr_scalar(get(n),n,1,TRUE)
  .kdr_scalar(kernel_continuations,"kernel_continuations",0,TRUE)
  if(length(confidence)!=1L||!is.finite(confidence)||confidence<=0||confidence>=1)stop("confidence must be between zero and one.")
  if(!is.logical(verbose)||length(verbose)!=1L||is.na(verbose)||!is.logical(keep_data)||length(keep_data)!=1L||is.na(keep_data))stop("verbose and keep_data must be logical scalars.")
  structure(as.list(environment())[names(formals(sys.function()))],class="kdrcc_control")
}

kdrcc_nuisance <- function(g_exposure="spline_kernel", f_method=c("basis","cells"),
    g_df=4L, kernel_rank=24L, f_df=3L, f_formula=NULL) {
  f_method<-match.arg(f_method)
  if(!(is.character(g_exposure)&&length(g_exposure)==1L&&g_exposure%in%c("spline_kernel","spline","kernel"))&&
      !inherits(g_exposure,"formula")&&!is.function(g_exposure))stop("g_exposure must specify spline_kernel, spline, kernel, a formula, or a matrix-producing function.")
  for(n in c("g_df","f_df")) .kdr_scalar(get(n),n,1,TRUE)
  .kdr_scalar(kernel_rank,"kernel_rank",0,TRUE)
  if(!is.null(f_formula)&&!inherits(f_formula,"formula"))stop("f_formula must be a formula.")
  structure(as.list(environment())[names(formals(sys.function()))],class="kdrcc_nuisance")
}

kdrcc_split <- function(data,strata,learning=.25,tuning=0,seed=27182818L) {
  .kdr_columns(data,strata)
  if(length(strata)!=1L||anyNA(data[[strata]]))stop("strata must identify one complete grouping column.")
  ids<-unique(as.character(data[[strata]]));n<-length(ids)
  count<-function(x,name){.kdr_scalar(x,name);if(x>0&&x<1)floor(n*x) else {if(x!=floor(x))stop(name," must be a fraction or an integer count.");as.integer(x)}}
  l<-count(learning,"learning");t<-count(tuning,"tuning")
  if(l<1L||l+t>=n)stop("The split requires learning and final-estimation groups.")
  ids<-with_seed(seed,sample(ids));roles<-rep("estimation",n)
  roles[seq_len(l)]<-"learning";if(t>0L)roles[l+seq_len(t)]<-"tuning"
  take<-function(role)data[as.character(data[[strata]])%in%ids[roles==role],,drop=FALSE]
  list(learning=take("learning"),tuning=take("tuning"),estimation=take("estimation"),
    assignment=data.frame(stratum=ids,role=roles,stringsAsFactors=FALSE),seed=seed)
}

.kdr_validate <- function(data,case,strata,exposures,adjust) {
  if(!inherits(adjust,"formula")||length(adjust)!=2L)stop("adjust must be a one-sided formula.")
  vars<-all.vars(adjust)
  if(length(case)!=1L||length(strata)!=1L||!length(exposures)||anyDuplicated(exposures))stop("Supply one case column, one strata column and unique exposure names.")
  if(length(intersect(c(case,strata),exposures))||length(intersect(vars,exposures)))stop("Exposure, adjustment and identification columns must have distinct roles.")
  .kdr_columns(data,unique(c(case,strata,exposures,vars)))
  used<-data[,unique(c(case,strata,exposures,vars)),drop=FALSE]
  if(anyNA(used))stop("Analysis columns contain missing values; supply complete matched sets.")
  if(!all(vapply(data[exposures],is.numeric,logical(1))))stop("Exposure columns must be numeric.")
  A<-as.matrix(data[exposures]);storage.mode(A)<-"double"
  if(any(!is.finite(A)))stop("Exposure values must be finite.")
  y<-data[[case]];if(!is.numeric(y)&&!is.logical(y))stop("case must contain numeric zero and one.")
  if(any(!y%in%c(0,1)))stop("case must contain zero and one.")
  ids<-as.character(data[[strata]]);g<-match(ids,unique(ids));o<-order(g,seq_along(g));g<-g[o]
  n<-max(g);counts<-tabulate(g)
  if(any(counts<2L)||any(drop(group_sum(y[o],g))!=1))stop("Every matched set must contain one event and at least one reference date.")
  list(data=data[o,,drop=FALSE],A=A[o,,drop=FALSE],event=as.integer(y[o]),group=g,n=n,
    ids=unique(ids),order=o,counts=counts,adjust_vars=vars)
}

.kdr_formula_fit <- function(formula,data) {
  mf<-model.frame(formula,data=data,na.action=na.fail,drop.unused.levels=FALSE)
  tt<-terms(mf);X<-model.matrix(tt,mf)
  list(terms=tt,xlevels=stats::.getXlevels(tt,mf),contrasts=attr(X,"contrasts"),columns=colnames(X),X=X)
}
.kdr_formula_new <- function(spec,data) {
  mf<-model.frame(spec$terms,data=data,na.action=na.fail,xlev=spec$xlevels)
  X<-model.matrix(spec$terms,mf,contrasts.arg=spec$contrasts)
  if(!identical(colnames(X),spec$columns))stop("The covariate design columns differ from the learned specification.")
  X
}
.kdr_basis_spec <- function(X,df) {
  if(!ncol(X))return(list())
  lapply(seq_len(ncol(X)),function(j){
    x<-X[,j];u<-sort(unique(x))
    if(length(u)<=df+1L)return(list(kind="linear",center=mean(x),scale=sd(x)))
    list(kind="spline",knots=unique(as.numeric(quantile(x,seq_len(max(1L,df-1L))/df))),
      boundary=range(x))
  })
}
.kdr_basis_new <- function(X,spec) {
  if(!length(spec))return(matrix(0,nrow(X),0L))
  do.call(cbind,lapply(seq_along(spec),function(j){s<-spec[[j]]
    if(s$kind=="linear")matrix((X[,j]-s$center)/if(is.finite(s$scale)&&s$scale>0)s$scale else 1,ncol=1L)
    else splines::ns(X[,j],knots=s$knots,Boundary.knots=s$boundary,intercept=FALSE)}))
}
.kdr_pair_indices <- function(d) {
  cases<-which(d$event==1L);controls<-which(d$event==0L);g<-d$group[controls]
  ev<-cases[g];plus<-c(ev,controls);minus<-c(controls,ev);gg<-c(g,g)
  list(plus=plus,minus=minus,group=gg,D=c(rep(1,length(g)),rep(0,length(g))),
    w=1/(2*(d$counts[gg]-1)),DS=d$S[plus,,drop=FALSE]-d$S[minus,,drop=FALSE])
}
.kdr_row_key <- function(data) {
  if(!ncol(data))return(rep("all",nrow(data)))
  vapply(seq_len(nrow(data)),function(i)paste(vapply(data,function(x){v<-as.character(x[i]);paste0(nchar(v),":",v)},character(1)),collapse="|"),character(1))
}
.kdr_pair_data <- function(d,pa,vars) {
  if(!length(vars))return(data.frame(row.names=seq_along(pa$plus)))
  x<-d$data[pa$plus,vars,drop=FALSE];y<-d$data[pa$minus,vars,drop=FALSE]
  names(x)<-paste0("plus.",vars);names(y)<-paste0("minus.",vars)
  cbind(x,y)
}
.kdr_design_fit <- function(d,s,adjust,nu) {
  a<-.kdr_formula_fit(adjust,d$data);X<-a$X
  keep<-which(colnames(X)!="(Intercept)" & colSums(kdr_center(X,d$group)^2)>1e-20)
  if(length(keep)) {
    qrS<-qr(kdr_center(X[,keep,drop=FALSE],d$group),tol=1e-12)
    keep<-sort(keep[qrS$pivot[seq_len(qrS$rank)]])
  }
  S<-X[,keep,drop=FALSE];a$keep<-keep;a$X<-NULL
  spec<-list(adjust=a,nu=nu,exposures=colnames(d$A),vars=d$adjust_vars)
  if(inherits(nu$g_exposure,"formula")) {
    spec$g_kind<-"formula";spec$g_formula<-.kdr_formula_fit(nu$g_exposure,d$data)
    spec$g_formula$X<-NULL
  } else if(is.function(nu$g_exposure)) {spec$g_kind<-"function";spec$g_fun<-nu$g_exposure}
  else {spec$g_kind<-nu$g_exposure;spec$g_splines<-.kdr_basis_spec(d$A,nu$g_df)}
  spec$thin<-s;spec$thin$B<-s$B[,seq_len(min(ncol(s$B),nu$kernel_rank)),drop=FALSE]
  G<-.kdr_g_raw(d,spec);T<-cbind(S,G);nc<-ncol(S)
  score<-kdr_center(T,d$group);scale<-sqrt(colMeans(score^2));eligible<-which(scale>1e-14)
  if(!length(eligible))stop("The auxiliary design has no within-set information.")
  Q<-qr(sweep(score[,eligible,drop=FALSE],2,scale[eligible],"/"),tol=1e-12)
  sel<-sort(eligible[Q$pivot[seq_len(Q$rank)]])
  if(nc && !all(seq_len(nc)%in%sel))stop("The exposure and adjustment auxiliary columns are aliased.")
  spec$g_keep<-sel;spec$g_B<-kdr_whiten(T[,sel,drop=FALSE],d$group)
  spec$g_C<-matrix(0,nc,ncol(spec$g_B))
  if(nc)spec$g_C<-spec$g_B[match(seq_len(nc),sel),,drop=FALSE]
  d$S<-S;d$counts<-tabulate(d$group);pa<-.kdr_pair_indices(d)
  pd<-.kdr_pair_data(d,pa,spec$vars)
  spec$f_kind<-nu$f_method
  if(nu$f_method=="cells") {
    spec$f_levels<-sort(unique(.kdr_row_key(pd)))
  } else if(!is.null(nu$f_formula)) {
    spec$f_formula<-.kdr_formula_fit(nu$f_formula,pd);spec$f_formula$X<-NULL
  } else {
    spec$f_raw_formula<-if(ncol(pd))reformulate(names(pd)) else ~1
    fraw<-.kdr_formula_fit(spec$f_raw_formula,pd);XX<-fraw$X[,colnames(fraw$X)!="(Intercept)",drop=FALSE]
    fraw$keep<-which(colnames(fraw$X)!="(Intercept)");fraw$X<-NULL
    spec$f_raw<-fraw;spec$f_splines<-.kdr_basis_spec(XX,nu$f_df)
  }
  RF<-.kdr_f_raw(d,pa,spec);qrF<-qr(RF[pa$D==0,,drop=FALSE],tol=1e-12)
  spec$f_keep<-sort(qrF$pivot[seq_len(qrF$rank)])
  spec$adjust_columns<-colnames(S)
  spec
}
.kdr_g_raw <- function(d,spec) {
  if(spec$g_kind=="formula")X<-.kdr_formula_new(spec$g_formula,d$data)
  else if(spec$g_kind=="function")X<-as.matrix(spec$g_fun(d$data))
  else {
    X<-matrix(0,nrow(d$A),0L)
    if(spec$g_kind%in%c("spline","spline_kernel"))X<-cbind(X,.kdr_basis_new(d$A,spec$g_splines))
    if(spec$g_kind%in%c("kernel","spline_kernel")&&ncol(spec$thin$B))X<-cbind(X,kdr_basis(d$A,spec$thin))
  }
  if(nrow(X)!=nrow(d$A)||!is.numeric(X)||any(!is.finite(X)))stop("g_exposure must return a finite numeric matrix with one row per date.")
  X
}
.kdr_f_raw <- function(d,pa,spec) {
  pd<-.kdr_pair_data(d,pa,spec$vars)
  if(spec$f_kind=="cells") {
    id<-match(.kdr_row_key(pd),spec$f_levels)
    if(anyNA(id))stop("A covariate-pair cell was absent from the learning specification.")
    X<-matrix(0,length(id),length(spec$f_levels));X[cbind(seq_along(id),id)]<-1;X
  } else if(!is.null(spec$f_formula)) .kdr_formula_new(spec$f_formula,pd)
  else {
    X<-.kdr_formula_new(spec$f_raw,pd)[,spec$f_raw$keep,drop=FALSE]
    cbind(1,.kdr_basis_new(X,spec$f_splines))
  }
}
