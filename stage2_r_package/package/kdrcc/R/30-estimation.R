.kdr_prepare <- function(d,s) {
  spec<-s$design
  S<-.kdr_formula_new(spec$adjust,d$data)[,spec$adjust$keep,drop=FALSE]
  d$S<-S;d$counts<-tabulate(d$group)
  pa<-.kdr_pair_indices(d);X<-kdr_basis(d$A,s)
  pa$Z<-X[pa$plus,,drop=FALSE]-X[pa$minus,,drop=FALSE]
  RF<-.kdr_f_raw(d,pa,spec)[,spec$f_keep,drop=FALSE]
  Traw<-cbind(S,.kdr_g_raw(d,spec))[,spec$g_keep,drop=FALSE]
  list(d=d,pa=pa,X=X,RF=RF,Traw=Traw,V=Traw%*%spec$g_B)
}
.kdr_auxiliary <- function(z,s,fold,control) {
  K<-max(fold);ff<-fold[z$pa$group];aux<-vector("list",K)
  cfg<-list(solver=list(maxit=control$maxit,tolerance=control$tolerance,
    step_tolerance=control$step_tolerance,max_backtrack=20L))
  for(k in seq_len(K)) {
    kdr_log("Auxiliary fold",paste0(k,"/",K))
    cand<-which(fold[z$d$group]!=k);ids<-unique(z$d$group[cand]);ng<-length(ids)
    train<-function(FUN)FUN(list(V=z$V[cand,,drop=FALSE],event=z$d$event[cand],
      group=match(z$d$group[cand],ids)))
    gf<-kdr_cl_fit(train,ncol(z$V),cfg)
    C<-s$design$g_C;B<-s$design$g_B
    g<-list(gamma=drop(C%*%gf$beta),raw_beta=drop(B%*%gf$beta),
      influence_map=B%*%kdr_inverse(gf$H)%*%t(C),ntrain=ng,
      history=gf$history,iterations=gf$iterations)
    ix<-which(ff!=k & z$pa$D==0);w<-z$pa$w[ix];RF<-z$RF[ix,,drop=FALSE]
    gram<-crossprod(RF,RF*w)/ng
    Qi<-kdr_inverse(gram);Pi<-t(Qi%*%crossprod(RF,z$pa$Z[ix,,drop=FALSE]*w)/ng)
    aux[[k]]<-list(g=g,f=list(Pi=Pi,gram_inverse=Qi,ntrain=ng,
      minimum_mass=if(s$design$f_kind=="cells")min(diag(gram)) else NA_real_))
  }
  list(folds=fold,models=aux)
}
.kdr_state <- function(beta,z,aux,derivatives=TRUE) {
  M<-length(beta);N<-z$d$n;pa<-z$pa;ff<-aux$folds[pa$group]
  U<-numeric(M);J<-matrix(0,M,M);parts<-vector("list",length(aux$models))
  for(k in seq_along(aux$models)) {
    ix<-which(ff==k);a<-aux$models[[k]]
    Z<-pa$Z[ix,,drop=FALSE];R<-Z-z$RF[ix,,drop=FALSE]%*%t(a$f$Pi)
    eta<-drop(Z%*%beta+pa$DS[ix,,drop=FALSE]%*%a$g$gamma)
    e<-numeric(length(ix));yes<-pa$D[ix]==1
    if(any(-eta[yes]>650))stop("The exponential moment exceeds its numerical range.")
    e[yes]<-exp(-eta[yes]);tt<-e-(1-pa$D[ix]);w<-pa$w[ix]
    U<-U+colSums(R*(w*tt))
    if(derivatives)J<-J-crossprod(R,Z*(w*e))
    parts[[k]]<-list(ix=ix,Z=Z,R=R,e=e,t=tt,w=w)
  }
  list(U=U/N,J=J/N,parts=parts)
}
.kdr_objective <- function(beta,z,aux,s,lambda,derivatives=TRUE) {
  a<-.kdr_state(beta,z,aux,derivatives);W<-crossprod(s$T);u<-drop(s$T%*%a$U)
  list(value=sum(u^2)/2+lambda*drop(crossprod(beta,s$penalty%*%beta))/2,
    state=a,gradient=if(derivatives)drop(crossprod(a$J,W%*%a$U)+lambda*s$penalty%*%beta) else NULL,
    gauss_newton=if(derivatives)crossprod(s$T%*%a$J)+lambda*s$penalty else NULL)
}
.kdr_solve <- function(z,aux,s,control) {
  b<-numeric(s$dimension);history<-list();converged<-FALSE
  for(it in seq_len(control$maxit)) {
    q<-.kdr_objective(b,z,aux,s,control$lambda);ok<-FALSE
    for(damp in c(0,1e-10,1e-8,1e-6,1e-4,.01)) {
      H<-q$gauss_newton+diag(damp*max(1,mean(diag(q$gauss_newton))),length(b))
      step<-tryCatch(drop(solve(H,-q$gradient)),error=function(e)NULL)
      if(is.null(step))next
      if(max(abs(q$gradient))<control$tolerance&&max(abs(step))/(1+max(abs(b)))<control$step_tolerance) {converged<-TRUE;break}
      slope<-sum(q$gradient*step);if(!is.finite(slope)||slope>=0)next
      rate<-1
      for(bt in 0:20) {
        v<-tryCatch(.kdr_objective(b+rate*step,z,aux,s,control$lambda,FALSE)$value,error=function(e)Inf)
        if(is.finite(v)&&v<=q$value+1e-4*rate*slope+1e-15){ok<-TRUE;break}
        rate<-rate/2
      }
      if(ok)break
    }
    history[[it]]<-data.frame(iteration=it,objective=q$value,gradient=max(abs(q$gradient)),
      residual=max(abs(s$T%*%q$state$U)),step_size=if(converged)0 else if(ok)rate else NA_real_,damping=damp)
    kdr_log("Final estimation",paste0(it,"/",control$maxit),"gradient",signif(max(abs(q$gradient)),3))
    if(converged)break
    if(!ok)stop("The regularized criterion line search failed.")
    b<-b+rate*step
  }
  if(!converged)stop("Final estimation reached its iteration limit.")
  list(beta=b,history=kdr_bind(history),converged=TRUE,objective=q$value,
    gradient=max(abs(q$gradient)),residual=max(abs(s$T%*%q$state$U)),iterations=it,
    gamma=Reduce(`+`,lapply(aux$models,function(a)a$g$gamma))/length(aux$models))
}
.kdr_curvature <- function(beta,z,aux,s,lambda) {
  a<-.kdr_state(beta,z,aux);N<-z$d$n;M<-length(beta);W<-crossprod(s$T)
  r0<-drop(W%*%a$U);JW<-t(a$J)%*%W;E<-matrix(0,M,M);sensitivity<-list()
  for(k in seq_along(a$parts)) {
    q<-a$parts[[k]];ix<-q$ix;DS<-z$pa$DS[ix,,drop=FALSE];RF<-z$RF[ix,,drop=FALSE]
    rr<-drop(q$R%*%r0);ae<-q$w*q$e/N
    E<-E+crossprod(q$Z,q$Z*(ae*rr))
    Gg<--crossprod(q$R,DS*ae)
    Ag<-JW%*%Gg+crossprod(q$Z,DS*(ae*rr))
    v<-colSums(RF*(q$w*q$t))/N
    U<-crossprod(q$Z,RF*ae)
    sensitivity[[k]]<-list(Ag=Ag,v=v,U=U)
  }
  H<-kdr_sym(crossprod(s$T%*%a$J)+E+lambda*s$penalty)
  list(H=H,state=a,r0=r0,JW=JW,sensitivity=sensitivity)
}
.kdr_covariance <- function(beta,z,aux,s,lambda) {
  a<-.kdr_curvature(beta,z,aux,s,lambda);N<-z$d$n;M<-length(beta)
  ev<-eigen(a$H,symmetric=TRUE,only.values=TRUE)$values
  if(min(ev)<=0)stop("Exact local curvature is nonpositive.")
  Hi<-kdr_inverse(a$H);rows<-matrix(0,length(z$pa$D),M)
  for(q in a$state$parts)rows[q$ix,]<-(q$R*(q$w*q$t))%*%t(a$JW)-q$Z*(q$w*q$e*drop(q$R%*%a$r0))
  S<-group_sum(rows,z$pa$group);ff<-aux$folds[z$pa$group]
  for(k in seq_along(aux$models)) {
    ax<-aux$models[[k]];g<-ax$g
    eta<-drop(z$Traw%*%g$raw_beta);mx<-as.numeric(tapply(eta,z$d$group,max));ee<-exp(eta-mx[z$d$group])
    pr<-ee/drop(group_sum(ee,z$d$group))[z$d$group]
    IG<-group_sum((z$Traw%*%g$influence_map)*(z$d$event-pr),z$d$group)
    IG[aux$folds==k,]<-0
    S<-S+(N/g$ntrain)*IG%*%t(a$sensitivity[[k]]$Ag)
    ix<-which(ff!=k & z$pa$D==0)
    RF<-z$RF[ix,,drop=FALSE];R<-z$pa$Z[ix,,drop=FALSE]-RF%*%t(ax$f$Pi)
    Q<-ax$f$gram_inverse;v<-a$sensitivity[[k]]$v;U<-a$sensitivity[[k]]$U
    # Contraction of the full finite-dimensional f sensitivity with its influence.
    h<-drop(Q%*%v);B<-U%*%Q;w<-z$pa$w[ix]
    Frows<--(R%*%t(a$JW))*(w*drop(RF%*%h))+(RF%*%t(B))*(w*drop(R%*%a$r0))
    FS<-group_sum(Frows,z$pa$group[ix]);ids<-sort(unique(z$pa$group[ix]))
    S[ids,]<-S[ids,,drop=FALSE]+(N/ax$f$ntrain)*FS
  }
  centered<-sweep(S,2,colMeans(S));meat<-crossprod(centered)/(N*(N-1))
  V<-kdr_sym(Hi%*%meat%*%t(Hi));kdr_vroot(V)
  list(matrix=V,status="available",target="Regularized population target conditional on independent learning and selection",
    hessian=a$H,minimum_curvature=min(ev),condition=kappa(a$H,exact=TRUE))
}

.kdr_learn <- function(d,adjust,nu,control) {
  a<-.kdr_formula_fit(adjust,d$data);X<-a$X
  keep<-which(colnames(X)!="(Intercept)" & colSums(kdr_center(X,d$group)^2)>1e-20)
  if(length(keep)){q<-qr(kdr_center(X[,keep,drop=FALSE],d$group),tol=1e-12);keep<-sort(keep[q$pivot[seq_len(q$rank)]])}
  d$S<-X[,keep,drop=FALSE]
  cfg<-list(M=control$centers,kernel=list(coverage_fraction=control$coverage,kmeans_starts=2L,
    kmeans_maxit=25L,initial_scale=control$initial_scale,rho_bounds=control$rho_bounds,kappa=control$kappa,
    start_factors=1,maxit=control$kernel_maxit,continuations=control$kernel_continuations,
    projected_gradient_tol=control$kernel_tolerance,penalty_rows=control$penalty_rows),
    solver=list(maxit=control$maxit,tolerance=control$tolerance,step_tolerance=control$step_tolerance,max_backtrack=20L))
  s<-kdr_geometry(d,cfg,control$seed);s$learning<-kdr_learning_fit(d,s,cfg);s$rho<-s$learning$rho
  s$B<-kdr_kernel_whiten(kdr_kernel_raw(d$A,s),d$group);s$dimension<-ncol(s$B)
  kdr_log("Kernel basis",s$dimension,"of",s$M,"candidate centers retained")
  s$penalty<-kdr_kernel_penalty(d,s,cfg,control$seed+1L)
  s$design<-.kdr_design_fit(d,s,adjust,nu)
  z<-.kdr_prepare(d,s);C<-s$design$g_C
  gf<-kdr_cl_fit(function(FUN)FUN(list(V=z$V,event=d$event,group=d$group)),ncol(z$V),cfg)
  ix<-which(z$pa$D==0);w<-z$pa$w[ix];RF<-z$RF[ix,,drop=FALSE]
  Q<-crossprod(RF,RF*w)/d$n;Pi<-t(kdr_inverse(Q)%*%crossprod(RF,z$pa$Z[ix,,drop=FALSE]*w)/d$n)
  ax<-list(folds=rep(1L,d$n),models=list(list(g=list(gamma=drop(C%*%gf$beta)),f=list(Pi=Pi))))
  # For the common learning adjustment design, use its profile-likelihood coefficient.
  if(ncol(d$S)==nrow(C))ax$models[[1]]$g$gamma<-tail(s$learning$beta,ncol(d$S))
  s$T<-kdr_inverse(.kdr_state(numeric(s$dimension),z,ax)$J)
  s$learning_N<-d$n;s$exposures<-colnames(d$A);s$adjust<-adjust;s$nuisance<-nu
  n<-min(control$evaluation_points,nrow(d$A));s$evaluation<-as.data.frame(d$A[seq_len(n),,drop=FALSE])
  s$reference<-s$center;s$kappa<-control$kappa
  class(s)<-c("kdrcc_kernel","list");s
}

kdrcc <- function(data,case="case",strata="set_id",exposures,adjust=~1,
    learning_data=NULL,nuisance=kdrcc_nuisance(),control=kdrcc_control(),
    reference=NULL,kernel=NULL,fold_id=NULL) {
  call<-match.call();started<-proc.time()[3]
  if(!inherits(control,"kdrcc_control")||!inherits(nuisance,"kdrcc_nuisance"))stop("Use kdrcc_control() and kdrcc_nuisance() for configuration.")
  old<-options(kdrcc.trace=control$verbose);on.exit(options(old),add=TRUE)
  if(is.null(kernel)&&is.null(learning_data)) {
    parts<-kdrcc_split(data,strata,learning=.25,seed=control$seed)
    learning_data<-parts$learning;data<-parts$estimation;assignment<-parts$assignment
  } else assignment<-NULL
  d<-.kdr_validate(data,case,strata,exposures,adjust)
  if(d$n<=control$folds)stop("The estimation sample requires more matched sets than cross-fitting folds.")
  if(is.null(kernel)) {
    ld<-.kdr_validate(learning_data,case,strata,exposures,adjust)
    if(identical(d$A,ld$A)&&identical(d$event,ld$event))stop("Learning and estimation samples must use distinct matched sets.")
    kdr_log("Independent kernel learning",ld$n,"matched sets")
    kernel<-.kdr_learn(ld,adjust,nuisance,control)
  } else {
    if(!inherits(kernel,"kdrcc_kernel")||!identical(kernel$exposures,exposures))stop("The learned kernel exposure specification differs from exposures.")
    if(!identical(deparse(kernel$adjust),deparse(adjust)))stop("The learned adjustment specification differs from adjust.")
    if(!missing(nuisance)&&!identical(nuisance,kernel$nuisance))stop("The learned kernel carries its auxiliary design; use its nuisance specification.")
    control$centers<-nrow(kernel$centers);control$kappa<-kernel$kappa
  }
  if(is.null(fold_id))fold_id<-with_seed(control$seed+2L,sample(rep(seq_len(control$folds),length.out=d$n)))
  if(length(fold_id)!=d$n||anyNA(fold_id)||!setequal(unique(fold_id),seq_len(control$folds)))stop("fold_id must provide all cross-fitting fold labels in matched-set order.")
  z<-.kdr_prepare(d,kernel);aux<-.kdr_auxiliary(z,kernel,fold_id,control)
  point<-.kdr_solve(z,aux,kernel,control)
  covar<-tryCatch(.kdr_covariance(point$beta,z,aux,kernel,control$lambda),error=function(e)
    list(matrix=NULL,status="unavailable",error=conditionMessage(e),target="Regularized population target conditional on independent learning and selection"))
  if(covar$status=="unavailable")warning("Covariance calculation: ",covar$error,call.=FALSE)
  if(is.null(reference))reference<-kernel$reference
  if(!is.null(names(reference)))reference<-reference[exposures]
  if(length(reference)!=length(exposures)||any(!is.finite(reference)))stop("reference must give one finite value per exposure.")
  names(point$beta)<-paste0("kernel",seq_along(point$beta));names(reference)<-exposures
  out<-list(call=call,version=list(package="0.1.0",object="1.0"),
    specification=list(case=case,strata=strata,exposures=exposures,adjust=adjust,nuisance=kernel$nuisance,control=control),
    data_summary=list(N=d$n,L=kernel$learning_N,dates=nrow(d$A),matched_set_sizes=table(d$counts),
      adjustment_columns=kernel$design$adjust_columns,candidate_centers=nrow(kernel$centers),
      kernel_dimension=kernel$dimension),
    reference=reference,kernel=kernel,nuisance=aux,coefficients=point$beta,covariance=covar,
    tuning=kernel$tuning_record,diagnostics=list(estimation=point$history,kernel=kernel$learning$diagnostics,
      objective=point$objective,gradient=point$gradient,residual=point$residual,converged=point$converged,
      rho_lower=sum(kernel$rho<=kernel$base_rho*control$rho_bounds[1]*(1+1e-6)),
      rho_upper=sum(kernel$rho>=kernel$base_rho*control$rho_bounds[2]*(1-1e-6)),
      elapsed_seconds=unname(proc.time()[3]-started)),evaluation=kernel$evaluation,
    reproducibility=list(seed=control$seed,fold_id=fold_id,strata=d$ids,assignment=assignment,R=as.character(getRversion()),platform=R.version$platform),
    data=if(control$keep_data)data else NULL)
  class(out)<-"kdrcc";out
}
