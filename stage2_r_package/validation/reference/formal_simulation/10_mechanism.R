kdr_anchor <- function(x,a) {
  f<-splinefun(a$x,a$percent,method="monoH.FC");lo<-min(a$x);hi<-max(a$x);tau<-(hi-lo)/3
  y<-f(pmax(lo,pmin(hi,x)));l<-x<lo;u<-x>hi
  y[l]<-f(lo)+f(lo,deriv=1)*tau*(exp((x[l]-lo)/tau)-1)
  y[u]<-f(hi)+f(hi,deriv=1)*tau*(1-exp(-(x[u]-hi)/tau))
  if(any(y<=-100))stop("Invalid response-surface tail")
  log1p(y/100)
}
kdr_truth <- function(A,cfg,scenario="baseline") {
  t<-cfg$baseline$truth
  H<-vapply(1:6,function(j)kdr_anchor(A[,j],t$anchors[[j]]),numeric(nrow(A)))
  tt<-(A[,7]-t$temperature$optimum)/t$temperature$distance
  value<-drop(H%*%t$pollutant_coefficients)+t$temperature$cold_log_coefficient*pmin(tt,0)^2+
    t$temperature$hot_log_coefficient*pmax(tt,0)^2
  Z<-sweep(sweep(A[,1:8,drop=FALSE],2,cfg$baseline$reference$center),2,cfg$baseline$reference$scale,"/")
  for(i in seq_len(nrow(t$interactions))) {
    a<-t$interactions[i,];value<-value+a$coefficient*Z[,a$q]*Z[,a$r]
  }
  value
}
kdr_generate <- function(N,p,cfg,seed,events=TRUE) {
  d<-generate_candidates(N,cfg$baseline$dgp,seed)
  for(j in 1:6)d[[exposure_names[j]]]<-exp(cfg$baseline$mapping[j,1]+cfg$baseline$mapping[j,2]*log(d[[exposure_names[j]]]))
  A<-as.matrix(d[,exposure_names]);g<-d$group
  if(!is.null(cfg$baseline$dgp$copula_eta))A<-sim_D_apply_copula(A,g,d$date,cfg,seed)
  h<-kdr_truth(A,cfg);eta<-h+outcome_confounding(d$W,cfg$baseline$dgp)
  mx<-as.numeric(tapply(eta,g,max));ep<-exp(eta-mx[g]);prob<-ep/drop(group_sum(ep,g))[g]
  y<-integer(length(g))
  if(events) {
    counts<-tabulate(g);ends<-cumsum(counts);starts<-c(0L,head(ends,-1L));cs<-cumsum(prob)
    local<-cs-c(0,cs[head(ends,-1L)])[g]
    u<-with_seed(stage_seed(seed,"case_date"),runif(N))
    choice<-pmin(drop(group_sum(as.integer(local<u[g]),g)),counts-1L)+1L
    y[starts+choice]<-1L
    stopifnot(all(drop(group_sum(y,g))==1L))
  }
  ans<-list(A=A,W=d$W,group=g,date=d$date,event=y,truth=h,probability=prob,n=as.integer(N),
    fold=with_seed(stage_seed(seed,"fold"),sample(rep(seq_len(cfg$nuisance$folds),length.out=N))))
  sim_add_null(ans,p,cfg,seed)
}
sim_add_null <- function(d,p,cfg,seed) {
  stopifnot(p>=8L)
  if(p>8L) {
    E<-vapply(seq_len(p-8L),function(j) with_seed(stage_seed(seed,"independent_null",j),{
      u<-rnorm(nrow(d$A));v<-rnorm(d$n)
      sqrt(1-cfg$noise$group_fraction)*u+sqrt(cfg$noise$group_fraction)*v[d$group]
    }),numeric(nrow(d$A)))
    d$A<-cbind(d$A[,1:8,drop=FALSE],matrix(E,nrow=nrow(d$A)))
  } else d$A<-d$A[,1:8,drop=FALSE]
  colnames(d$A)<-kdr_names(p);d
}
sim_subset <- function(d,n) {
  if(n==d$n)return(d)
  stopifnot(n>=1,n<=d$n);take<-d$group<=n
  for(k in c("A","W","date","event","truth","probability"))
    d[[k]]<-if(is.matrix(d[[k]]))d[[k]][take,,drop=FALSE] else d[[k]][take]
  d$group<-d$group[take];d$fold<-d$fold[seq_len(n)];d$n<-as.integer(n);d$pairs<-NULL;d
}
sim_pool_file <- function(cache,role,r,b)file.path(cache,"raw",paste0(role,"_r",r,"_b",b,".rds"))
sim_raw <- function(cfg,cache,N,p,r,role="estimation") {
  memo<-new.env(parent=emptyenv());used<-0
  limit<-if(role=="estimation")cfg$execution$raw_memory_MB*1024^2 else 0
  if(!length(limit))limit<-0
  function(FUN,progress=FALSE) {
    bs<-cfg$execution$block_groups;nb<-ceiling(N/bs)
    for(b in seq_len(nb)) {
      key<-as.character(b)
      if(exists(key,memo,inherits=FALSE))d<-memo[[key]] else {
        d<-readRDS(sim_pool_file(cache,role,r,b))
        seed<-kdr_seed(cfg,"data_pool",role,r,b)
        d<-sim_add_null(d,p,cfg,seed)
        d<-sim_subset(d,min(bs,N-(b-1L)*bs))
        d$offset<-(b-1L)*bs;d$pairs<-kdr_pairs(d)
        bytes<-as.numeric(object.size(d))
        if(used+bytes<=limit){memo[[key]]<-d;used<<-used+bytes}
      }
      FUN(d)
      if(progress)kdr_bar(paste("Data blocks:",role),b,nb)
      kdr_heartbeat(paste("Reading",role,"block",b,"of",nb))
    }
    invisible(NULL)
  }
}
sim_learning_data <- function(cfg,cache,L,p,r) {
  parts<-list();sim_raw(cfg,cache,L,p,r,"learning")(function(d)parts[[length(parts)+1L]]<<-d)
  ans<-list()
  for(n in c("A","W","date","event","truth","probability","fold"))
    ans[[n]]<-if(n=="A")do.call(rbind,lapply(parts,`[[`,n)) else do.call(c,lapply(parts,`[[`,n))
  offsets<-c(0L,head(cumsum(vapply(parts,function(d)d$n,integer(1))),-1L))
  ans$group<-unlist(Map(function(d,o)d$group+o,parts,offsets),use.names=FALSE)
  ans$n<-as.integer(sum(vapply(parts,function(d)d$n,integer(1))));ans$pairs<-kdr_pairs(ans);ans
}
sim_evaluation <- function(cfg,p,path) {
  old<-sim_read(path);if(!is.null(old))return(old)
  d<-kdr_generate(cfg$evaluation$groups,p,cfg,kdr_seed(cfg,"evaluation"),events=FALSE)
  ref<-apply(d$A,2,median);qq<-cfg$evaluation$quantiles
  joint<-t(vapply(qq,function(q)apply(d$A,2,quantile,q),numeric(p)))
  rows<-match(seq_len(cfg$evaluation$function_points),d$group);A<-d$A[rows,,drop=FALSE]
  truth0<-kdr_truth(matrix(ref,1,p),cfg)
  out<-list(reference=ref,quantile=qq,joint=joint,individual=A,
    true_joint=kdr_truth(joint,cfg)-truth0,true_individual=kdr_truth(A,cfg)-truth0,
    coverage_ids=seq_len(cfg$evaluation$coverage_points))
  out$interpretation<-sim_interpretation_context(out,cfg)
  kdr_save(out,path);out
}
