# 08/09
# Tests to set illegal links
library(EMtree)
library(PLNmodels)
#----------------
#-- Functions needed but not exported yet:
F_Vec2Sym <- function(A.vec){
  # Makes a symmetric matrix from the vector made of its lower tirangular part
  n = (1+sqrt(1+8*length(A.vec)))/2
  A.mat = matrix(0, n, n)
  A.mat[lower.tri(A.mat)] = A.vec
  A.mat = A.mat + t(A.mat)
  return(A.mat)
}
F_Sym2Vec <- function(A.mat){
  # Makes a vector from the lower triangular par of a symmetric matrix
  return(A.mat[lower.tri(A.mat)])
}

#-- modified gradient and likelihood functions to handle null weights:
F_NegLikelihood <- function(beta.vec, log.psi, P){
  M = Meila(F_Vec2Sym(beta.vec))
  lambda = SetLambda(P, M)
  return(- sum(F_Sym2Vec(P)*(log(beta.vec+(beta.vec==0))+F_Sym2Vec(log.psi))) +
           log(SumTree(F_Vec2Sym(beta.vec)))+
           lambda*(sum(beta.vec)-0.5))
}
 
F_NegGradient_Trans <- function(gamma, log.psi, P){
  beta=exp(gamma)
  beta[gamma==0]=0
  M = Meila(F_Vec2Sym(beta))
  lambda = SetLambda(P, M)
  return(- F_Sym2Vec(P)+ beta*(F_Sym2Vec(M) + lambda))
}
 
#-- new EMtree function :
new_EMtree<-function(PLN.Cor,unlinked=NULL, maxIter=30, cond.tol=1e-10, eps1 = 1e-6,eps2=1e-4, plot=FALSE){
  CorY=cov2cor(PLN.Cor$model_par$Sigma)
  alpha.psi = Psi_alpha(CorY, n, cond.tol=1e-10)
  psi = alpha.psi$psi
  # new beta.init with the unliked parameter
  beta.init = matrix(1, p, p); diag(beta.init) = 0; 
  if(!is.null(unlinked)) beta.init[unlinked, unlinked]=0
  beta.init = beta.init / sum(beta.init)
  
  # new FitBeta function
  beta.tol = 1e-4
  beta.min = 1e-16
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi+(psi==0))
  iter = 0
  logpY = rep(0, maxIter)
  beta.diff = diff.loglik = 2 * eps2
  T1<-Sys.time()
  stop=FALSE
  while (((beta.diff > eps1) || (diff.loglik>eps2) ) && iter < maxIter && !stop){
    iter = iter+1
    P = EdgeProba(beta.old*psi)
    init=F_Sym2Vec(beta.old)
    long=length(F_Sym2Vec(beta.old))
    #modification for 0 in gamma_init
    gamma_init=log(init+(init==0))
    gamma_init[init==0]=0
    gamma = stats::optim(gamma_init, F_NegLikelihood_Trans, gr=F_NegGradient_Trans,method='BFGS',
                         log.psi, P)$par
    beta=exp(gamma)
    beta[which(beta< beta.min)] = beta.min
    beta=F_Vec2Sym(beta)
    logpY[iter] = -F_NegLikelihood(F_Sym2Vec(beta),log.psi,P)
    beta.diff = max(abs(beta.old-beta))
    beta.old = beta
    diff.loglik=logpY[iter]-logpY[iter-1]
    if(iter > 1){diff.loglik =  abs(diff.loglik)}else{diff.loglik=1}
  }
  time<-difftime(Sys.time(),T1)
  logpY = logpY[1:iter]
  if(plot){
    g<-tibble(p=logpY) %>% rowid_to_column() %>%
      ggplot(aes(rowid,p))+geom_point()+geom_line()+theme_minimal()+labs(x="Iter",y="Likelihood")
    print(g)
  }
  P = EdgeProba(beta*psi)
  cat("\nConvergence took",round(time,2), attr(time, "units")," and ", iter," iterations.")
  return(list(edges_prob=P, edges_weight=beta, logpY=logpY,maxIter=iter, norm.cst = SumTree(beta), timeEM=time))
}

#-- new Resampling
new_ResampleEMtree <- function(counts, unlinked=NULL, covar_matrix=NULL  , O=NULL, v=0.8, S=1e2, maxIter=30, cond.tol=1e-14,cores=3){
  cat("Computing ",S,"probability matrices with", cores, "core(s)... ")
  t1=Sys.time()
  counts=as.matrix(counts)
  n = nrow(counts);  p = ncol(counts)
  P = p * (p - 1) / 2 ; V = round(v * n)
  Pmat = matrix(0, S, P)
  #- offsets and covariates
  if(is.null(O)){ O=matrix(1, n, p)}
  if(is.null(covar_matrix)){#default intercept
    X=matrix(1,nrow=n,ncol=1)
  }else{X=as.matrix(covar_matrix)}
  #- parallel computation of S fits of new_EMtree
  obj<-parallel::mclapply(1:S,function(b){
    set.seed(b)
    sample = sample(1:n, V, replace = F)
    counts.sample = counts[sample,]
    X.sample = data.frame(X[sample,])
    O.sample = O[sample,]
    try({
      suppressWarnings(
      PLN.sample <- PLNmodels::PLN(counts.sample ~ -1  + offset(log(O.sample)) + ., data=X.sample, control = list("trace"=0))
    )
      inf<-new_EMtree( PLN.sample,unlinked, maxIter=maxIter, cond.tol=cond.tol,
                    plot=FALSE)[c("edges_prob","maxIter","timeEM")]
    },silent=TRUE )
    if(!exists("inf")) inf=NA #depending on the sample drawn, it is possible that computation fail
    # because of bad conditioning of the Laplacian matrix of the weights beta.
    # This is the price for setting some weights to zero.
    return(inf)
  }, mc.cores=cores)
  bad_samples=which(do.call(rbind, lapply(obj, length))!=3)
  time=difftime(Sys.time(), t1)
  cat(round(time,2),  attr(time, "units"))
  if(length(bad_samples)!=0){
    cat("\n",length(bad_samples), " failed samples.\n")
    obj=obj[-bad_samples]
  } 
  Pmat<-do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$edges_prob)}))
  summaryiter = do.call(c,lapply(obj,function(x){x$maxIter}))
  times<-do.call(c,lapply(obj,function(x){x$timeEM}))
  return(list(Pmat=Pmat,maxIter=summaryiter,times=times))
}