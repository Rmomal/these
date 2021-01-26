# 13/10
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
F_NegLikelihood <- function(beta.vec, log.psi, P,sum.constraint){
  M = Meila(F_Vec2Sym(beta.vec))
  lambda = SetLambda(P, M,sum.constraint)
  return(- sum(F_Sym2Vec(P)*(log(beta.vec+(beta.vec==0))+F_Sym2Vec(log.psi))) +
           log(SumTree(F_Vec2Sym(beta.vec)))+
           lambda*(sum(beta.vec)-sum.constraint/2))
}

F_NegGradient_Trans <- function(gamma, log.psi, P,sum.constraint){
  beta=exp(gamma)
  beta[gamma==0]=0
  M = Meila(F_Vec2Sym(beta))
  lambda = SetLambda(P, M,sum.constraint)
  return(- F_Sym2Vec(P)/beta + (F_Sym2Vec(M) + lambda))
}
F_NegLikelihood_Trans <- function(gamma, log.psi, P,sum.constraint, trim=TRUE){
  if(trim){
    gamma=gamma-mean(gamma)
    gamma[which(gamma<(-30))]=-30
    gamma[which(gamma>(40))]=40
  }
  M = Meila(F_Vec2Sym(exp(gamma)))
  lambda = SetLambda(P, M,sum.constraint)

  suppressWarnings(
    res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+
      log(SumTree(F_Vec2Sym(exp(gamma))))+
      lambda*(sum(exp(gamma))-sum.constraint/2))

  if(is.nan(res)){
    #  cat(max(gamma),": higher bound ")
    gamma[which(gamma>(20))]=20
    gamma[which(gamma<(-20))]=-20
    M = Meila(F_Vec2Sym(exp(gamma)))
    lambda = SetLambda(P, M,sum.constraint)
    res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+
      log(SumTree(F_Vec2Sym(exp(gamma))))+
      lambda*(sum(exp(gamma))-sum.constraint/2)
    if(is.nan(res))   cat("\nbeta optimization failed: ",SumTree(F_Vec2Sym(exp(gamma))),
                          sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) ,
                          sum(F_Sym2Vec(P)),"\n")
  }

  return( res)
}
binf.constraint<-function(p,min.order=308){
  round( p*(p-1)*10^(-min.order/(p-1)))+1
}
SetLambda <- function(P, M,sum.constraint=1, eps = 1e-6, start=1){
  # F.x has to be increasing. The target value is 0
  F.x <- function(x){
    if(x!=0){
      sum.constraint - sum(P / (x+M))
    }else{
      sum.constraint - (2*sum(P[upper.tri(P)] / M[upper.tri(M)]))
    }
  }
  x.min = ifelse(F.x(0) >0,-start,1e-4);
  t1<-Sys.time()
  while(F.x(x.min)>0 ){
    x.min = x.min -1
    if(difftime(Sys.time(),t1)>1) stop("Could not set lambda.")}
  x.max = start
  while(F.x(x.max)<0){x.max = x.max * 2}
  x = (x.max+x.min)/2
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  while(abs(x.max-x.min) > eps){
    if(f > 0) {
      x.max = x
      f.max = f
    } else{
      x.min = x
      f.min = f
    }
    x = (x.max+x.min)/2;
    f = F.x(x)
  }
  return(x)
}
binf.constraint<-function(p,min.order=308){
  round(p*(p-1)*(10^(-min.order/(p-1))))+1
}
#-- new EMtree function :
new_EMtree<-function(PLN.Cor,unlinked=NULL,n=NULL,unif=TRUE, maxIter=30, cond.tol=1e-10,
                     eps1 = 1e-6,eps2=1e-4, plot=FALSE){
  if(inherits(PLN.Cor, "PLNfit")){
    CorY=cov2cor(PLN.Cor$model_par$Sigma)
    n=PLN.Cor$n
  }else if(inherits(PLN.Cor, "matrix") & nrow(PLN.Cor) == ncol(PLN.Cor)){
    CorY = PLN.Cor
  }else{
    stop("PLN.Cor must be a PLN object or a squarred gaussian correlation matrix")
  }
  p=ncol(CorY)
  # set the tempering parameter alpha, for a good conditioning of the psi matrix
  alpha.psi = Psi_alpha(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi

  # beta.init with adaptative mean value depending on the network dimensions
  sum.weights=binf.constraint(p)
  mean.val=sum.weights/(p*(p-1))
  beta.init = matrix(mean.val, p, p);
  if(!is.null(unlinked)) beta.init[unlinked, unlinked]=0 #unliked nodes

  if(!unif){ # try different starting points for this EM
    beta.init = matrix(runif(n=p*p, min=0.9*mean.val,max=1.1*mean.val ), p,p)
    beta.init=t(beta.init)%*%beta.init/2}
  diag(beta.init)=0

  # new FitBeta function
  beta.tol = 1e-4
  beta.min = 1e-16
  beta.old = beta.init# / sum(beta.init)
  log.psi = log(psi+(psi==0))
  iter = 0
  logpY = rep(0, maxIter)
  beta.diff = diff.loglik = 2 * eps2
  T1<-Sys.time()
  stop=FALSE
  while (((beta.diff > eps1) || (diff.loglik>eps2) ) && iter < maxIter && !stop){
    iter = iter+1
    P=Kirshner(W=beta.old*psi)
    if(is.nan(sum(P))) browser()
    init=F_Sym2Vec(beta.old)

    #gradient ascent in log scale
    gamma_init=log(init+(init==0))
    gamma_init[init==0]=0
    gamma = stats::optim(gamma_init, F_NegLikelihood_Trans, gr=F_NegGradient_Trans,method='BFGS',
                         log.psi, P,sum.weights, control=list(maxit=300,  abstol=1e-5, reltol=1e-5))$par

    beta=exp(gamma)
    beta[which(beta< beta.min)] = beta.min # numerical zeros
    beta=F_Vec2Sym(beta)
    if(!is.null(unlinked)) beta[unlinked, unlinked]=0
    logpY[iter] = -F_NegLikelihood(F_Sym2Vec(beta),log.psi,P,sum.weights)
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
  P = Kirshner(beta*psi)
  cat("\nConvergence took",round(time,2), attr(time, "units")," and ", iter," iterations.")
  return(list(edges_prob=P, edges_weight=beta, logpY=logpY,maxIter=iter, norm.cst = SumTree(beta), timeEM=time))
}

#-- new Resampling
new_ResampleEMtree <- function(counts, unlinked=NULL, covar_matrix=NULL  ,
                               O=NULL, v=0.8, S=1e2, maxIter=30, cond.tol=1e-14,cores=3){
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
      inf<-new_EMtree( PLN.sample,unlinked,n=n, maxIter=maxIter, cond.tol=cond.tol,
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

Kirshner<-function(W){
  # W = squared weight matrix
  # Kirshner (07) formulas
  p = nrow(W)
  L = Laplacian(W)[-1,-1]
  #improve numerical capacity with gmp on error with solve
  Q=tryCatch({solve(L)},
             error=function(e){inverse.gmp(L)})
  Q = rbind(c(0, diag(Q)),
            cbind(diag(Q), (diag(Q)%o%rep(1, p-1) + rep(1, p-1)%o%diag(Q) - 2*Q)))
  Q = .5*(Q + t(Q))
  P = W * Q
  P = .5*(P + t(P))
  return(P)
}

