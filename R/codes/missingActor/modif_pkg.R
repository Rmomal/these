generator_PLN<-function(Sigma,covariates=NULL, n=50){
  # ajout d'une constante, par rapport à EMtree::generator_PLN
  p<-ncol(Sigma)
  if(!is.null(covariates)){
    n<-nrow(covariates)
    c<-ncol(covariates)
    string<-paste0("~", paste(colnames(covariates), collapse=" + "))
    formula<-as.formula(string)
    m<- model.matrix(formula,covariates)[,-1]
    mc<-ncol(m)
    beta<-matrix(runif(p*mc),mc,p)
    prod=m %*% beta+2
  }else{
    prod=2
  }
  D<-diag(Sigma)
  R<-cov2cor(Sigma)
  U<- rmvnorm(n, rep(0,p), R)
  matsig=(matrix(rep(sqrt(D),n),n,p, byrow = TRUE))
  Y = matrix(rpois(n*p, exp(U*matsig+prod )), n, p)
 
  return(list(Y=Y, U=U))
}
generator_param<-function(G,signed=FALSE,v=0){
  lambda = 1;  p=ncol(G);  sumlignes=rowSums(matrix(G,p,p));  D=diag(sumlignes+v)
  if(signed){
    Gsign = F_Vec2Sym(F_Sym2Vec(G * matrix(2*rbinom(p^2, 1, .3)-1, p, p)))
    omega = lambda*D - Gsign
    while(min(eigen(omega)$values) < 1e-10 & lambda<1e3){
      lambda = 1.1*lambda
      omega = lambda*D - Gsign
    }
  }else{
    omega = lambda*D + G
    while (min(eigen(omega)$values) < 1e-10){
      lambda = 1.1*lambda
      omega =lambda*D + G
    }
  }
  sigma = (solve(omega))
  sim=list(sigma=sigma,omega=omega,cste=lambda)
  return(sim)
}

library(huge)
library(Matrix)
#qu'il se taise
generator_graph<-function(p = 20, graph = "tree", prob = 0.1, dens=0.3, r=5){
  theta = matrix(0, p, p)
  if (graph == "cluster") {
    theta<-SimCluster(p,3,dens,r)
  }
  if (graph == "scale-free") {
    theta = huge.generator(d=p,graph="scale-free", verbose=FALSE)$theta
    
  }
  if(graph=="tree"){
    theta<-SpannTree(p)
  }
  if(graph=="erdos"){
    theta<- erdos(p=p,prob=prob)
  }
  return(theta = Matrix(theta, sparse = TRUE))
}
data_from_scratch<-function(type, p=20,n=50, ratio=5, covariates=NULL, prob=log(p)/p,
                            dens=log(p)/p, signed=FALSE,v=0,draw=FALSE){
  graph<- generator_graph(graph=type,p=p,prob=prob,dens=dens,r=ratio)
  param<-generator_param(G=as.matrix(graph),signed=signed,v=v)
  data<-generator_PLN(param$sigma,covariates,n)
  Y=data$Y
  U=data$U
  if(draw){
    g=as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="nicely")+
      geom_edge_link()+
      geom_node_point(size=2, color="blue")
    print(g)
  }
  return(list(Y=Y,U=U,omega= param$omega))
}
F_Vec2Sym <- function(A.vec){
  # Makes a symmetric matrix from the vector made of its lower tirangular part
  # the diagonal is nul
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
EMtree_corZ<-function(CovZ,n,  maxIter=30, cond.tol=1e-10, verbatim=FALSE, plot=FALSE){
  CorZ=cov2cor(CovZ)
  p = ncol(CorZ)
  alpha.psi = Psi_alpha(CorZ, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, maxIter = maxIter,
                        verbatim=verbatim, plot=plot)
  return(FitEM)
}

computeAlpha<-function(omegainitO,default=0.6, MO, SO, plot=TRUE){
  n=nrow(MO)
  alpha.grid=(1:n)/n
  phi=CorOmegaMatrix(omegainitO)
  vec.cond<-sapply(alpha.grid ,function(alpha){
    expr=F_Sym2Vec(alpha*(n*0.5*log(phi)-0.5*omegainitO*(t(MO)%*%MO)))
    expr=F_Vec2Sym(expr)+diag(colSums(SO))
    expr=exp(expr)
    expr=(expr-mean(expr))
    lambda = svd(expr)$d 
    cond = min(abs(lambda))/max(abs(lambda))
  })
  alpha=alpha.grid[max(which(vec.cond>1e-10))]
  if(is.na(alpha)){ alpha=default
  message("default alpha")}
  datacond<-data.frame(alpha=alpha.grid,test=log(vec.cond))
  if(plot){
    g<-datacond %>% gather(key, value, -alpha) %>% ggplot(aes(alpha,value, color=key))+geom_point()+mytheme+
      geom_hline(yintercept = log(1e-10))+guides(color=FALSE)+labs(y="Conditioning", title=round(alpha, 3))
    print(g)
  }
  return(alpha)
}