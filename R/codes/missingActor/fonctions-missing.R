#############
# initialization
initEM <- function(Sigma = NULL,n=1e6,cst=1.1,pca=TRUE,cliqueList) {
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Initialize Sigma and Omega by taking the principal component of cliques as initial value for the hidden variables.
  #   The PCA is done on the correaltion matrix, the variance of hidden variables are set to the empirical variances of the pca.
  #   The corresponding precision term in omega is set to 1 for identifiability reasons.
  # INPUT
  #   Sigma     :  variance-covariance matrix (p x p)
  #   cliqueList : list of found cliques, given as vectorss of nodes indices
  # OUTPUT
  #   Sigma0    : initial value of complete covariance matrix ((p+r)x(p+r) matrix)
  #   K0        : initial value of complete precision matrix ((p+r)x(p+r) matrix)
  #   clique    : vector containing the indices of the nodes in the clique
  # -----------------------------------------------------------------------------------------------------------------------------
  
  r=length(cliqueList) ; p=ncol(Sigma);H=(p+1):(p+r)
  # code sans simulation de données
  Corr <- cov2cor(Sigma); sigma <- sqrt(diag(Sigma))
  coef <- matrix(0, p, r) 
  sapply(seq_along(cliqueList), function(c){ 
    coef[, c] <<- rep(0, p); 
    if(length(cliqueList[[c]])>1){ # ne pas faire d'acp si la clique est un noeud
      pca <-eigen(cov2cor(Corr[cliqueList[[c]], cliqueList[[c]]]))
      coef[cliqueList[[c]], c] <<- pca$vectors[, 1]
    }else{
      coef[,c]<<-Corr[cliqueList[[c]],]
    }   
  }) 
  
  # Recontructing Sigma
  CorrFull <- rbind(cbind(Corr, Corr%*%coef), cbind(t(coef)%*%Corr, cst*t(coef)%*%Corr%*%coef))
  sigmaFull <- c(sigma,1/sqrt(diag(CorrFull)[(p+1):(p+r)])) 
  SigmaFull <- diag(sigmaFull) %*% CorrFull %*% diag(sigmaFull)
  
  # Initialising Omega
  OmegaFull <- tryCatch({solve(SigmaFull)},
                        error=function(e){browser()},finally={})
  # coefDiag <- c(rep(1, p), 1/sqrt(diag(OmegaFull)[p+(1:r)]))
  # OmegaFull <- diag(coefDiag) %*% OmegaFull %*% diag(coefDiag)
  if(r>1){
    OmegaFull[H,H]<-diag(diag(OmegaFull[H,H]))
  }
  
  return(list( Sigma0 = SigmaFull, K0 = OmegaFull, cliquelist = cliqueList))
}
#Find the clique created by the missing actor
findCliques<-function(covX,k){
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  # Find k groups of nodes which are 'clique likes"
  # INPUT
  #   covX   : variance-covariance matrix of data  (pxp)
  #   k      : number of groups
  #            is determined with hierarchical clustering.
  # OUTPUT
  #   cliquelist    : a list of groups of indices
  # -----------------------------------------------------------------------------------------------------------------------------
  
  Y = 1- abs(cov2cor(covX))
  res.kmeans<-kmeans(Y,k)
  # Compute the mean distance within each cluster (clique should have the smallest mean distances)
  scores<-res.kmeans$withinss / ((res.kmeans$size)*(res.kmeans$size-1)/2)
  # Renumbering of the cluster according the score (cluster one will have the smallest score)
  reordered.cluster<-order(scores)[res.kmeans$cluster]
  return(cliquelist = split(1:length(reordered.cluster),reordered.cluster))
}
FitSparsePCA <- function(Y, r=1, alphaGrid=10^(seq(-4, 0, by=.1))){
  # Fit sparse PCA on a grid of alpha
  # Needs an estimate o Sigma: empirical variance of the rotated scores 
  #    + diagonal risidual variances (why not?)
  n <- nrow(Y); p <- ncol(Y); alphaNb <- length(alphaGrid); nbDifferentH=-1
  while(nbDifferentH!=r){
    # Fits all sparsePCA
    sPCA <- list()
    for(a in 1:alphaNb){
      sPCA[[a]] <- spca(Y, k=r, alpha=alphaGrid[a], verbose=FALSE)
      sPCA[[a]]$Sigma <- cov(sPCA[[a]]$scores%*%t(sPCA[[a]]$transform))
      
      resVar <- (n-1)*apply(Y - sPCA[[a]]$scores %*% t(sPCA[[a]]$transform), 2, var)/n
      sPCA[[a]]$Sigma <- sPCA[[a]]$Sigma  + diag(resVar)
      sPCA[[a]]$df <- 1 + sum(sPCA[[a]]$loadings!=0)
      sPCA[[a]]$loglik <- sum(dmvnorm((Y), sigma=sPCA[[a]]$Sigma, log=TRUE))
      sPCA[[a]]$bic <- sPCA[[a]]$loglik - log(n)*sPCA[[a]]$df/2
    }
    #2 neighbors minimum
    df<-unlist(lapply(sPCA, function(sPca){sPca$df}))
    
    good<-do.call(rbind,lapply(sPCA, function(spca){
      vec_col<-apply(spca$loadings, 2,function(col){
        sum(col!=0)>1})
      return(sum(vec_col)==r)
    }))
    # Selects alpha via pseudo-BIC
    loglik <- unlist(lapply(sPCA, function(sPca){sPca$loglik}))
    bic <- unlist(lapply(sPCA, function(sPca){sPca$bic}))
    aOpt <- which.max(bic[good])
    # Find the cliques
    alphaOpt <- alphaGrid[aOpt]
    sPcaOpt <- sPCA[[aOpt]]
    sPcaOpt$loadings
    cliques <- list(); 
    sapply(1:ncol(sPcaOpt$loadings), function(j){
      cliques[[j]] <<- which(sPcaOpt$loadings[, j]!=0)
    })
    nbDifferentH<-length(unique(cliques))
  }
  return(list(sPcaOpt=sPcaOpt, alphaGrid=alphaGrid, alphaOpt=alphaOpt, 
              loglik=loglik, bic=bic, cliques=cliques))
}

boot_FitSparsePCA<-function(Y, B,r, cores=3){
  cliqueList<-mclapply(1:B, function(x){
    n=nrow(Y); v=0.8; n.sample=round(0.8*n, 0)
    ech=sample(1:n,n.sample,replace = FALSE)
    Y.sample=Y[ech,]
    c=FitSparsePCA(Y.sample,r=r)$cliques
    return(c)
  }, mc.cores=cores)
  
  nb_occ<-tabulate(match(cliqueList,unique(cliqueList)))
  cliqueList<-unique(cliqueList)
  return(list(cliqueList=cliqueList,nb_occ=nb_occ) )
}

init.mclust<-function(S,nb.missing=1, n.noise=50,plot=TRUE, title="",trueClique=NULL){
  Scomp=prcomp(S,scale. = TRUE)
  data=data.frame(Scomp$rotation[,1:2]%*%diag(Scomp$sdev[1:2]))
  datapolar=cart2pol(x=data[,1],y=data[,2])[,1:2] # recup r et theta
  datapolar_half=datapolar %>% mutate(theta2=ifelse(theta>pi,theta-pi,theta)) %>%
    dplyr::select(r,theta2)
  
  # noise in polar coords
  r <- sqrt(runif(n.noise))
  theta2 <- runif(n.noise, 0, pi)
  datapolarall=rbind(datapolar_half,cbind(r,theta2)) 
  newdata=pol2cart(datapolarall$r,datapolarall$theta2)[,1:2]#plot(newdata,col=col, xlim=c(-1,1),ylim=c(-1,1), pch=20)
  noiseInit<-sample(c(T,F), size=ncol(S), replace=T, prob=c(3, 1)) # noiseInit<-c(rep(F,ncol(S)),rep(T,n.noise))
  
  clust= tryCatch({
    Mclust(data=newdata, initialization = list(noise=noiseInit),  G=nb.missing, verbose = FALSE)
  }, error = function(e) {
    message("new noise")
    r <- sqrt(runif(n.noise))
    theta2 <- runif(n.noise, 0, pi)
    datapolarall=rbind(datapolar_half,cbind(r,theta2)) 
    newdata=pol2cart(datapolarall$r,datapolarall$theta2)[,1:2]
    Mclust(data=newdata, initialization = list(noise=noiseInit),  G=nb.missing, verbose = FALSE)
  }, finally = { })
  groups<-mclust::map(clust$z)[1:nrow(S)]
  res<-lapply(1:nb.missing, function(c){
    which(groups==c)
  })
  
  return(list(init=res, data=data))
}

computeFPN<-function(cliqueList, trueClique,p){
  N=setdiff(1:p,trueClique)
  FP=unlist(lapply(cliqueList, function(init){sum(init%in%N)/length(N)}))
  FN=unlist(lapply(cliqueList, function(init){
    sum(setdiff(1:p,init)%in%trueClique)/length(trueClique)}))
  return(data.frame(FP=FP, FN=FN))
}
plotInitMclust<-function(res, title){
  tmp=sapply(res$init, function(c){
    res=1*(1:p) %in% c
  })
  tmp=tmp %>% as_tibble() %>% mutate(null=ifelse(rowSums(tmp)==0,1,0))
  colors= as.matrix((tmp))%*%matrix(1:ncol(tmp), ncol(tmp), 1)
  g= ggplot(res$data,aes(X1,X2, label=1:p,color=as.factor(colors)))+geom_point(size=0.1)+
    theme_light()+geom_text()+labs(x="eig vect 1",y="eig vect 2", title=title)+
    guides(color=FALSE)+scale_color_brewer(palette="Dark2")+
    geom_hline(yintercept=0, color="gray50")+geom_vline(xintercept=0, color="gray50")
  print(g)
  
}
# Initialization 
missing_from_scratch<-function(n,p,r,type,plot){
  data=data_from_scratch(type = type,p = p+r,n = n,signed = FALSE,prob = 2/p,v = 0)
  omega=data$omega
  data=generator_PLN(solve(omega),covariates = NULL,n=n)
  hidden=which(diag(omega)%in%sort(diag(omega), decreasing = TRUE)[1])[1] # on cache les r plus gros
  trueClique=lapply(hidden, function(h){
    which(omega[h,-h]!=0)})  
  group=1*(diag(omega)==diag(omega)[hidden][1])
  if(plot){
    G=draw_network(1*(omega==1),groupes=group, 
                   layout="nicely",curv=0,nb=2,pal="black",nodes_label =1:(p+r))$G
    print(G)
  }
  if(r!=0){
    Kh  <- omega[hidden,hidden]
    Ko  <- omega[-hidden,-hidden]
    Koh <- omega[-hidden,hidden]
    Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
    sigmaO=solve(Km)
    counts=data$Y[,-hidden]
    ZH=data$Z[,hidden]
  }else{
    group=NULL
    counts=data$Y 
    ZH=NULL
    sigmaO=solve(omega)
  }
  return(list(Y=counts, ZH=ZH, Sigma=sigmaO, Omega=omega, TC=trueClique, H=hidden))
}


initVEM<-function(counts,initviasigma,sigma_obs,r){
  p=ncol(counts)
  n=nrow(counts)
  # Tree
  Wginit <- matrix(1, p+r, p+r); Wginit =Wginit / sum(Wginit)
  Winit <- matrix(1, p+r, p+r); Winit =Winit / sum(Winit)
  Winit[1:p,1:p] <- EMtree_corZ(cov2cor(sigma_obs),n = n,maxIter = 20,verbatim = FALSE)$edges_weight
  diag(Wginit) = 1;diag(Winit) = 1
  # Z
  if(r!=0){
    initial.param<-initEM(sigma_obs,n=n,cliqueList = (initviasigma),cst=1.05, pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
    omegainit=initial.param$K0
    MHinit<-sapply(initviasigma, function(clique){
      pr=prcomp(t(counts[,clique]),scale. = FALSE)
      res = matrix(pr$rotation[,1]*pr$sdev[1],nrow=n,ncol=1)
      return(res)
    })
  }else{#init with no missing actors
    omegainit=solve(sigma_obs)
    MHinit=NULL
  }
  return(list(Wginit= Wginit, Winit= Winit, omegainit=omegainit,MHinit=MHinit))
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
#################
# For OPTIM

#=====
# VE step

computeWg<-function(phi,omega,W,MH,MO,S, alpha,hidden=TRUE, hist=FALSE, verbatim=FALSE ){
  p=ncol(MO); q=ncol(omega); O = 1:p ; n=nrow(MO) ;  binf=exp(-20) ; bsup=exp(30)
  Wg<-matrix(0,ncol(omega),ncol(omega))
  logWg<-matrix(0,ncol(omega),ncol(omega))
  if(hidden){ r=ncol(MH)
  H = (p+1):q  
  M=cbind(MO,MH)
  }else{ r=0 
  M=MO}
  SigmaTilde= (t(M)%*%M + diag(colSums(S)))/n
  
  ## update
  logWg<-log(W+exp(-20)*(W==0))+alpha*(n*0.25*log(phi+exp(-20)*(phi==0))-0.5*n*(omega*SigmaTilde))
  diag(logWg) = 0
  #--- centrage O
  mean.gamma=F_Sym2Vec(logWg)
  gammaO=F_Sym2Vec(logWg[O,O])
  gammaO=gammaO-mean(gammaO) 
  Wg[O,O]=exp(F_Vec2Sym(gammaO))
  
  if(hidden){
    #--- centrage OH
    gammaOH=logWg[O,H]
    gammaOH=gammaOH-mean(gammaOH)
    rangeO=max(gammaO)-min(gammaO)
    rangeOH=max(gammaOH)-min(gammaOH)
    if(verbatim)  cat(paste0(", ratio rOH/rOO=", round(rangeOH/rangeO,5)))
    Wg[O,H]=exp(gammaOH)
    Wg[H,O]<-t(Wg[O,H])
  }
  #--- trimming
  
  if(hist){
    par(mfrow=c(2,1))
    hist(gammaO, breaks=20,main="O")
    hist(gammaOH, breaks=20,main="OH")
  }
  # q99=quantile(Wg[O,O],0.99)
  # Wg[Wg>min(bsup,q99)]<-min(bsup,q99)
  # Wg[Wg<binf]<-binf
  
  Wg[W==0]=0
  Wg[phi==0]=0
  if(hidden) Wg[H,H]=0
  diag(Wg)=0
  norm1=10^(round(log10(min(Wg[Wg!=0]))-(-15))) #division maximale pour ne pas tuer les poids les plus petits
  norm2=max(10^(round(log10(max(Wg))-1)),1)
  Wg=Wg/(min(norm1, norm2))
  return(list(Wg=Wg ))
}


CorOmegaMatrix<-function(omega){ 
  q=ncol(omega)
  CorOmega = matrix(0,q,q)
  sapply(1:q-1,
         function(j){
           sapply((j+1):q, #visits all combinations of observed nodes
                  function(k){
                    CorOmega[k, j] <<- 1-(omega[k,j]^2/(omega[k,k]*omega[j,j]))
                    CorOmega[j, k] <<- CorOmega[k, j]
                  }
           )
         })
  CorOmega[abs(CorOmega)>0 & abs(CorOmega)<1e-15]<-1e-15
  return(CorOmega)
}

#=====
# M step


computeOffDiag<-function(omegaDiag,SigmaTilde,p){ 
  q=length(omegaDiag)   
  hidden=(q!=p)
  omega=matrix(0,q,q)
  sapply(1:(q-1),
         function(j){
           sapply((j+1):q,
                  function(k){
                    omega[k, j] <<- (1 - ( 1+4*SigmaTilde[k,j]^2*omegaDiag[j]*omegaDiag[k])^0.5)/(2*SigmaTilde[k,j])
                    omega[j, k] <<- omega[k, j]
                  }
           )
         }
  )
  diag(omega)=omegaDiag
  if(hidden){
    H=(p+1):q
    if(length(H)>1){
      omega[H,H]<-diag(diag(omega[H,H]))
    }
  }
  return(omega)
}

omegaDiag<-function(Pg,omega,SigmaTilde){ #closed form
  quantity<-Pg*omega*SigmaTilde
  quantity[lower.tri(quantity,diag=FALSE)]<-NA
  vectorSuml=colSums(quantity, na.rm=TRUE)
  omegaDiag = (1-vectorSuml)/diag(SigmaTilde)
  return(omegaDiag)
}
#=====
#Lower bound

LowerBound<-function(Pg ,omega, M, S, W, Wg,p, logSTW, logSTWg){
  n=nrow(M) ; q=nrow(omega) ; O=1:p ; r=q-p
  phi=CorOmegaMatrix(omega)
  hidden = (q!=p)
  if(hidden) H=(p+1):q
  
  #Egh lop (Z |T)
  t1<- sum(F_Sym2Vec(n*0.5* Pg * log (phi+(phi==0))))
  t2<-(- 0.5)* sum( ((Pg+diag(q))*omega)*(t(M)%*%M + diag(colSums(S))) ) 
  t3<- n*0.5* sum(log(diag(omega)))  - q*n*0.5*log(2*pi)
  T1<-t1+t2+t3

  # Eglog(p) - Eg log(g)
  T2<-sum(F_Sym2Vec(Pg) * (log(F_Sym2Vec(W)+(F_Sym2Vec(W)==0)) - log(F_Sym2Vec(Wg)+(F_Sym2Vec(Wg)==0)) )) - logSTW+ logSTWg
  
  #Eh log h(Z), reste constant car omegaH fixé à 1 pour identifiabilité 
  T3<- (0.5*sum(log(S))) + n*q*0.5*(1+log(2*pi))
  
  J=T1+T2+T3
  if(is.nan(J)) browser()
  return(c(J=J, T1=T1, T2=T2, t1=t1, t2=t2))
}


########################################
# GENERAL
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
EdgeProba <- function(W, verbatim=FALSE, p.min=1e-16){
  logWcum = logSumTree(W)$det
  if(!isSymmetric(W))  cat('Pb: W non symmetric!')
  p = nrow(W); P = matrix(0, p, p)
  #core of computation
  sapply(1:(p-1),
         function(j){
           sapply((j+1):p,
                  function(k){
                    W_jk = W; W_jk[j, k] = W_jk[k, j] = 0 #kills kj edge in W_kj
                    P[k, j] <<- 1 - exp(logSumTree(W_jk)$det - logWcum )
                    P[j, k] <<- P[k, j]
                  })
         })
  # if(sum(is.na(P))!= 0) browser()
  if(length(which(P<p.min))!=0){
    
    P[which(P<p.min)]= 0
  }
  
  diag(P)=0
  return(P)
}
Kirshner <- function(W){
  # W = squared weight matrix
  # Kirshner (07) formulas
  # W = beta.unif*phi
  p = nrow(W)
  if(sum(colSums(W)==0 )==0){
    sums=apply(W,2,function(x){ sum(x!=0)})
    suspects_paires=which(sums==1)
    if(length(suspects_paires)<2){
      index=1
    }else{
      mate=sapply(suspects_paires, function(suspect){
        which(W[,suspect]!=0)
      })
      if(length(intersect(mate, suspects_paires))!=0){
        index=intersect(mate, suspects_paires)[1]
      }else{
        index=1
      }
    }
  }else{
    index=which(colSums(W)==0)
  }
  L = Laplacian(W)[-index,-index]
 # if(!is.finite(sum(L))) browser()
  # Leigen = eigen(L)
  # Q = (Leigen$vectors) %*% diag(1/Leigen$values) %*% t(Leigen$vectors)
 #  Q = inverse.fractional(L)
  Q = inverse.gmp(L)
  Q = rbind(c(0, diag(Q)),
            cbind(diag(Q), (diag(Q)%o%rep(1, p-1) + rep(1, p-1)%o%diag(Q) - 2*Q)))
  Q = .5*(Q + t(Q))
  P = W * Q
  P = .5*(P + t(P))
  return(P)
}

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
  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Y = matrix(rpois(n*p, exp(Z+prod )), n, p)
  return(list(Y=Y, Z=Z))
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
  sigma = cov2cor(solve(omega))
  sim=list(sigma=sigma,omega=omega,cste=lambda)
  return(sim)
}
#rewrite data from scratch
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
data_from_scratch<-function(type, p=20,n=50, r=5, covariates=NULL, prob=log(p)/p,
                            dens=log(p)/p, signed=FALSE,v=0,draw=FALSE){
  graph<- generator_graph(graph=type,p=p,prob=prob,dens=dens,r=r)
  param<-generator_param(G=as.matrix(graph),signed=signed,v=v)
  data<-generator_PLN(param$sigma,covariates,n)
  if(draw){
    g=as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="nicely")+
      geom_edge_link()+
      geom_node_point(size=2, color="blue")
    print(g)
  }
  return(list(data=data,omega= param$omega))
}
mytheme.light <- list(theme_light(), scale_color_brewer("",palette="Set3"),guides(color=FALSE),
                      theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                            plot.title = element_text(hjust = 0.5)))

mytheme.dark <-function(legend){list= list(theme_light(), scale_color_brewer(legend,palette="Dark2"), scale_fill_brewer(legend,palette="Dark2"),
                                           theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                                                 plot.title = element_text(hjust = 0.5)))
return(list)}
# mypal<-c(brewer.pal(3, "Blues"),brewer.pal(3, "Reds"),brewer.pal(3, "Greens"))
# mypal<-c(brewer.pal(8, "Dark2"),"blue","red")
mytheme <- list(theme_light(), scale_fill_brewer("",palette="Dark2"),#scale_colour_hp_d(option = "RavenClaw", name = ""), #scale_color_manual("",values=mypal),
                theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                      plot.title = element_text(hjust = 0.5)))

################
# diagnostic functions
auc<-function(pred,label){ #require(ROCR)
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}
accppvtpr<-function(probs,omega,h, seuil=0.5){
  Acc=round(mean(1*(probs>seuil)==omega),2) #(TP+TN)/(TP+FP+TN+FN)
  AccH=round(mean(1*(probs[h,]>seuil)==omega[h,]),2)
  AccO=round(mean(1*(probs[-h,-h]>seuil)==omega[-h,-h]),2)
  PPV=round(sum((omega!=0)*(probs>seuil))/(sum((omega!=0)*(probs>seuil))+ sum((omega==0)*(probs>seuil))),2)#TP/(TP+FP)
  PPVH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/(sum((omega[h,]!=0)*(probs[h,]>seuil))+ sum((omega[h,]==0)*(probs[h,]>seuil))),2)
  PPVO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))+  sum((omega[-h,-h]==0)*(probs[-h,-h]>seuil))),2)
  
  TPR=round(sum((omega!=0)*(probs>seuil))/sum(omega!=0), 2)
  TPRH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/sum(omega[h,]!=0), 2)
  TPRO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/sum(omega[-h,-h]!=0), 2)
  
  return(c(Acc,AccH,AccO,PPV,PPVH,PPVO,TPR,TPRH,TPRO))
}
courbes_seuil<-function(probs,omega,h,seq_seuil){
  tmp=sapply(seq_seuil,function(x)  accppvtpr(seuil=x,probs=probs,omega=omega,h=h))
  res=data.frame(cbind(t(tmp),seq_seuil))
  colnames(res)=c("Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","seuil")
  return(as_tibble(res))
}
compute_nSNR<-function(K, indexmissing){
  H=indexmissing ; p=ncol(K)
  O=(1:p)[-H]
  num=norm(K[O,H]%*%solve(K[O,O])%*%K[H,O],type='F')^2
  denom=norm(K[H,H], type='F')^2
  return(num/denom)
}
################
# graphics
ggimage<-function(data){
  melted_data <- melt(data)
  ggplot(melted_data, aes(x=Var1, y=Var2, fill=value)) + theme_light()+labs(x="",y="")+
    geom_tile() +guides(fill=FALSE)+ theme(plot.title = element_text(size=10, hjust=0.5))+ coord_fixed()
}
plotVEM<-function(probs,omega,r,seuil){
  # plots heatmaps for the chosen threshold and print verdicts as title
  q=ncol(omega)
  h=(q-r):q
  performance=accppvtpr(probs,omega,h,seuil)
  Acc=performance[1] ;AccH=performance[2] ;AccO=performance[3] 
  PPV=performance[4] ;PPVH=performance[5] ; PPVO=performance[6]
  TPR=performance[7] ;TPRH=performance[8] ;TPRO=performance[9] 
  p1<-ggimage(probs)+labs(title=paste0("G hat"))
  p2<-ggimage(omega)+labs(title="True G")
  auc<-round(auc(pred = probs, label = omega),3)
  grid.arrange(p1,p2,ncol=2, top=paste0("Tpr=",TPR," (TprO=",TPRO," , TprH=",TPRH,
                                        ")\n Ppv=",PPV," (PpvO=",PPVO," , PpvH=",PPVH,")",
                                        "\n AUC=",auc))
}


plotVerdict<-function(values,colonne){
  # plots verdict curves along threshold from VEM result
  colonne<-enquo(colonne)
  values %>%as_tibble() %>% gather(key,value, - !!colonne)%>%
    mutate(status=ifelse(substr(key,4,4)=="","all",substr(key,4,4))) %>% 
    mutate(key=substr(key,1,3)) %>% spread(key,value)%>% 
    mutate(PPV=ifelse(PPV<0,NA,PPV)) %>% 
    ggplot(aes(TPR,PPV,color=status))+
    geom_rect(aes(xmin=0.5, xmax=1, ymin=0.5, ymax=1), fill="gray90", color="gray90",alpha=0.1)+
    geom_point()+  geom_line()+
    facet_wrap(~status)+mytheme.dark("")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#===========
VE<-function(MO,SO,SH,omega,W,Wg,MH,Pg,logSTW,logSTWg,eps, alpha, filterWg=FALSE,plot=FALSE,verbatim, hist=FALSE){
  #--Setting up
  t1=Sys.time()
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); trim=FALSE
  hidden=(ncol(omega)!=ncol(MO))
  if(hidden){
    H=(p+1):ncol(omega);r=length(H)
    SH <-matrix(1,n,r) # SH fixed to 1 for identifiability
    S<-cbind(SO,SH)# all SHi have same solution, depending only on Mstep
    M=cbind(MO,MH) 
  }else{
    S=SO
    r=0
    M=MO
  }
  LB0=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p, logSTW=logSTW,logSTWg=logSTWg)[1]
  phi=CorOmegaMatrix(omega)
  
  #-- Updates
  #--- MH 
  if(hidden){
    MH.new<- (-MO) %*% (Pg[O,H] * omega[O,H])/diag(omega)[H]
    #if(sum(MH.new!=0)==0) browser()
    MH=MH.new
    M=cbind(MO,MH)
    LB1=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p, logSTW=logSTW,logSTWg=logSTWg),"MH")
  }else{LB1=LB0}
  
  #--- Wg
  compWg= computeWg(phi,omega,W,MH,MO,S=S,alpha,hidden=hidden, hist=hist, verbatim=verbatim)  
  
  Wg.new = compWg$Wg
  #if(sum(colSums(Wg.new)==0)!=0) browser()
 # Pg.new = suppressWarnings( EdgeProba(Wg.new))
  
  logSTWg.tot=logSumTree(Wg.new)
  logSTWg.new=logSTWg.tot$det
  max.prec=logSTWg.tot$max.prec
  Pg.new=Kirshner(Wg.new)

  # if(sum(is.na(Pg.new))!= 0){# si pb pas de pb
  #   message("NaNs in Pg")
  #   Wg.new=Wg
  #   Pg.new=Pg
  # }
  
  LB2=LowerBound(Pg = Pg.new, omega=omega, M=M, S=S,W=W, Wg=Wg.new,p, logSTW=logSTW,logSTWg=logSTWg.new)
  if(filterWg){ 
    bool=as.numeric(LB2[1])>as.numeric(LB1[1])
    if(bool ){
      Wg=Wg.new
      Pg=Pg.new
      logSTWg=logSTWg.new
    }else{message("no Wg update")
      LB2=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p, logSTW=logSTW,logSTWg=logSTWg)}
  }else{
    Wg=Wg.new
    Pg=Pg.new
    logSTWg=logSTWg.new}
  
  LB2=c(LB2,"Wg")
  
  #-- end
  if(hidden){ LB= rbind(LB1, LB2)
  }else{ LB=LB2 }
  
  res=list(Gprobs=Pg,Gweights=Wg,Hmeans=M,Hvar=S,LB=LB,logSTWg=logSTWg,max.prec=max.prec )
  return(res)
}

#===========
Mstep<-function(M,S,Pg, omega,W,logSTW,logSTWg, plot=FALSE,eps, verbatim=FALSE,Wg,p, filterDiag,iterVEM){
  n=nrow(S)  ; O=1:p ; q=ncol(omega) ; iterM=0 ; diff=1
  hidden=(q!=p)
  if(hidden) H=(p+1):q
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  
  #--- Omega
  omega2=omega
  for(i in 1:2){
    omDiag1 <- diag(omega) # omegaH bouge aussi
    omega1=computeOffDiag(omDiag1, SigmaTilde,p)
    #if(sum(omega1[O,H]!=0)==0) browser()
    omDiag2 <- omegaDiag(Pg,omega2,SigmaTilde) # omegaH bouge aussi
    omega2=computeOffDiag(omDiag2, SigmaTilde,p)
    #if(sum(omega2[O,H]!=0)==0) browser()
  }
 
  LB11=LowerBound(Pg = Pg, omega=omega1, M=M, S=S,W=W, Wg=Wg,p,logSTW=logSTW,logSTWg=logSTWg)
  LB12=LowerBound(Pg = Pg, omega=omega2, M=M, S=S,W=W, Wg=Wg,p,logSTW=logSTW,logSTWg=logSTWg)
  
  if(iterVEM>3 || !hidden  ){
    if(filterDiag){
       if(as.numeric(LB11[1])>as.numeric(LB12[1])){
      message("no diag")
      omega=omega1
      LB1=c(LB11,"omega")
    }else{
      omega=omega2
      LB1=c(LB12,"omega")}
    }else{
      omega=omega2
      LB1=c(LB12,"omega")}
  }else{
    omega=omega2
    LB1=c(LB12,"omega")
  }
  

  #--- Beta
  #  while((diff>eps) && (iterM <= 2)){
  #  iterM=iterM+1
  for(i in 1:2){
    Mei=Meila(W)
    if(sum(!is.finite(F_Sym2Vec(Meila(W))))!=0 ||  sum(F_Sym2Vec(Meila(W))<=0)!=0){
      Wtest=W*100
      if( sum(F_Sym2Vec(Meila(Wtest))==0)==0){
        Mei=Meila(Wtest)
      }else{
        Mei = Mei+1e-14
      }
    }
    # c=0.1
    # while(min(eigen(Laplacian(W)[-1,-1])$values)<0){
    #   W=W+diag(c,ncol(W))
    #   c=c*1.1
    # }
    W.new= Pg/(Mei )
   # if(sum(!is.finite(F_Sym2Vec(W.new)))) browser()
    if(length(which(W.new< 1e-16))!=0){
      W.new[W.new< 1e-16] = 0 # pas de bmax
    }
    if(hidden) W.new[H,H]=0
    diag(W.new)=0
    #if(is.na(logSumTree(W.new)$det )){
    norm1=10^(round(log10(min(W[W!=0]))-(-14))) #division maximale pour ne pas tuer les poids les plus petits
    norm2=max(10^(round(log10(max(W))-2)),1)
    
    W.new=W.new/(min(norm1, norm2))
    #}
    # diff=max(abs(F_Sym2Vec(W.new)-F_Sym2Vec(W)))
    W=W.new
    # }
  }
  logSTW.tot=logSumTree(W)
  logSTW=logSTW.tot$det
  max.prec=logSTW.tot$max.prec
  LB2=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p,logSTW=logSTW,logSTWg=logSTWg),"W")
  res=list(W=W, omega=omega,  LB=rbind(LB1,LB2), logSTW=logSTW,max.prec=max.prec) 
  return(res)
}

Mstep_nobeta<-function(M,S,Pg, omega,W,Wg,p){
  n=nrow(S) 
  O=1:p ; q=ncol(omega)
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  for(i in 1:2){
    omDiag <- omegaDiag(Pg,omega,SigmaTilde) # omegaH bouge aussi
    omega.new=computeOffDiag(omDiag, SigmaTilde,p)
    omega=omega.new
  }
  
  LB1=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p),"omega")
  res=list(W=W, omega=omega,  LB=LB1) 
  return(res)
}

#===========
VEMtree<-function(counts,MO,SO,MH,ome_init,W_init,Wg_init, maxIter=20,eps=1e-2, alpha,
                  verbatim=TRUE, plot=FALSE, print.hist=FALSE, filterWg=FALSE,filterDiag=FALSE,nobeta=FALSE){
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO)
  hidden=!is.null(MH)
  if(hidden){  
    H=(p+1):ncol(ome_init);r=length(H)
    SH=matrix(1,nrow=n,ncol=r)
    #   SH <-matrix(1/(diag(omega)[H]),n,r, byrow = TRUE)
    M=cbind(MO,MH) ; S=cbind(SO, SH)
  }else{
    r=0
    SH= NULL
    M=MO ; S=SO
  }
  omega=ome_init;  W=W_init;  Wg=Wg_init; Pg=matrix(0.5, ncol(W),ncol(W)) 
  iter=0 ; lowbound=list()
  diffW=c();diffOm=c();diffWg=c();diffPg=c();diffJ=c();diffOmDiag=c(0);diffMH=c(); diffquantile=c() ;diffPhi<-c(); J<-c()
  diffWiter=1 ; diffJ=1 ; J=c()
  max.prec=FALSE
  t1=Sys.time()
  logSTW=logSumTree(W)$det
  logSTWg=logSumTree(Wg)$det
  #(diffW[iter] > eps) ||
  while(  (diffOm[iter] > eps) && (diffJ>eps)  && (iter < maxIter) || iter<2 ){ 
    iter=iter+1 
    # if(iter==28) browser()
    if(verbatim) cat(paste0("\n Iter n°", iter))
    #--- VE
    resVE<-VE(MO=MO,SO=SO,SH=SH,omega=omega,W=W,Wg=Wg,MH=MH,Pg=Pg,logSTW,logSTWg,eps=1e-3,verbatim=verbatim,
              alpha=alpha,filterWg=filterWg, plot=FALSE, hist=print.hist)
    M=resVE$Hmeans 
    S=resVE$Hvar 
    Pg.new=resVE$Gprobs
    Wg.new=resVE$Gweights 
    if(hidden){ 
      SH<-matrix(S[,H],n,r)
      MH.new<-matrix(M[,H],n,r)
      diffMH[iter]<-abs(max(MH.new-MH))
      MH=MH.new
    }
    if(resVE$max.prec) max.prec=TRUE
    diffWg[iter]<-abs(max(Wg.new-Wg))
    diffPg[iter]<-abs(max(Pg.new-Pg))
    Wg=Wg.new
    Pg=Pg.new
    logSTWg=resVE$logSTWg
    #--- M
    if(nobeta){
      resM<-Mstep_nobeta(M=M,S=S,Pg=Pg, omega=omega,W=W,Wg=Wg, p=p)
      omega.new=resM$omega
      diffOm[iter]=abs(max(F_Sym2Vec(omega.new)-F_Sym2Vec(omega)))
      if(iter>1)diffOmDiag[iter]=sum(diag(omega.new)-diag(omega))
      omega=omega.new
      diffW[iter]= diffOm[iter]
      Jiter=as.numeric(tail(resM$LB[1],1))
    }else{
      resM<-Mstep(M=M,S=S,Pg=Pg, omega=omega,W=W,logSTW,logSTWg,plot=FALSE, eps=1e-3 ,Wg=Wg, p=p,filterDiag=filterDiag,
                  iterVEM=iter)
      W.new=resM$W
      #diffW[iter]=abs(max(log(W.new+(W.new==0))-log(W+(W==0))))
      diffW[iter]=abs(max(W.new-W))
      diffWiter=diffW[iter]
      W=W.new
      omega.new=resM$omega
      logSTW=resM$logSTW
      if(resM$max.prec) max.prec=TRUE
      diffOm[iter]=abs(max((omega.new)-(omega)))
      if(iter>1)diffOmDiag[iter]=sum(diag(omega.new)-diag(omega))
      omega=omega.new
      Jiter=as.numeric(tail(resM$LB[,1],1))
    }
    lowbound[[iter]] = rbind( resVE$LB, resM$LB)
    
    diffJ=Jiter - tail(J,1)
    J<-c(J, Jiter)
  }
  ########################
  resVE<-VE(MO=MO,SO=SO,SH=SH,omega=omega,W=W,Wg=Wg,logSTW,logSTWg,MH=MH,Pg=Pg,eps=1e-3,verbatim=verbatim,
            alpha=alpha,filterWg=filterWg, plot=FALSE, hist=print.hist)
  M=resVE$Hmeans 
  S=resVE$Hvar 
  Pg=resVE$Gprobs
  Wg=resVE$Gweights 
  logSTWg=resVE$logSTWg
  if(resVE$max.prec) max.prec=TRUE
  if(hidden){ 
    SH<-matrix(S[,H],n,r)
    MH<-matrix(M[,H],n,r)
  }
  #--- M
  if(nobeta){
    resM<-Mstep_nobeta(M=M,S=S,Pg=Pg, omega=omega,W=W,Wg=Wg, p=p)
    omega=resM$omega
  }else{
    resM<-Mstep(M=M,S=S,Pg=Pg, omega=omega,W=W,logSTW,logSTWg,plot=FALSE, eps=1e-3 ,Wg=Wg, p=p,filterDiag=filterDiag, iterVEM=iter)
    if(resM$max.prec) max.prec=TRUE
    W=resM$W
    omega=resM$omega
  }
  ##########################
 
  lowbound[[iter+1]] = rbind( resVE$LB, resM$LB) 
  lowbound=data.frame(do.call(rbind,lowbound))
  lowbound[,-ncol(lowbound)]<-apply(lowbound[,-ncol(lowbound)],2,function(x) as.numeric(as.character(x)))
  colnames(lowbound)[ncol(lowbound)] = "parameter"
  if(hidden){
    if(nobeta){ features<-data.frame(diffPg=diffPg, diffOm=diffOm, diffWg=diffWg)
    }else{ features<-data.frame(diffPg=diffPg, diffW=diffW, diffOm=diffOm, diffWg=diffWg)
    }
  }else{
    if(nobeta){ features<-data.frame(diffPg=diffPg,  diffOm=diffOm, diffWg=diffWg)
    }else{ features<-data.frame(diffPg=diffPg, diffW=diffW, diffOm=diffOm, diffWg=diffWg)
    }
  }
  t2=Sys.time()
  time=t2-t1
  if(verbatim) cat(paste0("\nVEMtree ran in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal weights difference: ",round(diffW[iter],7)))
  if(plot){
    g1<-features  %>%  rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line() + facet_wrap(~key, scales="free")+
      labs(x="",y="", title="Parameters")+ mytheme.dark("")+guides(color=FALSE)
    g2<- lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>% 
      ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(parameter)), size=3)+geom_line()+
      facet_wrap(~key, scales="free")+
      labs(x="iteration",y="", title="Lower bound and components")+mytheme+
      scale_color_discrete("")
    grid.arrange(g1,g2, ncol=1)
  }
  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,omega=omega, lowbound=lowbound, features=features, finalIter=iter, time=time,max.prec=max.prec))
}

True_lowBound<-function(Y, M,S,theta,X, W, Wg, Pg, omega){
  p=ncol(Y)
  logSTW=logSumTree(W)$det
  logSTWg=logSumTree(Wg)$det
  partJ<-LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p, logSTW, logSTWg)[1]
  partY<-sum(
    -exp(X%*%t(theta) + M[,1:p]+S[,1:p]/2)+ Y*(X%*%t(theta)+M[,1:p])-
      lgamma(Y+1)
  )
  TrueJ<-as.numeric(partJ)+partY
  return(TrueJ)
}

VBIC<-function(TrueJ,p,r,d,n){
  q=p+r
  nbparam<-p*(d) + (p*(p+1)/2 +r*p)+(q*(q-1)/2 - 1) #d comprends l'intercept
  vBIC=TrueJ - nbparam*log(n)/2
  return(vBIC)
}

ICL_T<-function(TrueJ, Wg,Pg,n,r){
  pen_T=-( sum( Pg * log(Wg) ) - logSumTree(Wg)$det )
  # if(pen_T<0) browser()
  ICL=TrueJ - pen_T 
  return(ICL)
}
ICL_ZH<-function(TrueJ,S,n,r){
  q=ncol(S) ; p=q-r
  H=(p+1):q
  if(r!=0){
    pen_ZH= 0.5*sum(log(S[,H])) + n*r*0.5*(1+log(2*pi))
  }else{
    pen_ZH= 0
  }
  ICL=TrueJ - pen_ZH
  return(ICL)
}

ICL<-function(TrueJ, Pg,Wg ,S,n,r,d, omega){
  q=ncol(S) ; p=q-r
  H=(p+1):q
  if(r!=0){
    pen_ZH= 0.5*sum(log(S[,H])) + n*r*0.5*(1+log(2*pi))
  }else{
    pen_ZH= 0
  }
  Pg=Kirshner(Wg)
  pen_T=-( sum( Pg * log(Wg+(Wg==0)) ) - logSumTree(Wg)$det) 
  pen_r<-p*(d) + (p*(p+1)/2 +r*p+r)+(q*(q-1)/2 - 1) #d comprends l'intercept
 # norm=n*q*(q-1)/2
 
  
  ICL=TrueJ   - (pen_T + pen_ZH + pen_r*log(n)/2)
  return(ICL)
}

criteria<-function(List.vem,counts,theta, matcovar,r){
  n=nrow(counts);p= ncol(counts)
  data<-lapply(List.vem,function(vem){
    J<-True_lowBound(counts,vem$M,vem$S, theta, matcovar,vem$W, vem$Wg, 
                     vem$Pg, vem$omega )
    vBIC<-VBIC(J,p,r=r, d=ncol(matcovar), n=n)
    # ICLT<-ICL_T(J, vem$Wg,vem$Pg,n,r)
    # ICLZH<-ICL_ZH(J,vem$S, n,r)
    ICL<-ICL(J, vem$Pg,vem$Wg,vem$S, n,r,d=ncol(matcovar), omega=vem$omega)
    res=data.frame(vBIC=vBIC, ICL=ICL,J)
    return(res)
  })
  res=do.call(rbind,data)
  res$r=r
  return(res)
}

List.VEM<-function(cliquesObj, counts, sigma_obs, MO,SO,r,alpha, cores,maxIter,eps, nobeta){
  p=ncol(counts) ; O=1:p ; n=nrow(counts)
  
  #--- run all initialisations with parallel computation
  list<-mclapply(seq_along(cliquesObj$cliqueList), function(num){
    #init
    c=cliquesObj$cliqueList[[num]]
    init=initVEM(counts = counts, initviasigma=c, sigma_obs,r = r)
    Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
    #run VEMtree
    VEM<-VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=eps, alpha=alpha,maxIter=maxIter, 
                 verbatim = TRUE, print.hist=FALSE, filterWg = TRUE, 
                 filterDiag = TRUE, nobeta=nobeta)
    VEM$clique=c
    VEM$nbocc=cliquesObj$nb_occ[num]
    return(VEM)
  }, mc.cores=cores)
  #--- keep only converged vems
  # converged<-do.call(rbind,lapply(list,function(vem){
  #   diffW=vem$features$diffW
  #   conv=(diffW[length(diffW)]<=eps)}))
  # if(sum(converged)!=0){
  #   list=list[converged]
  # } 
  return(list)
}


logSumTree<-function(W){
  if(sum(colSums(W)==0 )==0){
    sums=apply(W,2,function(x){ sum(x!=0)})
    suspects_paires=which(sums==1)
    if(length(suspects_paires)<2){
      index=1
    }else{
      mate=sapply(suspects_paires, function(suspect){
        which(W[,suspect]!=0)
      })
      if(length(intersect(mate, suspects_paires))!=0){
        index=intersect(mate, suspects_paires)[1]
      }else{
        index=1
      }
    }
  }else{
    index=which(colSums(W)==0)
  }
  mat=Laplacian(W)[-index, -index]
  max.prec=FALSE
  #output=suppressWarnings(log(det(mat)))
  output=det.fractional(mat, log=TRUE)
  if(output==log(.Machine$double.xmax )){
     max.prec=TRUE
     message("max.prec!")
  }
  # output=tryCatch({det.fractional(mat, log=TRUE)},# exact computation
  #                 error=function(e){ # if error, trim matrix 
  #                   trim=TRUE
  #                   W=W/10
  #                   W[W<1e-16]=0
  #                   output=log(det(Laplacian(W)[-1,-1]))
  #                   return(output)}, finally={})
  # if(!is.finite(output)){
  #   if(det(mat)<0) message("exact det na")
  #   if(!is.finite(det(mat))) message("exact det inf")
  #   output=  tryCatch({det.fractional(mat, log=TRUE)},# exact computation
  #                     error=function(e){ # if error, trim matrix 
  #                       trim=TRUE
  #                       W=W/10
  #                       W[W<1e-16]=0
  #                       output=log(det(Laplacian(W)[-1,-1]))
  #                       return(output)}, finally={})
  # }
  return(list(det=output,max.prec=max.prec))
}
