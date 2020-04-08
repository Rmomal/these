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
  
  r=length(cliqueList) ; p=ncol(Sigma)
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
  sigmaFull <- c(sigma, rep(1, r)) 
  SigmaFull <- diag(sigmaFull) %*% CorrFull %*% diag(sigmaFull)
 
  # Initialising Omega
  OmegaFull <- tryCatch({solve(SigmaFull)},
                        error=function(e){browser()},finally={})
  coefDiag <- c(rep(1, p), 1/sqrt(diag(OmegaFull)[p+(1:r)]))
  OmegaFull <- diag(coefDiag) %*% OmegaFull %*% diag(coefDiag)
  OmegaFull[(p+1):(p+r),(p+1):(p+r)]<-diag(1,r)
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
  Wg_init <- matrix(1, p+r, p+r); Wg_init =Wg_init / sum(Wg_init)
  W_init <- matrix(1, p+r, p+r); W_init =W_init / sum(W_init)
  W_init[1:p,1:p] <- EMtree_corZ(cov2cor(sigma_obs),n = n,maxIter = 20,verbatim = FALSE)$edges_weight
  diag(Wg_init) = 1;diag(W_init) = 1
  # Z
  if(r!=0){
    initial.param<-initEM(sigma_obs,n=n,cliqueList = (initviasigma),cst=1.05, pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
    omega_init=initial.param$K0
    MHinit<-sapply(initviasigma, function(clique){
      pr=prcomp(t(counts[,clique]),scale. = FALSE)
      res = matrix(pr$rotation[,1]*pr$sdev[1],nrow=n,ncol=1)
      return(res)
    })
  }else{#init with no missing actors
    omega_init=solve(sigma_obs)
    MHinit=NULL
  }
  return(list(Wginit= Wg_init, Winit= W_init, omegainit=omega_init,MHinit=MHinit))
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

argminKL <- function(gamma, Pg, M,S,omega,phi,W,p ){
  r=ncol(omega)-p ; O = 1:p ; n=nrow(M)
  hidden= (r!=0)
  if(hidden){
    H=(p+1):(p+r)
  }
  EhZoZo = t(M[,O])%*%M[,O]+ diag(colSums(S[,O]))
  if(hidden){
    KL <- 0.5 * sum(diag( (t(M[, H]) %*% M[, H] + sum(S[, H])))) + # omegaH is 1
      sum(diag((Pg[O, H] * omega[O, H]) %*% (t(M[, H]) %*% M[, O]))) +
      0.5 * sum(diag((Pg[O, O] * omega[O, O]) %*% EhZoZo)) -
      n * 0.5 * sum(log(diag(omega))) -
      n * 0.5 * 2 * sum(F_Sym2Vec(Pg) * log(F_Sym2Vec(phi))) +
      2 * sum(F_Sym2Vec(Pg) * (gamma - log(F_Sym2Vec(W)))) -
      log(SumTree(F_Vec2Sym(exp(gamma)))) +
      log(SumTree(W)) - 
      0.5 * sum(log(S[, H]))#-lambda*(sum(exp(gamma)) - 0.5)
  }else{
    KL <-    0.5 * sum(diag((Pg[O, O] * omega[O, O]) %*% EhZoZo)) -
      n * 0.5 * sum(log(diag(omega))) -
      n * 0.5 * 2 * sum(F_Sym2Vec(Pg) * log(F_Sym2Vec(phi))) +
      2 * sum(F_Sym2Vec(Pg) * (gamma - log(F_Sym2Vec(W)))) -
      log(SumTree(F_Vec2Sym(exp(gamma)))) +
      log(SumTree(W)) 
  }
  return(KL)
}
findAlpha<-function(log_expr,cond.tol=1e-10){
  alpha.grid = (1:n)/n
  #----
  condtest<-tryCatch({
    expr=sapply(alpha.grid,function(alpha){
      test=F_Sym2Vec(log_expr*alpha)
      test=test-mean(test)
      test=F_Vec2Sym(exp(test))
      lambda = svd(test)$d
      cond = min(abs(lambda))/max(abs(lambda))
    })}, error=function(e){browser()}, finally={})
  
  if(condtest[length(condtest)]>cond.tol){
    alpha=1
  }else{
    alpha=alpha.grid[max(which(condtest>cond.tol))]
  }
  return(alpha)
}
computeWg<-function(phi,omega,W,MH,MO,S, alpha,hidden=TRUE, hist=FALSE, verbatim=FALSE ){
  
  p=ncol(MO); O = 1:p ; n=nrow(MO)
  if(hidden){  r=ncol(MH)
  H = (p+1):(p+r)  
  M=cbind(MO,MH)
  }else{
    r=0
  }
  
  logWg<-matrix(0,(p+r),(p+r)) ; Wg<-matrix(0,(p+r),(p+r))
  psi=phi^(n*0.5)
  
  logWg[O,O]<-log(W[O,O])+alpha*(n*0.5*log(phi[O,O])-0.5*omega[O,O]*(t(MO)%*%MO))
  diag(logWg) = 0
  #---
  
  # alpha.grid=(1:n)/n
  # condsans1<-sapply(alpha.grid ,function(alpha){
  #   test3=F_Sym2Vec(alpha*(n*0.5*log(phi)-0.5*omega*(t(M)%*%M)))
  #   test3=F_Vec2Sym(exp(test3))
  #    test3=test3-mean(test3)
  #   lambda = svd(test3)$d
  #   cond3 = min(abs(lambda))/max(abs(lambda))
  # })
  # datacond<-data.frame(alpha=alpha.grid,   SansDiag1=log(condsans1),   Diag1=log(condDiag1),
  #                      SansDiag2=log(condsans2),   Diag2=log(condDiag2))
  # g=datacond %>% gather(key, value, -alpha) %>% ggplot(aes(alpha,value, color=key))+geom_point()+mytheme+
  #   geom_hline(yintercept = -23)
  # print(g)
  binf=-20
  bsup=20
  # shrinking and centering
  gammaO=F_Sym2Vec(logWg[O,O])
  if(verbatim)  cat(paste0(", rangeO=", round(max(gammaO)-min(gammaO),5)))
  if(hist){ par(mfrow=c(2,2))
    hist(gammaO, breaks=20, main="O in")
  }
  gammaO=gammaO-mean(gammaO)
  gammaO[which(gammaO<(binf))]=binf
  gammaO[which(gammaO>bsup)]=bsup
  if(hist) hist(gammaO, breaks=20, main="O out")
  
  Wg[O,O]=exp(F_Vec2Sym(gammaO))
  diag(Wg) = 1
  
  if(hidden){
    
    logWg[O,H]<-log(W[O,H])+alpha*(n*0.5*log(phi[O,H])-omega[O,H]*t((t(MH)%*%MO)))
    logWg[H,O]<-t(logWg[O,H])
    
    gammaOH=logWg[O,H]
    if(verbatim)  cat(paste0(", rangeH=", round(max(gammaOH)-min(gammaOH),5)))
    if(hist){
      hist(gammaOH, breaks=20, main="OH in")
    } 
    gammaOH=gammaOH-mean(gammaOH)
    gammaOH[which(gammaOH<(binf))]=binf
    gammaOH[which(gammaOH>bsup)]=bsup
    if(hist) hist(gammaOH, breaks=20, main="OH out")
    Wg[O,H]=exp(gammaOH)
    Wg[H,O]<-t(Wg[O,H])
  }  
  return(Wg)
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
  return(CorOmega)
}

#=====
# M step

argmaxJ<-function(gamma,Pg,omega,sigmaTilde,phi,n, trim=TRUE, verbatim=TRUE){
  test= log(SumTree(exp(gamma)))
  if(is.nan(test)) gamma<-log(exp(gamma)+1)
  maxJ <- sum(Pg*F_Vec2Sym(F_Sym2Vec(gamma) + n*0.5*F_Sym2Vec(log(phi)))) + 
    n*0.5*sum(log(diag(omega))) - log(SumTree(exp(gamma))) - 
    0.5*n*sum(diag( Pg * omega %*% sigmaTilde))# - n*0.5*sum(diag(omega*sigmaTilde)) #termes diagonaux
  # + lambda*(sum(W)-0.5)
  return(maxJ)
}

# Diangonal no longer updated 
optimDiag<-function(Pg,omega,SigmaTilde){ #closed form
  quantity<-Pg*omega*SigmaTilde
  diag(quantity)=0
  vectorSuml<-colSums(quantity)
  vectorSuml[vectorSuml>1]<-1-1e-6
  omegaDiag = (1-vectorSuml)/diag(SigmaTilde)
  return(omegaDiag)
}
optimDiag2 <- function(omega_ii,i,  omega, SigmaTilde, Pg) { #for numerical optimization
  q <- nrow(omega)
  phi = CorOmegaMatrix(omega)
  newphi=(1- phi)/phi
  diag(newphi)=0
  quotient = diag(Pg %*% newphi)[i]
  
  grad=(1/omega_ii)*(quotient+1) - SigmaTilde[i, i]
  return( grad)
}
dichotomie <- function(a, b, F.x, eps, bool=FALSE){
  
  x.min = ifelse(F.x(0) >0,b,1e-4);
  while(F.x(x.min)<0){x.min = x.min -x.min/2}
  x.max = a
  while(F.x(x.max)>0){x.max = x.max * 2}
  #  cat(paste0("\n [xmin ; xmax] = [",round(x.min,3),",",round(x.max,3),"]"))
  x = (x.max+x.min)/2
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  iter=0
  while(abs(x.max-x.min) > eps){
    iter=iter+1
    
    if(f < 0) {
      x.max = x
      f.max = f
    }else{
      x.min = x
      f.min = f
    }
    x = (x.max+x.min)/2;
    f = F.x(x)
  }
  #  cat(paste0("f=",round(f,8),", om_ii=",round(x,3),"\n"))
  return(x)
}

computeOffDiag<-function(omegaDiag,SigmaTilde,p){ 
 # browser()
  q=length(omegaDiag)   ;
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
    omega[H,H]<-diag(1,length(H))
  }

  return(omega)
}

computeOmegaOH<-function(omegaDiag,SigmaTilde,p){
  q=length(omegaDiag)
  omegaOH=matrix(0, q-1, q-p)
  sapply(1:(q-1),
         function(k){
           sapply((p+1):q,
                  function(j){
                    omegaOH[k,(j-p)] <<- (1 - ( 1+4*SigmaTilde[k,j]^2*omegaDiag[j]*omegaDiag[k])^0.5)/(2*SigmaTilde[k,j])
                  })
         } )
  
  return(omegaOH)
}


#=====
#Lower bound

LowerBound<-function(Pg ,omega, M, S, W, Wg,p){
  n=nrow(M)
  q=nrow(omega)
  O=1:p ; r=q-p
  psi=CorOmegaMatrix(omega)
  
  #Egh lop (Z |T)
  t1<- sum(F_Sym2Vec(n*0.5* Pg * log (psi )))+ n*0.5* sum(log(diag(omega)))  - q*n*0.5*log(2*pi)
  t2<-(- 0.5)* sum( (Pg*omega)*(t(M)%*%M + diag(colSums(S))) )
  T1<-t1+t2
  # Eg log(p) - Eg log(g)
# browser()
  test= log(SumTree(Wg))
  if(is.nan(test)) Wg=Wg+1
  test= log(SumTree(W))
  if(is.nan(test)) W=W+1
  T2<-sum(F_Sym2Vec(Pg) * (log(F_Sym2Vec(W)) - log(F_Sym2Vec(Wg)) )) - log(SumTree(W))+ log(SumTree(Wg))

  #Eh log h(Z), reste constant car omegaH fixé à 1 pour identifiabilité 
  T3<- 0.5*sum(log(S)) + n*q*0.5*(1+log(2*pi))
  
  J=T1+T2+T3
  if(is.nan(J)) browser()
  return(c(J=J, T1=T1, T2=T2, t1=t1, t2=t2))
}



#####################
# For double hidden proba, currently not useful

HiddenEdgeProba <- function(W,r=1, verbatim=FALSE){
  #computes the probability for two links to a hidden covariates to NOT be there
  #Coded for 1 hidden covariate
  it=-1
  Wcum = SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  q = nrow(W); P = matrix(0, q, q)
  p=q-1 #one h
  h=q # the hidden covariate is stored in last row/column
  #core of computation
  sapply(1:p-1,
         function(j){
           sapply((j+1):p, #visits all combinations of observed nodes
                  function(k){
                    W_jk = W 
                    W_jk[h, k] = W_jk[k, h] = 0 #kills kh edge in W_kj
                    W_jk[h, j] = W_jk[j, h] = 0 #kills jh edge in W_kj
                    P[k, j] <<- 1 - SumTree(W_jk) / Wcum
                    P[j, k] <<- P[k, j]
                  }
           )
         }
  )
  P[which(P<1e-10)]=1e-10  
  diag(P)=0
  return(P)
}
CgMatrix<-function(Pg,Pghkl,omega,p){ # code for 1 hidden covariate. Otherwise, must be a sum on h
  h=ncol(omega)
  omegaH=omega[h,h]
  Cg = matrix(1/omegaH,p,p)
  
  sapply(1:p-1,
         function(j){
           sapply((j+1):p, #visits all combinations of observed nodes
                  function(k){
                    Cg[k, j] <<- omega[h,k]*omega[h,j] * (Pg[h,k]+Pg[h,j]+Pghkl[k,j] -1)
                    Cg[j, k] <<- Cg[k, j]
                  }
           )
         })
  return(Cg)
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
EdgeProba <- function(W, verbatim=FALSE, p.min=1e-10){
  it=-1
  Wcum = SumTree(W)
  test= (SumTree(W))
  if(is.nan(test)) W=W+1
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  p = nrow(W); P = matrix(0, p, p)
  #core of computation
  sapply(1:(p-1),
         function(j){
           sapply((j+1):p,
                  function(k){
                    W_jk = W; W_jk[j, k] = W_jk[k, j] = 0 #kills kj edge in W_kj
                    P[k, j] <<- 1 - SumTree(W_jk) / Wcum
                    P[j, k] <<- P[k, j]
                  }
           )
         }
  )
  P[which(P<p.min)]=p.min
  diag(P)=0
  return(P)
}

generator_PLN<-function(Sigma,covariates=NULL, n=50){
  # ajout d'une constante, pa rapport à EMtree::generator_PLN
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

courbes_seuil<-function(probs,omega,h,seq_seuil){
  tmp=sapply(seq_seuil,function(x)  accppvtpr(seuil=x,probs=probs,omega=omega,h=h))
  res=data.frame(cbind(t(tmp),seq_seuil))
  colnames(res)=c("Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","seuil")
  return(as_tibble(res))
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
    facet_wrap(~status)+mytheme.dark+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

compute_nSNR<-function(K, indexmissing){
  H=indexmissing ; p=ncol(K)
  O=(1:p)[-H]
  num=norm(K[O,H]%*%solve(K[O,O])%*%K[H,O],type='F')^2
  denom=norm(K[H,H], type='F')^2
  return(num/denom)
}

#===========
VE<-function(MO,SO,SH,omega,W,Wg,MH,Pg,maxIter,minIter,eps, alpha,beta.min=exp(-20),
             plot=FALSE,verbatim=FALSE, hist=FALSE, filterPg=TRUE, filterWg=TRUE){
  #--Setting up
  t1=Sys.time()
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); 
  hidden=(ncol(omega)!=ncol(MO))
  if(hidden){
    H=(p+1):ncol(omega);r=length(H)
    # SH <-diag(1,r) #omegaH is set to diag 1
    S<-cbind(SO, matrix(rep(1,n),n,r))# all SHi have same solution, depending only on Mstep
    M=cbind(MO,MH) 
  }else{
    S=SO
    r=0
    M=MO
  }
  
  LB=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)[1]
  phi=CorOmegaMatrix(omega)
  KL0 = argminKL(F_Sym2Vec(log(Wg)), Pg, M,S,omega,phi,W,p)
  KL.new=KL0
  #-- Probabilities estimates
  Pg.new = EdgeProba(Wg)
  #  LB0=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)[1]
  LB0=LowerBound(Pg = Pg.new, omega=omega, M=M, S=S,W=W, Wg=Wg,p)
  if(filterPg){
    if(as.numeric(LB0[1])>as.numeric(LB)){
      Pg=Pg.new
    }
  }else{ Pg=Pg.new}
  LB0=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p),"Pg")
  #browser()
  #-- Updates
  #--- MH 
  if(hidden){
    MH.new<- (-MO) %*% (Pg[O,H] * omega[O,H]) #omegaH=diag(1)
    diffMH<-max(abs(MH-MH.new))
    M.new=cbind(MO,MH.new)
    LB1=LowerBound(Pg = Pg, omega=omega, M=M.new, S=S,W=W, Wg=Wg,p)
    # if(LB1[1]>LB0){
    MH=MH.new
    M=cbind(MO,MH)
    KL.new<-argminKL(F_Sym2Vec(log(Wg)), Pg, M,S,omega,phi,W,p)
    # } 
    LB1=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p),"MH")
  }else{LB1=LB0}
  #--- Wg
  Wg.new= computeWg(phi,omega,W,MH,MO,S=S,alpha,hidden=hidden, hist=hist, verbatim=FALSE)  
  diag(Wg.new)=1
  Wg.new[which(Wg.new< beta.min)] = beta.min
  diffW=max(abs((F_Sym2Vec(Wg.new))-(F_Sym2Vec(Wg))))
  LB2=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg.new,p)
  if(filterWg){ 
    bool=tryCatch({as.numeric(LB2[1])>as.numeric(LB0[1])},
             error=function(e){browser()}, finally={})
    if(bool){
      Wg=Wg.new
      KL.new<-argminKL(F_Sym2Vec(log(Wg.new)), Pg, M,S,omega,phi,W,p)
    }
  }else{
    Wg=Wg.new}
  LB2=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p),"Wg")
  diffKL = KL.new - KL0
  #-- end
  
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nVE step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal W difference: ",round(diffW[iter],5),
                          "\nFinal KL difference: ",round(diff[iter],4)))
  if(plot){
    if(hidden){
      plotdata=data.frame(diff.W=diffW,  diff.KL=diffKL,diff.MH=diffMH, PartofKL=KL.new)
    }else{
      plotdata=data.frame(diff.W=diffW,  diff.KL=diffKL, PartofKL=KL.new)
    }
    g=plotdata %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+ facet_wrap(~key, scales="free")+ labs(x="iter",y="") +
      mytheme.dark
    print(g)
  }
  if(hidden){
    LB= rbind(LB0,LB1, LB2)
    res=list(Gprobs=Pg,Gweights=Wg,Hmeans=M,Hvar=S, KL=KL.new, diff=diffMH, diffW=diffW,LB=LB)
  }else{
    LB= rbind(LB0, LB2)
    res=list(Gprobs=Pg,Gweights=Wg,Hmeans=M,Hvar=S, KL=KL.new, diffW=diffW,LB=LB)
  }
  
  return(res)
}

#===========
Mstep<-function(M,S,Pg, omega,W,maxIter, beta.min,beta.max, trim=TRUE,plot=FALSE,eps, verbatim=FALSE,Wg,p){
  t1=Sys.time()
  diffJ=(1); diff.J=c(1);  diff.W=c(0.05);diffW=1
  maxJ=c()
  n=nrow(S) 
  O=1:p ; q=ncol(omega)
  hidden=(q!=p)
  if(hidden) H=(p+1):q
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  iter=0
  omegaDiag=diag(omega)
  LB0=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)[1]
  
  omega.new=computeOffDiag(omegaDiag, SigmaTilde,p)
  diffinvSigO= max(abs(solve(omega.new)[O,O] - SigmaTilde[O,O]))
  # diffinvSigO= max(abs(solve(SigmaTilde[O,O])-((omega.new[O,O]-1/omega.new[H,H]*matrix(omega.new[O,H], p, length(H))%*%matrix(omega.new[H,O],length(H),p)))))
  if(verbatim) cat(paste0(", diffinvSigO=",round(diffinvSigO,5)))
  LB1=LowerBound(Pg = Pg, omega=omega.new, M=M, S=S,W=W, Wg=Wg,p)
  #if(LB1[1]>LB0){
  omega=omega.new
  #}
  LB1=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p),"omega")
  phi=CorOmegaMatrix(omega) # de omega ou omega.new ?
  Wsave=W
  while((diffW>eps) & (iter <= maxIter)){
    iter=iter+1
    if(verbatim) cat(paste0("\n Mstep:",iter," diffW=",round(diffW,4)))
    
    Mei=Meila(W)  # c'est Meila qui justifie de faire un while
    W.new= Pg/(Mei) 
    
    diag(W.new)=1
    W.new[which(W.new< beta.min)] = beta.min
    W.new[which(W.new> beta.max)] = beta.max
    diffW=max(abs(F_Sym2Vec(W.new)-F_Sym2Vec(W)))
    
    maxJ[iter]<-argmaxJ(log(W.new),Pg,omega,SigmaTilde,phi,n)
    if(is.nan(maxJ[iter])) browser()
    if(iter>1){
      diffJ = (maxJ[iter] - maxJ[iter-1])
      diff.J = c(diff.J,diffJ)
      diff.W = c(diff.W,diffW)
    } 
    W=W.new
  }
  
  LB2=LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)
  #la mise à jour des W est non négociable
  #if(LB2[1]>max(LB0, LB1[1])){
  W=W.new
  diffW=max(abs(F_Sym2Vec(W)-F_Sym2Vec(Wsave)))
  # }else{W=Wsave}
  
  LB2=c(LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p),"W")
  
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nM step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal maxJ difference: ",round(diff.J[iter],4)))
  
  if(plot){
    g=data.frame(Diff.W=diff.W,Diff.J=diff.J, Jbound=maxJ) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+ mytheme.dark+
      facet_wrap(~key, scales="free")+labs(x="iter",y="") 
    print(g)
  }
  
  res=list(W=W, omega=omega, diff=diff, diffW=diffW, finalJ=maxJ[iter], LB=rbind(LB1,LB2))
  
  return(res)
}

#===========
VEMtree<-function(counts,MO,SO,MH,ome_init,W_init,Wg_init, verbatim=TRUE,maxIter=20, 
                  plot=TRUE, eps=1e-2, alpha, vraiOm, print.hist=FALSE, filterPg=TRUE,
                  filterWg=TRUE){
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO);
  hidden=!is.null(MH)
  
  if(hidden){  
    H=(p+1):ncol(ome_init);r=length(H)
    SH=matrix(1,nrow=n,ncol=r)
  }else{
    r=0
    SH= NULL
  }
  
  omega=ome_init;  W=W_init;  Wg=Wg_init
  iter=0 ; lowbound=list()
  KL=c() ; J=c();diffW=c();diffOm=c();diffWg=c();diffPg=c();diffMH=c();diffOmDiag=c();rvalue=c();diffquantile=c()
  diffWiter=1
  t1=Sys.time()
  Pg=matrix(0.5, ncol(W),ncol(W))
  
  
  while((((diffW[iter] > eps) || (diffOm[iter] > eps)) && (iter < maxIter))|| iter<1){ #(diffW[iter] > eps && diffOm[iter] > eps && iter < maxIter) || iter<1
    
    iter=iter+1 
    if(verbatim) cat(paste0("\n Iter n°", iter))
    #VE
    
    resVE<-VE(MO,SO,SH,omega,W,Wg,MH=MH,Pg=Pg,maxIter=1,minIter=1,eps=1e-3, plot=FALSE, 
              alpha=alpha, verbatim=FALSE, hist=print.hist,filterPg=filterPg, filterWg=filterWg)
    KL[iter]=resVE$KL
    M=resVE$Hmeans ; 
    S=resVE$Hvar ;
    
    Pg.new=resVE$Gprobs
    Wg.new=resVE$Gweights
    if(hidden){ 
      SH<-matrix(S[,H],n,r)
      MH.new<-matrix(M[,H],n,r)
      diffMH[iter]<-abs(max(MH.new-MH))
      MH=MH.new
    }
    diffWg[iter]<-abs(max(Wg.new-Wg))
    diffPg[iter]<-abs(max(Pg.new-Pg))
    Wg=Wg.new
    Pg=Pg.new
    
    #M
    resM<-Mstep(M,S,Pg, omega,W,maxIter=2, beta.min=exp(-20),beta.max=exp(20), eps=1e-3 ,plot=FALSE, verbatim=FALSE,
                Wg=Wg, p=p)
    W.new=resM$W
    
    diffW[iter]=abs(max(W.new-W))
    diffWiter=diffW[iter]
    W=W.new
    omega.new=resM$omega
    diffOm[iter]=abs(max(F_Sym2Vec(omega.new)-F_Sym2Vec(omega)))
    diffOmDiag[iter]=abs(max(diag(omega.new)-diag(omega)))
    omega=omega.new
    if(!is.null(vraiOm)){
      diffquantile[iter]=quantile(F_Sym2Vec(omega)[F_Sym2Vec(vraiOm)==1], 0.25) - quantile(F_Sym2Vec(omega)[F_Sym2Vec(vraiOm)==0], 0.75)
    }
    
    lowbound[[iter]] = rbind( resVE$LB, resM$LB)#LowerBound(Pg ,omega, M, S, W, Wg,p)
  }
  
  lowbound=data.frame(do.call(rbind,lowbound))
  lowbound[,1:5]<-apply(lowbound[,1:5],2,function(x) as.numeric(as.character(x)))
  if(hidden){
    features<-data.frame(diffMH=diffMH, diffWg=diffWg, diffPg=diffPg, diffW=diffW, diffOm=diffOm)
  }else{
    features<-data.frame(diffWg=diffWg, diffPg=diffPg, diffW=diffW, diffOm=diffOm)
  }
  
  if(!is.null(vraiOm)){
    features$diffquantile=diffquantile
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
 
    g2<- lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-V6) %>% 
      ggplot(aes(rowid,value, group=key))+geom_point(aes(color=as.factor(V6)), size=3)+geom_line()+
      facet_wrap(~key, scales="free")+
      labs(x="iteration",y="", title="Lower bound and components")+mytheme+
      scale_color_discrete("")
    
    grid.arrange(g1,g2, ncol=1)
  }
  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,omega=omega, lowbound=lowbound, features=features, finalIter=iter, time=time))
}

True_lowBound<-function(Y, M,S,theta,X, W, Wg, Pg, omega){
  p=ncol(Y)
  partJ<-LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p)[1]
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

ICL_T<-function(TrueJ, W,n,r){
  P=EdgeProba(W)
  pen_T=-( sum( P * log(W) ) - log(SumTree(W)) )
  if(pen_T<0) browser()
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
ICL<-function(TrueJ, Pg,W ,S,n,r,d, omega){
  q=ncol(S) ; p=q-r
  H=(p+1):q
  if(r!=0){
    pen_ZH= 0.5*sum(log(S[,H])) + n*r*0.5*(1+log(2*pi))
  }else{
    pen_ZH= 0
  }
  P=Pg
  pen_T=-( sum( P * log(W) ) - log(SumTree(W)) )
  pen_r<-p*(d) + (p*(p+1)/2 +r*p+r)+(q*(q-1)/2 - 1) #d comprends l'intercept
  norm=n*q*(q-1)/2
  #browser()
  ICL=TrueJ-norm  - (pen_T + pen_ZH+pen_r*log(n)/2)
  return(ICL)
}

criteria<-function(List.vem,counts,theta, matcovar,r){
  n=nrow(counts);p= ncol(counts)
  data<-lapply(List.vem,function(vem){
    J<-True_lowBound(counts,vem$M,vem$S, theta, matcovar,vem$W, vem$Wg, vem$Pg, vem$omega )
    vBIC<-VBIC(J,p,r=r, d=ncol(matcovar), n=n)
    ICLT<-ICL_T(J, vem$W,n,r)
    ICLZH<-ICL_ZH(J,vem$S, n,r)
    ICL<-ICL(J, vem$Pg,vem$W,vem$S, n,r,d=ncol(matcovar), omega=vem$omega)
    res=data.frame(vBIC=vBIC, ICL=ICL,J)
    return(res)
  })
  res=do.call(rbind,data)
  res$r=r
  return(res)
}
FitSparsePCA <- function(Y, r=1, alphaGrid=10^(seq(-4, 0, by=.1))){
  # Fit sparse PCA on a grid of alpha
  # Needs an estimate o Sigma: empirical variance of the rotated scores 
  #    + diagonal risidual variances (why not?)
  n <- nrow(Y); p <- ncol(Y); alphaNb <- length(alphaGrid); 
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
  return(list(sPcaOpt=sPcaOpt, alphaGrid=alphaGrid, alphaOpt=alphaOpt, 
              loglik=loglik, bic=bic, cliques=cliques))
}

boot_FitSparsePCA<-function(Y, B,r){
  cliqueList<-list() ; a=1
  lapply(1:B, function(x){
    n=nrow(Y); v=0.8; n.sample=round(0.8*n, 0)
    ech=sample(1:n,n.sample,replace = FALSE)
    Y.sample=Y[ech,]
    c=FitSparsePCA(Y.sample,r=r)$cliques
    if(length(unique(c))==r){
      cliqueList[[a]] <<- c
      a<<-a+1
    }
  })
  cliqueList<-unique(cliqueList)
  return(cliqueList)
}



List.VEM<-function(cliqueList, counts, sigma_obs, MO,SO,alpha,r, cores,maxIter,eps){
  p=ncol(counts) ; O=1:p
  list<-mclapply(cliqueList, function(c){
    #    browser()
    #init
    init=initVEM(counts = counts, initviasigma=c, sigma_obs,r = r)
    Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
    
    #run
    # alpha<-tryCatch(expr={computeAlpha(omegainit[O,O],default =0.3, MO, SO, plot=FALSE)},
    #                 error=function(e){message("sythetic alpha")
    #                   return(0.3)})
    VEM<-VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, 
                 eps=eps, alpha=alpha,verbatim = FALSE,
                 maxIter=maxIter, plot=FALSE,vraiOm=NULL, print.hist=FALSE, filterPg=TRUE,
                 filterWg = TRUE)
    
    return(VEM)
  }, mc.cores=cores)
  converged<-do.call(rbind,lapply(list,function(vem){
    diffW=vem$features$diffW
    conv=(diffW[length(diffW)]<=eps)}))
  if(sum(converged)!=0){
    list=list[converged]
  } 
  return(list)
}
