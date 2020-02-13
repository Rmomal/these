#############
# initialization
#Find the clique created by the missing actor
findCliques<-function(covX,k){
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  # Find k groups of nodes which are 'clique likes"
  # INPUT
  #   X      : data matrix (nxp)
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
  # browser()
  X<-data.frame(t(t(eigen(S)$vectors[,1:2])*sqrt(eigen(S)$values[1:2])))
  
  p=nrow(X)+nb.missing
  b=apply(X, 2, range)
  poissonNoise<-apply(b, 2, function(x,n=n.noise){
    runif(n,min=x[1]-0.1, max=x[2]+0.1)
  })
  data=rbind(X, poissonNoise)
  noiseInit<-sample(c(T,F), size=nrow(X)+n.noise, replace=T, prob=c(3, 1))
  
  datapolar=cart2pol(x=data[,1],y=data[,2])[,1:2]
  datapolar=datapolar %>% mutate(theta2=ifelse(theta>pi,theta-pi,theta)) %>% dplyr::select(r,theta2)
  newdata=pol2cart(datapolar$r,datapolar$theta2)[,1:2]
  
  clust=Mclust(data=newdata,initialization = list(noise=noiseInit),G=nb.missing)
  groups<-mclust::map(clust$z)
  res<-which(groups==1)[which(groups==1)<=nrow(X)] 
  if(plot){
    if(!is.null(trueClique)){
      #False positives rate
      N=setdiff(1:p,trueClique)
      FP=sum(res%in%N)/length(N)
      # False negatives rate
      FN=sum(setdiff(1:p,res)%in%trueClique)/length(trueClique)
      title=paste0(title,":"," FN=", round(FN,2),",FP=",round(FP,2))
    }
    
    g= ggplot(X,aes(X1,X2, label=rownames(X),color=rownames(X)%in%res))+geom_point(size=0.1)+
      theme_light()+geom_text()+labs(x="eig vect 1",y="eig vect 2", title=title)+
      guides(color=FALSE)+scale_color_brewer(palette="Dark2")+
      geom_hline(yintercept=0, color="gray50")+geom_vline(xintercept=0, color="gray50")
    print(g)
  }
  return(res)
}

#################
# For OPTIM

#=====
# VE step
VE<-function(MO,SO,sigma_obs,omega,W,Wg,MH = matrix(1,n,r),maxIter,minIter,eps, 
             alpha,form,beta.min=1e-10, plot=FALSE,verbatim=FALSE){
  t1=Sys.time()
  
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); H=(p+1):ncol(omega);r=length(H); iter=0 ; 
  omegaH=omega[H,H]; 
  diffKL=-1 ; diff=c(1000); Wdiff=1; diffW=c(0.1); diffMH=1; diffM=c(0.1)
  
  M<-cbind(MO, MH)
  phi=CorOmegaMatrix(omega)
  SH <- 1/omegaH
  S<-cbind(SO, rep(SH,n)) # all SHi have same solution, depending only on Mstep
  # logWtree<-computeWtree(omega, W, Wg, MH, MO, SO,phi,alpha=alpha, trim=FALSE)  
  Wg= computeWg(phi,omega,W,MH,MO,alpha,form=form)
  
  diag(Wg)=1
  KL<-numeric(maxIter)
  
  while( ((diffKL < 0) && (iter < maxIter)) || (iter < minIter)){
    iter=iter+1
    #   cat(iter)
    #-- Probabilities estimates
    Pg = EdgeProba(Wg)
    if(sum(is.nan(Pg))!=0){
      browser()
      cat( ": adjust ", summary(c(logWtree)), "\n")
      logWtree=log(adjustW(exp(logWtree)))
      cat("to: ",summary(c(logWtree)), "\n")
      Pg = EdgeProba(exp(logWtree))
    } 
    #-- Updates
    #- MH et SH
    MH.new<-  (-MO) %*% (Pg[O,H] * omega[O,H]) / omegaH
    diffMH<-max(abs(MH-MH.new))
    
    # M<-cbind(MO, MH.new)
    # Mei=Meila(Wg) #justifier que Mei soit calculée avec Wg et pas Wgtree
    # lambda=SetLambda(Pg,Mei)
    Wg.new= computeWg(phi,omega,W,MH,MO,alpha, form=form) #Pg/(Mei+lambda)
    diag(Wg.new)=1
    Wg.new[which(Wg.new< beta.min)] = beta.min
    Wdiff=max(abs(F_Sym2Vec(Wg.new)-F_Sym2Vec(Wg)))
    if(is.nan(Wdiff)) browser()
    
    KL[iter]<-argminKL(F_Sym2Vec(log(Wg.new)), Pg, M,S,omega,phi,W,p)
    if(iter>1) diffKL =(KL[iter] - KL[iter-1])
    MH=MH.new
    Wg=Wg.new
    #-- end
    if(iter>1){
      diff = c(diff,diffKL)
      diffW=c(diffW,Wdiff)
      diffM=c(diffM,diffMH)
    } 
  }
  KL=KL[1:iter]
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nVE step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal W difference: ",round(diffW[iter],5),
                          "\nFinal KL difference: ",round(diff[iter],4)))
  if(plot){
    g=data.frame(Diff.W=diffW,  diff.KL=diff,diff.MH=diffM, PartofKL=KL) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+scale_color_brewer(palette="Dark2")+
      facet_wrap(~key, scales="free")+theme_light()+labs(x="iter",y="")+
      theme(strip.background=element_rect(fill="gray50",colour ="gray50"))
    print(g)
  }
  Pg = EdgeProba(Wg)
  M=cbind(MO,MH)
  
  res=list(Gprobs=Pg,Gweights=Wg,Hmeans=M,Hvar=S, KL=KL[iter], diff=diff, diffW=diffW )
  return(res)
}

argminKL <- function(gamma, Pg, M,S,omega,phi,W,p ){
  # browser()
  r=ncol(omega)-p
  O = 1:p
  H=(p+1):(p+r)
  omegaH=omega[H,H]
  
  EhZoZo = t(M[,O])%*%M[,O]+ diag(colSums(S[,O]))
  
  KL <- -n*0.5*sum(log(diag(omega)))+0.5*omegaH*( t(M[,H])%*%M[,H]+sum(S[,H]) ) +   
    sum(diag(Pg[O,H]*omega[O,H] %*% (t(M[,H])%*%M[,O]) )) + 
    0.5*sum(diag( (Pg[O,O]*omega[O,O]) %*% EhZoZo)) -
    n*0.5*2*sum(F_Sym2Vec(Pg)*log(F_Sym2Vec(phi)))+
    2*sum(F_Sym2Vec(Pg)*(gamma-log(F_Sym2Vec(W))))-log(SumTree(F_Vec2Sym(exp(gamma))))+
    log(SumTree(W)) -0.5*sum(log(S[,H]))#-lambda*(sum(exp(gamma)) - 0.5)
  
  
  return(KL)
}
computeWg<-function(phi,omega,W,MH,MO, alpha, form="theory"){
  p=ncol(MO) ; r=ncol(MH)
  O = 1:p ; H = (p+1):(p+r)
  psi=phi^(n*alpha*0.5)
  logWg<-matrix(0,(p+r),(p+r))
  Wg<-matrix(0,(p+r),(p+r))
  if(form=="theory"){
    #   browser()
    logWg[O,O]<-log(W[O,O])+log(psi[O,O])-0.5*alpha*omega[O,O]*(t(MO)%*%MO)
    logWg[O,H]<-log(W[O,H])+log(psi[O,H])-alpha*omega[O,H]*(t(MH)%*%MO)
  }
  if(form=="id1"){
    
    logWg[O,O]<-log(W[O,O])+log(psi[O,O])+abs(-0.5*alpha*omega[O,O]*(t(MO)%*%MO))
    logWg[O,H]<-log(W[O,H])+log(psi[O,H])+abs(-alpha*omega[O,H]*(t(MH)%*%MO))
  }
  
  if(form=="id2"){
    logWg[O,O]<- -0.5*alpha*omega[O,O]*(t(MO)%*%MO)-log(W[O,O])-log(psi[O,O])
    logWg[O,H]<- -alpha*omega[O,H]*(t(MH)%*%MO)-log(W[O,H]) -log(psi[O,H])
  }
  if(form=="id3"){
    logWg[O,O]<- abs(-0.5*alpha*omega[O,O]*(t(MO)%*%MO))-log(W[O,O])-log(psi[O,O])
    logWg[O,H]<- abs(-alpha*omega[O,H]*(t(MH)%*%MO))-log(W[O,H])-log(psi[O,H]) 
  }
  if(form=="id4"){
    #  browser()
    logWg[O,O]<-log(W[O,O]) -0.5*alpha*omega[O,O]*(t(MO)%*%MO)-log(psi[O,O])
    logWg[O,H]<-log(W[O,H]) -alpha*omega[O,H]*(t(MH)%*%MO)-log(psi[O,H])
  }
  
  logWg[H,O]<-t(logWg[O,H])
  diag(logWg) = 0
  # shrinking and centering
  gamma=logWg[O,H]
  gamma[which(gamma<(-30))]=-30
  gamma=gamma-mean(gamma)
  gamma[which(gamma>10)]=10
  Wg[O,H]=exp(gamma)
  Wg[H,O]<-t(Wg[O,H])
  gamma=F_Sym2Vec(logWg[O,O])
  gamma[which(gamma<(-30))]=-30
  gamma=gamma-mean(gamma)
  gamma[which(gamma>max(log(Wg[O,H])))]=max(log(Wg[O,H]))
  Wg[O,O]=exp(F_Vec2Sym(gamma))
  
  diag(Wg) = 1
  
  return(Wg)
}
CorOmegaMatrix<-function(omega){ # code for 1 hidden covariate. Otherwise, must be a sum on h
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
Mstep<-function(M,S,Pg, omega,W,maxIter, beta.min, trim=TRUE,plot=FALSE,eps, verbatim=FALSE){
  t1=Sys.time()
  diffJ=(1); diff.J=c(1);  diff.W=c(0.05);diffW=1
  maxJ=c()
  n=nrow(S) 
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  iter=0
  
  while((diffW>eps && iter < maxIter) || iter < 3){
    iter=iter+1
    
    if(verbatim)  cat(paste0("\nIter n°:",iter))
    #-- Updates
    # beta
    
    Mei=Meila(W)  
    # lambda=SetLambda(Pg,Mei)
    W.new= Pg/(Mei)#+lambda)
    diag(W.new)=1
    W.new[which(W.new< beta.min)] = beta.min
    diffW=max(abs(F_Sym2Vec(W.new)-F_Sym2Vec(W)))
    
    # W=trimW(W)
    # omega
    maxi = 1e+3
    mini = 1e-3
    bool=FALSE
    omegaDiag <- sapply(1:(p+r), function(i){
      # cat("\ndich. ",i," : ")
      #    if( iter==1 && i==1) bool=TRUE
      dichotomie(mini, maxi, function(omega_ii)
        optimDiag( omega_ii,i, omega, SigmaTilde, Pg), 1e-5, bool=bool)
    })
    
    omega.new=computeOffDiag(omegaDiag,SigmaTilde)
    phi=CorOmegaMatrix(omega)
    maxJ[iter]<-argmaxJ(log(W.new),Pg,omega.new,SigmaTilde,phi,n)
    if(is.nan(maxJ[iter])) browser()
    if(iter>1){
      diffJ = (maxJ[iter] - maxJ[iter-1])
      diff.J = c(diff.J,diffJ)
      diff.W = c(diff.W,diffW)
    } 
    omega=omega.new
    W=W.new
  }
  t2=Sys.time(); time=t2-t1
  if(verbatim) cat(paste0("\nM step converged in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal maxJ difference: ",round(diff.J[iter],4)))
  # browser()
  if(plot){
    g=data.frame(Diff.W=diff.W,Diff.J=diff.J, Jbound=maxJ) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+scale_color_brewer(palette="Dark2")+
      facet_wrap(~key, scales="free")+theme_light()+labs(x="iter",y="")+
      theme(strip.background=element_rect(fill="gray50",colour ="gray50"))
    print(g)
  }
  res=list(W=W, omega=omega, diff=diff, diffW=diffW, finalJ=maxJ[iter])
  return(res)
}

argmaxJ<-function(gamma,Pg,omega,sigmaTilde,phi,n, trim=TRUE, verbatim=TRUE){
  
  maxJ <- sum(Pg*F_Vec2Sym(F_Sym2Vec(gamma) + n*0.5*F_Sym2Vec(log(phi)))) + 
    n*0.5*sum(log(diag(omega))) - log(SumTree(exp(gamma))) - 
    0.5*n*sum( Pg * omega * sigmaTilde)# - n*0.5*sum(diag(omega*sigmaTilde)) 
  # + lambda*(sum(W)-0.5)
  
  return(maxJ)
}

optimDiag <- function(omega_ii,i,  omega, SigmaTilde, Pg) {
  # browser()
  q <- nrow(omega)
  phi = CorOmegaMatrix(omega)
  newphi=(1- phi)/phi
  diag(newphi)=0
  quotient = diag(Pg %*% newphi)[i]
  
  grad=(1/omega_ii)*(quotient+1) - SigmaTilde[i, i]
  return( grad)
}
dichotomie <- function(a, b, F.x, eps, bool=FALSE){
  if(bool) browser()
  x.min = ifelse(F.x(0) >0,a,1e-4);
  while(F.x(x.min)<0){x.min = x.min -x.min/2}
  x.max = b
  while(F.x(x.max)>0){x.max = x.max * 2}
  x = (x.max+x.min)/2
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  iter=0
  while(abs(x.max-x.min) > eps){
    iter=iter+1
    
    # if(iter%%10 == 0) cat(paste0(iter,", "))
    # if(iter > 30) browser()
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
  #  cat(paste0("f=",f,", omega=",x))
  return(x)
}
computeOffDiag<-function(omegaDiag,SigmaTilde){
  q=length(omegaDiag)
  omega=matrix(0,q,q)
  sapply(1:(q-1),
         function(j){
           sapply((j+1):q,
                  function(k){
                    omega[k, j] <<- (1 - sqrt( 1+4*SigmaTilde[k,j]^2*omegaDiag[j]*omegaDiag[k]))/(2*SigmaTilde[k,j])
                    omega[j, k] <<- omega[k, j]
                  }
           )
         }
  )
  diag(omega)=omegaDiag
  return(omega)
}



#####################
# For double hidden proba

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
EMtree_corZ<-function(CovZ,n,  maxIter=30, cond.tol=1e-10, verbatim=TRUE, plot=FALSE){
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
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  # while(!is.finite(Wcum)){
  #   #handles numerical issues with matrix tree theorem
  #   it=it+1
  #   borne=30-it
  #   if(verbatim) message(cat("W corrected, bound=",borne))
  #   
  #   W.log=log(F_Sym2Vec(W))
  #   W.center=W.log-mean(W.log)
  #   W.center[which(W.center<(-borne))]=-borne
  #   W=F_Vec2Sym(exp(W.center))
  #   Wcum = SumTree(W)
  # }
  
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
  return(Y)
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
data_from_scratch<-function(type, p=20,n=50, r=5, covariates=NULL,
                            prob=log(p)/p,dens=log(p)/p,
                            signed=FALSE,v=0,draw=FALSE){
  # make graph
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

################
# diagnostic functions
auc<-function(pred,obs){
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}
accppvtpr<-function(probs,omega,h, seuil=0.5){
  Acc=round(mean(1*(probs>seuil)==omega),2) #(TP+TN)/(TP+FP+TN+FN)
  AccH=round(mean(1*(probs[h,]>seuil)==omega[h,]),2)
  AccO=round(mean(1*(probs[-h,-h]>seuil)==omega[-h,-h]),2)
  PPV=round(sum((omega!=0)*(probs>seuil))/(sum((omega!=0)*(probs>seuil))+
                                             sum((omega==0)*(probs>seuil))),2)#TP/(TP+FP)
  PPVH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/(sum((omega[h,]!=0)*(probs[h,]>seuil))+
                                                      sum((omega[h,]==0)*(probs[h,]>seuil))),2)
  PPVO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))+
                                                            sum((omega[-h,-h]==0)*(probs[-h,-h]>seuil))),2)
  
  TPR=round(sum((omega!=0)*(probs>seuil))/sum(omega!=0), 2)
  TPRH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/sum(omega[h,]!=0), 2)
  TPRO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/sum(omega[-h,-h]!=0), 2)
  
  return(c(Acc,AccH,AccO,PPV,PPVH,PPVO,TPR,TPRH,TPRO))
}
ggimage<-function(data){
  melted_data <- melt(data)
  ggplot(melted_data, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +guides(fill=FALSE)+ theme(plot.title = element_text(size=10, hjust=0.5))
}
plotVEM<-function(probs,omega,r,seuil){
  q=ncol(omega)
  h=(q-r):q
  performance=accppvtpr(probs,omega,h,seuil)
  Acc=performance[1] ;AccH=performance[2] ;AccO=performance[3] 
  PPV=performance[4] ;PPVH=performance[5] ; PPVO=performance[6]
  TPR=performance[7] ;TPRH=performance[8] ;TPRO=performance[9] 
  p1<-ggimage(resVe$Gprobs)+labs(title=paste0("G hat (thresh=",seuil,")"))
  p2<-ggimage(ome)+labs(title="G")
  grid.arrange(p1,p2,ncol=2, top=paste0("Tpr=",TPR," (TprO=",TPRO," , TprH=",TPRH,
                                        ")\n Ppv=",PPV," (PpvO=",PPVO," , PpvH=",PPVH,")"))
  
}

plotM<-function(resM,ome,h,seuil){
  performance=accppvtpr(resVe$Gprobs,ome,h,seuil)
  Acc=performance[1] ;AccH=performance[2] ;AccO=performance[3] 
  PPV=performance[4] ;PPVH=performance[5] ; PPVO=performance[6]
  TPR=performance[7] ;TPRH=performance[8] ;TPRO=performance[9] 
  p1<-ggimage(resVe$Gprobs)+labs(title=paste0("G hat (thresh=",seuil,")"))
  p2<-ggimage(ome)+labs(title="G")
  grid.arrange(p1,p2,ncol=2, top=paste0("Tpr=",TPR," (TprO=",TPRO," , TprH=",TPRH,
                                        ")\n Ppv=",PPV," (PpvO=",PPVO," , PpvH=",PPVH,")"))
  
}


courbes_seuil<-function(probs,omega,h,seq_seuil){

  tmp=sapply(seq_seuil,function(x)  accppvtpr(seuil=x,probs=probs,omega=omega,h=h))
  res=data.frame(cbind(t(tmp),seq_seuil))
  colnames(res)=c("Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","seuil")
  return(res)
  
}
