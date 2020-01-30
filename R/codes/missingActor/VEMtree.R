# VEMtree
library(EMtree)
library(PLNmodels)
library(LITree)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

# simu parameters
set.seed(3)
n=200
p=14
r=1
type="erdos"
plot=TRUE


# Code for one hidden covariate


################
#----- DATA
# simulate graph and omega, then sigma0 and finally counts
data=data_from_scratch(type = type,p = p+r,n = n,signed = FALSE,prob = 5/p,v = 0.001)
omega=data$omega
hidden=which(diag(omega)%in%sort(diag(omega), decreasing = TRUE)[1:r])[1:r] # on cache les r plus gros
trueClique=which(omega[hidden,-hidden]!=0)
if(plot){
  G=draw_network(1*(omega==1),groupes=1*(diag(omega)==diag(omega)[hidden][1]), 
                 layout="nicely",curv=0,nb=2,pal="black",nodes_label = 1:(p+r))$G
  print(G)
}
Kh  <- omega[hidden,hidden]
Ko  <- omega[-hidden,-hidden]
Koh <- omega[-hidden,hidden]
Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
sigmaO=solve(Km)
counts=generator_PLN(sigmaO,covariates = NULL,n=n)

####################
#----- PLN on counts
# optimization of theta and h(Z_O)

PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M # MiO = ith row of MO
SO<-PLNfit$var_par$S # SiO = diag(ith row of SO)
sigma_obs=PLNfit$model_par$Sigma

####################
#-----  VE step
#--  Initialize

# whole Z
initviasigma=init.mclust(sigma_obs,title="Sigma",trueClique = NULL,n.noise=p*3+5)
initial.param<-initEM(sigma_obs,n=n,cliquelist = list(initviasigma),pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
omega_init=initial.param$K0
sigma_init=initial.param$Sigma0

# Tree
Wg_init <- matrix(1, p+r, p+r); diag(Wg_init) = 0; Wg_init =Wg_init / sum(Wg_init)
W_init <- matrix(1, p+r, p+r); diag(W_init) = 0; W_init =W_init / sum(W_init)
W_init[O,O] <- EMtree_corZ(cov2cor(sigma_obs),n = n,maxIter = 20)$edges_weight



VE<-function(MO,SO,sigma_obs,omega,W,Wg,maxIter,eps, beta.min=1e-6, plot=FALSE){
  t1=Sys.time()
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO); H=(p+1):ncol(omega); iter=0 ; diffKL=1 ;
  omegaH=omega[H,H]; MH = matrix(0,n,1)
  logWtree=log(Wg)
  KL<-rep(0,maxIter); diffKL=(100); diff=c(1000); Wtreediff=1; diffW=c(1)
  
  while( (Wtreediff > eps) && (iter < maxIter)){
    iter=iter+1
    #-- Probabilities estimates
    Pg = EdgeProba(exp(logWtree))
    Pghkl= HiddenEdgeProba(exp(logWtree),r=1, verbatim=FALSE)
    Cg = CgMatrix(Pg,Pghkl,omega,p)
    
    #-- Updates
    #- MH et SH
    MH<- (-MO) %*% (Pg[O,H] * omega[O,H]) / omegaH
    SH <- 1/omegaH
    M<-cbind(MO, MH)
    S<-cbind(SO, rep(SH,n)) # all SHi have same solution
    
    #- weights Wg (beta tildes) and final weights (Wtree)
    Mei=Meila(Wg) #justifier que Mei soit calculée avec Wg et pas Wgtree
    lambda=SetLambda(Pg,Mei)
    Wg= Pg/(Mei+lambda)
    Wg[which(Wg< beta.min)] = beta.min
    logWtree.new<-computeWtree(omega, W, Wg, MH, MO, SO)
    
    Wtreediff=max(abs(F_Sym2Vec(logWtree.new)-F_Sym2Vec(logWtree)))
    logWtree=logWtree.new
    # evaluer gradient en le vecteur donné 
    #regarder les P et les M dans le setlambda
    # test MH ressemble à ZH, est différent de ZH comme le dit SH
    
    #-- end
    KL[iter]<-argminKL(F_Sym2Vec(log(Wg)), Cg, Pg, M,S,omega,W,lambda,trim=TRUE)
    if(iter>1){
      diffKL =abs(KL[iter] - KL[iter-1])
      diff = c(diff,diffKL)
      diffW=c(diffW,Wtreediff)
    } 
  }
  
  Pg = EdgeProba(exp(logWtree))
  Pghkl= HiddenEdgeProba(exp(logWtree),r=1)
  Cg = CgMatrix(Pg,Pghkl,omega,p)
  KL=KL[1:iter]
  t2=Sys.time(); time=t2-t1
  cat(paste0("VE step converged in ",round(time,3), attr(test, "units"),
             "\nFinal log(Wtree) difference: ",round(diffW[iter],4)))
  if(plot){
    g=data.frame(Diff.W=resVe$diffW, PartofKL=resVe$KL) %>% rowid_to_column() %>% 
      pivot_longer(-rowid,names_to="key",values_to = "values") %>% 
      ggplot(aes(rowid,values, color=key))+
      geom_point()+geom_line()+scale_color_brewer(palette="Dark2")+
      facet_wrap(~key, scales="free")+theme_light()+labs(x="iter",y="")+
      theme(strip.background=element_rect(fill="gray50",colour ="gray50"))
    print(g)
  }
  res=list(Gprobs=Pg,Gweights=Wg,Gmeans=M,Gvar=S, Cg=Cg,KL=KL, diff=diff, diffW=diffW)
  return(res)
}
resVe=VE(MO,SO,sigma_obs,omega_init,W_init,Wg_init,maxIter=200,eps=1e-3, plot=TRUE)


####################
#-----  M steps

Mstep<-function(M,S,Pg, omega){
  n=nrow(S)
  SigmaTilde = (t(M)%*%M+ diag(colSums(S)) )/ n
  
  #-- Updates
  # beta
  while(diffW>eps && iter < maxIter){
    Mei=Meila(W)  
    lambda=SetLambda(Pg,Mei)
    W.new= Pg/(Mei+lambda)
    W[which(W< beta.min)] = beta.min
    diffW=max(abs(F_Sym2Vec(W.new)-F_Sym2Vec(W)))
    
    W=W.new
  
  
  # omega
  maxi = 1e50
  mini = 1e-50
  
  omegaDiag <- sapply(1:(p+r), function(i){
    dichotomie(mini, maxi, function(omega_ii)
      optimDiag(i, omega_ii, omega, SigmaTilde, Pg), 1e-5)
  })
  
  omega=computeOffDiag(omegaDiag,SigmaTilde)
  diff=c(diff,diffW)
  J=argmaxJ(F_Sym2Vec(log(W)),Pg,omega,sigmaTilde,n)
  #TODO diffJ
  }
  return(W, omega)
}





