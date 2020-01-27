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
hidden=which(diag(omega)%in%sort(diag(omega), decreasing = TRUE)[1:r]) # on cache les r plus gros
trueClique=which(omega[hidden,-hidden]!=0)
if(plot){
  G=draw_network(1*(omega==1),groupes=1*(diag(omega)==diag(omega)[hidden]), 
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
# Tree
Wg_init <- matrix(1, p+r, p+r); diag(Wg_init) = 0; Wg_init =Wg_init / sum(Wg_init)
W_init <- matrix(1, p+r, p+r); diag(W_init) = 0; W_init =W_init / sum(W_init)

# whole Z
initviasigma=init.mclust(sigma_obs,title="Sigma",trueClique = NULL,n.noise=p*3+5)
initial.param<-initEM(sigma_obs,n=n,cliquelist = list(initviasigma),pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
omega_init=initial.param$K0
sigma_init=initial.param$Sigma0


VE<-function(MO,SO,sigma_obs,omega,W,Wg,maxIter,eps, beta.min=1e-6){
  n=nrow(MO)
  p=ncol(MO)
  O=1:ncol(MO)
  H=(p+1):ncol(omega)
  omegaH=omega[H,H]
  KL<-rep(0,maxIter)
  iter=0 ; diffKL=1 ;
  
  while( (diffKL > eps) && (iter < maxIter)){
    iter=iter+1

     #  browser()  #--  Estimates
    Pg = EdgeProba(Wg)
   
    Pghkl= HiddenEdgeProba(Wg,r=1)
    Cg = CgMatrix(Pg,Pghkl,omega,p)
    
    #-- Updates
    # MH et SH
    MH<- (-MO) %*% (Pg[O,H] * omega[O,H]) / omegaH
    SH <- 1/omegaH
    
    M<-cbind(MO, MH)
    S<-cbind(SO, rep(SH,n)) # all SHi have same solution
    # beta g
  
    gamma = optim(par=log(F_Sym2Vec(Wg)),
                  fn=argminKL, gr=Grad_KL_Wg,
                  method='BFGS',
                  Cg, Pg, M,S,omega,W)$par #reminder : besoin d'optim à cause de la contrainte d'identifiabilité des beta
    beta=exp(gamma)

    beta[which(beta< beta.min)] = beta.min
    Wg=F_Vec2Sym(beta)
    # end
    KL[iter]<-argminKL(F_Sym2Vec(log(Wg)), Cg, Pg, M,S,omega,W,trim=TRUE)
    if(iter>1) diffKL = KL[iter] - KL[iter-1]
    print(diffKL)
  }
  Pg = EdgeProba(Wg)
  Pghkl= HiddenEdgeProba(Wg,r=1)
  Cg = CgMatrix(Pg,Pghkl,omega,p)
  KL=KL[1:iter]
  res=list(Gprobs=Pg,Gweights=Wg,Gmeans=M,Gvar=S, Cg=Cg,KL=KL)
  return(res)
}
resVe=VE(MO,SO,sigma_obs,omega_init,W_init,Wg_init,maxIter=10,eps=1e-5)
resVe$KL
resVe$Gprobs

####################
#-----  M steps

Mstep<-function(M,S){
  n=nrow(S)
  sigmaTilde = t(M)%*%M+ diag(colSums(S)) / n
}





