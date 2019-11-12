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

#Simu LiTree calcul pedges probas
library(PLNmodels)
library(EMtree)
library(mvtnorm)
library(MASS)
library(saturnin)
library(Matrix)
library(ggplot2)
library(tictoc)
# DATA
p=10
n=100
data=data_from_scratch(type = "erdos",p = p,n = n,signed = FALSE,prob = 5/p)
counts=data$data
hidden=which.max(diag(data$omega))
sigma=solve(data$omega)
counts_obs=counts[,-hidden]
sigma_obs=PLN(counts_obs~1)$model_par$Sigma
trueClique=which(data$omega[hidden,-hidden]!=0)
trueClique[which(trueClique>hidden)]=trueClique[which(trueClique>hidden)]-1

draw_network(data$omega,groupes=1*(diag(data$omega)==diag(data$omega)[hidden]), 
             layout="nicely",curv=0,nb=2,pal="black")
seq=1:p
omega_obs=solve(sigma_obs)
omega_obs[omega_obs<0]=0
draw_network(omega_obs,groupes=1*(1:(p-1) %in% trueClique), 
             layout="nicely",curv=0,nb=2,pal="black")


# LiTree
cliques=findCliques(sigma_obs,1)
list(trueClique)
#verif avec trueClique
initial.param<-initEM(sigma_obs,n=n,cliquelist = list(trueClique),pca=TRUE) # ajout d'une var manquante
K0=initial.param$K0
Sigma0=initial.param$Sigma0

tic()
inftreeAgr= EMtreeMissing(S=sigma_obs, k=1, K0 = K0, Sigma0 = Sigma0, pii=0.5, n=n,
                       max.iter = 200, eps = 0.1,method = "saturnin") #200 iterations 215.03 s n'a pas convergé
toc()

# test à faire pour vérifier que ça marche : treeAfre.EM sans noeud manquant, estce qu'on retrouve le arêtes

infEMtree=EMtree_corZ(CovZ=inftreeAgr$Sigma,n = n)$edges_prob
infEMtreeNaif=EMtree_corZ(CovZ=rbind(cbind(sigma_obs,rep(0,p-1)),c(rep(0,p-1),1)),n = n)$edges_prob
probdata=data.frame("sat"=F_Sym2Vec(inftreeAgr$alpha),"EMtree"=F_Sym2Vec(infEMtree),
                    "EMtreeNaif"=F_Sym2Vec(infEMtreeNaif))

probdata %>%as_tibble() %>% 
  ggplot(aes(EMtree,EMtreeNaif))+geom_point()+theme_light()

fun.auc.ggplot(Listpred = list( EMtreeNaif=infEMtreeNaif,EMtree=infEMtree,
                                saturnin=inftreeAgr$alpha),obs=data$omega, title="")
plot(inftreeAgr$likelihoods)
