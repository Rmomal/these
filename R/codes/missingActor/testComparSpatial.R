# DATA
nbgraph=10
valeur=20
type="erdos"
variable="d"

n=100

param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))
Z=rmvnorm(n,mean=rep(0, valeur), sigma=param$sigma)
Zmissing=Z[,-20]

# test missing variable
litree<-function(nb.missing.var, covX,n){
  cliques=findCliques(covX,nb.missing.var+1)[1:nb.missing.var]
  initial.param<-initEM(covX,n,cliquelist = cliques)

  res<-treeAgr.EM(S = covX,k=nb.missing.var, K0 = initial.param$K0,
                  Sigma0 = initial.param$Sigma0,
                  pii=0.5, n, max.iter = 100)
  return(res)
}
res=litree(1,cov(Zmissing),n)
ggplot(data.frame(lik=res$likelihoods,ind=1:length(res$likelihoods)), aes(ind,lik))+geom_point()

# comparaison sigma de res sur spatial et de sigma z originale

cov_matern <- function(distance, range = 1, sill = 1) {
  sill * sill * (1 + distance * sqrt(3) / range) * exp(-distance * sqrt(3) / range)
}
gener_Delta<-function(n){
  sqrt_n <- sqrt(n)
  grid <- data.frame(lon = rep(seq(0, 1, length.out = sqrt_n), times = sqrt_n),
                     lat = rep(seq(0, 1, length.out = sqrt_n), each = sqrt_n)
  )
  Delta <- as.matrix(cov_matern(distance = dist(grid, method = "euclidean")))
  diag(Delta) <- 1
  return(Delta)
}

Delta<-gener_Delta(n)
Epsi<- t(rmvnorm(1, rep(0, nrow(Delta)), Delta))

Zbind<-cbind(Z,Epsi)
Zadd<-Z+Epsi%*%rep(1,valeur)

#rnorm(valeur)

sigmaBind=cov(Zbind)
sigmaAdd=cov(Zadd)
plot(abs(apply(sigmaBind,2,median))) # on s'attend à ce que la covariance de epsi avec Z soit très faible

# on ajoute un effet spatial sur toutes les espèces
# détection effet ?
missingSpatial<-litree(1,sigmaAdd,n)
plot(abs(apply(missingSpatial$Sigma,2,median))) #litree trouve que l'acteur manquant est très peu corrélé
# aux autres données


# meilleure estimation du réseau sur le sigma "corrigé" ?
library(ROCR)
obs<-param$omega
pred_missSp<-missingSpatial$P[1:valeur,1:valeur]
prob_missSp<-missingSpatial$alpha[1:valeur,1:valeur]
Spatial<-litree(0,sigmaAdd,n)
pred_Sp<-Spatial$K[1:valeur,1:valeur]
prob_Sp<-Spatial$alpha[1:valeur,1:valeur]

auc=diagnost_auc(obs,prob_missSp)
auc=diagnost_auc(obs,prob_Sp)

#très bon sur les precisions, très mauvais sur les proba
# sur cet exemple la correction augmente l'auc de 6%


# recuperer le noeud manquant
# espérance conditionnelle
getCondEsp<-function(res,valeur){
  invSxx = solve(res$Sigma[-valeur, -valeur])
  Syx = res$Sigma[valeur, -valeur]
  # Syy = res$Sigma[valeur, valeur];condVar = (Syy - Syx%*%invSxx%*%Syx)[1, 1]
  condEsp = as.vector(Syx%*%invSxx)
  return(condEsp)
}
plotCondEsp<-function(res,valeur, Xmiss, Xobs){

  condEsp=getCondEsp(res, valeur)
  mux=mean(Xobs)
  Xpred = as.vector(condEsp%*%t(Xobs))
  plot(Xpred, Xmiss, pch=16,col="steelblue2",main=round(cor(Xpred, Xmiss), 2)); abline(0, 1)

}
plotCondEsp(missingSpatial,21,Epsi, Zadd)








