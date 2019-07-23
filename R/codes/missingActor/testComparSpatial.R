# DATA
nbgraph=10
valeur=20
type="cluster"
variable="d"

n=100

param<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds"))
Z=rmvnorm(n,mean=rep(0, valeur), sigma=param$sigma)
Zmissing=Z[,-c(18:20)]

# test missing variable
litree<-function(nb.missing.var, covX,n){
  cliques=findCliques(covX,nb.missing.var+1)[1:nb.missing.var]
  initial.param<-initEM(covX,n,cliquelist = cliques)

  res<-treeAgr.EM(S = covX,k=nb.missing.var, K0 = initial.param$K0,
                  Sigma0 = initial.param$Sigma0,
                  pii=0.5, n, max.iter = 50)
  return(res)
}
res=litree(3,cov(Zmissing),n)
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

sigmaBind=cov(Zbind)
sigmaAdd=cov(Zadd)
plot(abs(apply(sigmaBind,2,median))) # on s'attend à ce que la covariance de epsi avec Z soit très faible

# on ajoute un effet spatial sur toutes les espèces
# détection effet ?
missingSpatial<-litree(1,sigmaAdd,n)
plot(abs(apply(missingSpatial$Sigma,2,median))) #litree trouve que l'acteur manquant est très peu corrélé
# aux autres données

sigma_corr_spatial<-missingSpatial$Sigma[1:valeur,1:valeur]

# meilleure estimation du réseau sur le sigma "corrigé" ?
library(ROCR)
obs<-param$omega
pred_missSp<-missingSpatial$K[1:valeur,1:valeur]
prob_missSp<-missingSpatial$alpha[1:valeur,1:valeur]
Spatial<-litree(0,sigmaAdd,n)
pred_Sp<-Spatial$K[1:valeur,1:valeur]
prob_Sp<-Spatial$alpha[1:valeur,1:valeur]

auc=diagnost_auc(obs,pred_missSp)
auc=diagnost_auc(obs,pred_Sp)

#très bon sur les precisions, très mauvais sur les proba
# sur cet exemple la correction diminue l'auc de 1.5% sur les K, de 3% sur les proba
