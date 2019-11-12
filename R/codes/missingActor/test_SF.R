# inférer le réseau avec 1 variable manquante et voir si les moyennes correspondent à celles de la température
library(PLNmodels)
library(EMtree)
library(mvtnorm)
library(MASS)
library(saturnin)
library(Matrix)
library(ggplot2)
# DATA
G=generator_graph(p=30,graph="scale-free")
param=generator_param(G=G,signed=TRUE,v=0.01)
hidden=which.max(as.matrix(G))

data_tot=rmvnorm(n=200,sigma=param$sigma)

data_obs=data_tot[,-hidden]
sigma_obs=cov(data_obs)
# LiTree
cliques=findCliques(sigma_obs,1)
initial.param<-initEM(sigma_obs,n=200,cliquelist = cliques,pca=TRUE) # ajout de trois variables manquantes
K0=initial.param$K0
Sigma0=initial.param$Sigma0
infEM=EMtreeMissing(S = sigma_obs,k = 1,method="saturnin",K0 =K0 ,Sigma0 = Sigma0,n=200, XO=data_obs, condMeans=TRUE)

# on récupère les mus pour les plotter contre la variable manquante
res=data.frame(H=data_tot[,hidden],vecmu=t(infEM$mu),vecmu_scale=scale(t(infEM$mu),center = TRUE, scale=TRUE) )

ggplot(res,aes(vecmu,H))+
  geom_point()+theme_light()

#problème d'échelle
# problème de signe quand les données sont signées, forcément...