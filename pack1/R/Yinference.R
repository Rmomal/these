# Inference in the observed layer

library(sna)
library(ape)
library(poilog)
library(tidyverse)

####
# functions
####
# Tree sampling
tree <- function(p) {
  W = matrix(runif(p ^ 2), p, p)
  W = W + t(W)
  Tree = mst(W)
  return(Tree)
}

gplot(Tree, gmode = 'graph', label = 1:p)

# Omega
build_omega <- function(Tree) {
  coef = 1.1
  degree = colSums(Tree)
  Omega = Tree + coef * diag(degree)

  while (min(eigen(Omega)$values) < 0) {
    coef = 1.1 * coef
    Omega = Tree + coef * diag(degree)
  }
  return(Omega)
}

# Parms
params <- function(omega) {
  p<-ncol(omega)
  mu = rep(1, p)
  #browser()
  Sigma = solve(omega)
  Rho = cov2cor(Sigma)
  sigma = sqrt(diag(Sigma))
  Esp.Y = exp(mu + sigma ^ 2 / 2) #aitchison
  Var.Y = Esp.Y + Esp.Y ^ 2 * (exp(sigma ^ 2) - 1)
  Cov.Y = (Esp.Y %o% Esp.Y) * (exp(Sigma) - 1)
  size.Y = Esp.Y ^ 2 / (Var.Y - Esp.Y) # nm of successful trials for nbinom
  return(list(mu=mu,sigma= sigma, Rho=Rho, Esp=Esp.Y, Var= Var.Y, Cov=Cov.Y, size=size.Y))
}

# Une matrice phi par ligne d'observation
methode_phi <- function(Y, mu, sigma, Rho){
  p = length(Y)
  vect_py = sapply(1:p, function(j){dpoilog(Y[j], mu=mu[j], sig=sigma[j])})
  denum<-vect_py %o% vect_py
  matrice_phi<-matrix(0,p,p)
  #upper tri de phi
  sapply(1:(p - 1), function(j) {
    sapply((j + 1):p, function(k) {
      matrice_phi [j,k]<<-  dbipoilog(Y[j],Y[k],mu1 = mu[j],mu2 = mu[k],
                                      sig1 = sigma[j],sig2 = sigma[k],rho = Rho[j, k])
    })
  })
  matrice_phi<-.5*(matrice_phi+ t(matrice_phi))/denum
  return(matrice_phi)
}

####
# inference
####

build_param <- function(Tree,p) {
  #Tree <- tree(p)
  Omega <- build_omega( Tree)
  params <- params(Omega)
  #Y = rnbinom(p, size = params$size, mu = params$Esp)
  return(list( mu=params$mu,sigma=params$sigma, rho=params$Rho))
}
 matrice_phi<-function(Y,Tree){
  simul<-build_param(matrix(as.numeric(Tree),ncol(Y),ncol(Y)),ncol(Y))
  mu<-simul$mu
  sigma<-simul$sigma
  Rho<-simul$rho
  phi<-lapply(as.list(as_tibble(t(Y))),function(x) methode_phi(x,mu,sigma,Rho))
  return(phi)
 }
 matrice_phi_fromSigma<-function(Y,Sigma){
   simul<-build_param(matrix(as.numeric(Tree),ncol(Y),ncol(Y)),ncol(Y))
   mu<-simul$mu
   sigma<-simul$sigma
   Rho<-simul$rho
   phi<-lapply(as.list(as_tibble(t(Y))),function(x) methode_phi(x,mu,sigma,Rho))
   return(phi)
 }

TreeGGMpoisson<-function(phi){
  p<-ncol(phi)
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)
  #browser()
  FitEM = FitBetaStatic(beta.init=beta.unif, phi=phi,print=TRUE)
  return(FitEM)
}

###
# RUN
###
# Pour inférer le réseau Y à partir de données simulées.
# Pas de VEM car pas d'estimation des paramètres (paramètres tirés en amont pour pouvoir
# simuler)
# Y sont les comptages PLN. Ils sont simulés par metropolis hastings après avoir tiré l'arbre
#Tree . Y et Tree sont stockés; les paramètres sont regénérés à partir de Omega, qui n'a pas de
# raison de différer de l'omega utilisé lors de la simulation (mais ça peut être arrangé).
#
#  On charge les données simulées
# On créer les matrices phi pour chaque ligne d'observtaion des données
# A voir ce qu'on fait de toutes ces matrices : médiane, moyenne artihmétique / géométrique

# En dernière étape du VEM : oui, ici on est en données simulées donc on connait déjà
# les paramètres sigma et coef de régression (non utilisés ici).

Y<-readRDS("/home/momal/Git/these/pack1/R/Simu/PLNTree/Ysample.rds")
Tree<-readRDS("/home/momal/Git/these/pack1/R/Simu/PLNTree/Tree.rds")
phi<-matrice_phi(Y,Tree) # liste de matrices ncol*ncol et de taille nrow(Y)
somme=0*phi[[1]]
mean_phi<-Reduce("+",phi)/length(phi)
library(abind)
a <- do.call(abind, c(phi, list(along=3)))
median_phi<-apply(a, 1:2, median)

FitEM<-TreeGGMpoisson(mean_phi)
FitEM2<-TreeGGMpoisson(median_phi)

####
# Diagnostics
####
heatmap(Tree,Rowv=NA,Colv=NA,scale="none",main="Tree")
heatmap(FitEM$beta,Rowv=NA,Colv=NA,scale="none",main="mean")
heatmap(FitEM2$beta,Rowv=NA,Colv=NA,scale="none",main="median")

a1<-fun.auc.ggplot(FitEM$beta,Tree,"mean",seq(0,0.04,0.005))
a2<-fun.auc.ggplot(FitEM2$beta,Tree,"median",seq(0,0.04,0.005))
matrix<-FitEM$beta
heat<-function(matrix,title,legend=TRUE){
  fortile<-as_tibble(matrix) %>%
  rowid_to_column() %>%
  gather(key,value,-rowid)
position<-ifelse(legend,"right","none")
ggplot(fortile, aes(rowid,key, fill = value)) + geom_raster()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position=position)+
  labs(x = "", y = "",title=title)
}
h1<-heat(matrix(as.numeric(Tree),ncol(Y),ncol(Y)),"Original \ntree",FALSE)
h2<-heat(FitEM$beta,"mean")
h3<-heat(FitEM2$beta,"median")

library(cowplot)
plot_grid(h1,h2,h3, labels=c("Original \ntree","mean","median"), ncol = 1, nrow = 3)
ggdraw() +
  draw_plot(h1, .25, .7, .5, .3) +
  draw_plot(h2, 0, .4, .5, .3) +
  draw_plot(h3, .5, .4, .5, .3) +
  draw_plot(a1, 0, 0, .5, .4) +
  draw_plot(a2, .5, 0, .5, .4)


########
# petit test de performance de l'inverse simple de la matrice de precision de PLN
# vs inference par mixTreeGGM

# direct
library(PLNmodels)
PLN<-PLN(as.matrix(Y) ~ -1 )
sigmaZ<-PLN$model_par$Sigma
omegaZ<-abs(solve(sigmaZ))
fun.auc.ggplot(omegaZ,Tree,"inverse direct",seq(0.06,0.7,0.001))

# notre modele
phi<-1/sqrt(1 - cov2cor(sigmaZ)^2)
diag(phi)<-0
FitZ<-TreeGGMpoisson(phi)
fun.auc.ggplot(FitZ$beta,Tree,"EM",seq(0,0.7,0.001))
heatmap(Tree,Rowv=NA,Colv=NA,scale="none",main="Tree")
heatmap(FitZ$beta,Rowv=NA,Colv=NA,scale="none",main="mean")


# dans les Y
fun.auc.ggplot(FitEM$beta,Tree,"mean",seq(0,0.04,0.005))
