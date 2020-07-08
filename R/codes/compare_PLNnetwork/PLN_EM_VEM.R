# compare PLN EMtree and VEMtree is smth I wanna do for a long time.
# LEt's see what happens after 1h coding this.
# start 10:10

library(tidygraph)
library(ggraph)
library(tidyverse)
library(EMtree)
library(Matrix)
library(mvtnorm)
library(reshape2)
library(PLNmodels)
library(gridExtra)
# simulated cluster data with r=0
set.seed(4) ; n=200 ;p=15;r=0;type="cluster";plot=TRUE ;O=1:p
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y
sigmaO= missing_data$Sigma
upsilon=missing_data$Upsilon
PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma
G=(missing_data$G); 
g1<-ggimage(G)+labs(title="Truth")
g2<-draw_network(G, curv=0, layout="kk")$G
#-- normalize the PLN outputs
D=diag(sigma_obs)
matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
MO=MO*matsig
SO=SO*matsig^2

# RUN PLNnetwork
network_models<-PLNnetwork(counts~1)
PLNnet.path=lapply(network_models$penalties, function(pen){
  model_pen <- getModel(network_models, pen) 
  pen*(model_pen$latent_network()!=0)  
})
K.score <- Reduce("+",PLNnet.path)
scores<- as.matrix(K.score / max(K.score))
diag(scores)=0
g3<-ggimage(scores)+labs(title=paste0("PLNnetwork: AUC=",auc(pred =scores,label = G )))

# RUN EMtree
EMfit=EMtree(PLNfit)
g4<-ggimage(EMfit$edges_prob)+labs(title=paste0("EMtree: AUC=",auc(pred =EMfit$edges_prob,label = G )))

#RUN VEMtree
init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
VEMfit<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.1,
                  maxIter=30, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)
g5<-ggimage(VEMfit$Pg)+labs(title=paste0("VEMtree: AUC=",auc(pred =VEMfit$Pg,label = G )))

 
 
my_layout <-matrix(c(6,1,1,2,2,2,6,1,1,2,2,2,3,3,4,4,5,5,3,3,4,4,5,5), ncol=6, nrow=4, byrow = TRUE)

grid.arrange(g1, g2, g3, g4, g5, ncol=2,layout_matrix=my_layout)


