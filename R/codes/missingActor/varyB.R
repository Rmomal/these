library(EMtree)
library(PLNmodels)
library(LITree)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(mvtnorm)
library(useful)
library(mclust)
library(MASS)
library(parallel)
library(ROCR)
library(reshape2)#for ggimage
library(gridExtra)
library(harrypotter)
library(sparsepca)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-missing.R")

# simu parameters
set.seed(1) # major and  bad cor
n=200 ;p=14;r=1;type="scale-free";plot=TRUE
O=1:p
################
#----- DATA
# simulate graph and omega, then sigma0 and finally counts
missing_data<-missing_from_scratch(n,p,r,type,plot)
counts=missing_data$Y; sigmaO= missing_data$Sigma;omega=missing_data$Omega
trueClique=missing_data$TC;hidden=missing_data$H
PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  ;SO<-PLNfit$var_par$S  ;sigma_obs=PLNfit$model_par$Sigma
theta=PLNfit$model_par$Theta ;matcovar=matrix(1, n,1)
ome=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
diag(ome)=0

seqB=c(seq(10,200,20), seq(250,500,50), seq(600, 800, 100))

nbcliques=sapply(seqB, function(B){ 
  cat(paste0("B=",B," : "))
  t1<-Sys.time()
  res=length(boot_FitSparsePCA(scale(MO),B,r=1, cores=3)$cliqueList)
  t2<-Sys.time()
  time=difftime(t2, t1)
  cat(paste0(round(time,3), attr(time, "units"),"\n"))
  return(res)})
data.frame(nb=nbcliques, B=seqB) %>% 
  ggplot(aes(B, nb))+geom_point()+mytheme.dark("")

# best VEM  
seqB=c(seq(20,200,20), seq(250,400,50))

cliques_spca=boot_FitSparsePCA(scale(MO),B=800,r=1, cores=3)
List800<-List.VEM(cliquesObj=cliques_spca, counts, sigma_obs,
                  MO,SO,r=1,eps=1e-3, maxIter=200, alpha = 0.1,cores=3,
                  nobeta = FALSE)
J800=c(do.call(rbind, lapply(List800, function(vem){
  tail(vem$lowbound$J,1)
})) )
max.prec=c(do.call(rbind, lapply(List800, function(vem){
 vem$max.prec
})) )
data800<-data.frame(nb_init=1:length(List800), max.prec=max.prec) %>% as_tibble() %>% 
  mutate(maxJ = purrr::map(nb_init, function(x){max(J800[1:x])})) %>% unnest()
plot=data800 %>% ggplot(aes(nb_init, maxJ, color=max.prec))+geom_point()+geom_line()+
  facet_wrap(~ max.prec)+mytheme.dark("")+
  labs(title="max of J among the x first cliques of sPCA with B=800", x="x")
  
ggsave(plot, filename ="/Users/raphaellemomal/these/R/images/J_B800_s1.png", 
       width = 6, height=4 )
JvaryB<-lapply(seqB, function(B){
  cat(paste0("B=",B," : "))
  t1<-Sys.time()
  cliques_spca=boot_FitSparsePCA(scale(MO),B,r=1, cores=3)
  t2<-Sys.time()
  tmp_spca=difftime(t2, t1)
  t3<-Sys.time()
  ListVEM<-List.VEM(cliquesObj=cliques_spca, counts, sigma_obs,
                    MO,SO,r=1,eps=1e-3, maxIter=200, alpha = 0.1,cores=3,
                    nobeta = FALSE)
  t4<-Sys.time()
  tmp_vemtree=difftime(t4,t3)
  J=c(do.call(rbind, lapply(ListVEM, function(vem){
    tail(vem$lowbound$J,1)
  })) )
  max = which.max(J)
  maxJ=J[max]
  Pg=ListVEM[[max]]$Pg
  auc=round(auc(pred = Pg, label = ome),4)
  ppvh=  accppvtpr(Pg,ome,h=15,seuil=0.5)[5]
  time=difftime(t4, t1)
  cat(paste0(round(time,3), attr(time, "units"),"\n"))
  return(list(maxJ=maxJ,auc=auc,ppvh=ppvh, nbcliques=length(cliques_spca$cliqueList), tspca=tmp_spca,
              tvemtree=tmp_vemtree))
}) #calcul en 1h avec la seed 1

saveRDS(JvaryB2, file = "/Users/raphaellemomal/these/R/codes/missingActor/SimResults/JvaryB2.rds")
JvaryB=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/JvaryB.rds")

seqB=c(seq(10,200,20), seq(250,500,50), seq(600, 800, 100))
dataJvaryB=do.call(rbind,(JvaryB)) %>% as_tibble() %>% unnest()
dataJvaryB$B=seqB
dataJvaryB %>%  
  ggplot(aes(B, nbcliques))+geom_point()+mytheme.dark("")

dataJvaryB %>% 
  ggplot(aes(B, maxJ))+geom_point()+mytheme.dark("")

# dataJvaryB %>% 
#   ggplot(aes(B, auc))+geom_point()+mytheme.dark("")
# 
# dataJvaryB %>% 
#   ggplot(aes(B, ppvh))+geom_point()+mytheme.dark("")
