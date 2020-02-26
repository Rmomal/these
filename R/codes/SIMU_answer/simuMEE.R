# new simu answer MEE
library(ecoCopula)
library(MRFcov)
library(PLNmodels)
library(EMtree)
library(tidyverse)
library(rags2ridges)
theme_set(theme_bw())
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
vignette("Gaussian_Poisson_CRFs")
TPFN_compute(methods = c("ecoCopula","EMtree"),diffs=c("easy","hard"),
             types="erdos",B=30,S=30, cores=3)
build_TFPN_plots(type = "erdos",difficulty = "easy",FDR = TRUE)
E<-gridLine("erdos",TRUE)




type="erdos"
difficulty="easy"
nbgraph=1
dat<-readRDS(paste0(path,"TPFN/Data/dataTPFN_",type,"_",difficulty,nbgraph,".rds"))
Y<-dat[[1]]
edgesOrigin<-ifelse(abs(F_Sym2Vec(dat[[2]]))<1e-16,0,1)
covar <- readRDS(paste0(path,"TPFN/Data/covar_",difficulty,".rds"))
m<- model.matrix(~X1+X2+X3,covar)# on choisi de mettre la constante dans PLN
p=ncol(Y)
##
p=20
G=generator_graph(p=p,graph="tree")
edgesOrigin=F_Sym2Vec(as.matrix(G))
sigma=generator_param(G)$sigma
Y=generator_PLN(as.matrix(sigma))
m=rep(1,nrow(Y))
###-- ecoCopula --###
my_mod=manyglm(Y~1, family="negativ.binomial")
ecoinf =cgr(my_mod)
ecoinf=data.frame(prec=F_Sym2Vec(ecoinf$best_graph$prec),graph=F_Sym2Vec(ecoinf$best_graph$graph),
                  origine=edgesOrigin,cov=F_Sym2Vec(ecoinf$best_graph$cov))
ggplot(ecoinf, aes(y=cov, x=as.factor(origine)))+geom_boxplot()
table(ecoinf$graph, ecoinf$origine)



###-- MRFcov --###
# CRFmodstd <- MRFcov(data = cbind(Y,scale(m[,-1],scale=FALSE)), n_nodes = p, family = 'poisson')

CRFmod <- MRFcov(data = Y, n_nodes = 20, family = 'poisson', symmetrise = "min")

mrfinf=data.frame(inf=abs(F_Sym2Vec(CRFmod$graph)), origine=edgesOrigin)
mrfinf=mrfinf %>% mutate(binary = ifelse(inf!=0,1,0))
table(mrfinf$binary,mrfinf$origine)
ggplot(mrfinf, aes(y=inf, x=as.factor(origine)))+geom_boxplot()+theme_bw()

###-- PartCorr --###

par_corr_matrices <- lapply(seq_len(500), function(x){
  shuffle_data = data.frame(Y) %>%
    dplyr::sample_n(., nrow(.), TRUE)
  
  # regularization penalty increases sparsity, reduces spurious 'interactions'
  par_fit <- -suppressWarnings(rags2ridges::ridgeS(cov(shuffle_data), 
                                                   lambda = 1 / (5 * nrow(shuffle_data))))
})
par_corr_matrix <- Reduce(`+`,par_corr_matrices) / length(par_corr_matrices)
diag(par_corr_matrix) <- 0
rownames(par_corr_matrix) <- colnames(par_corr_matrix) <- colnames(Y)[1:p]

parcor=data.frame(inf=(F_Sym2Vec(par_corr_matrix)), origine=edgesOrigin)

ggplot(parcor, aes(y=inf, x=as.factor(origine)))+geom_beeswarm()

###-- EMtree --###
EMinf=EMtree(PLN(Y~1))$edges_prob
EMinf=data.frame(inf=F_Sym2Vec(EMinf), origine=edgesOrigin) 
EMinf=EMinf %>% mutate(binary = ifelse(inf>2/20,1,0))
table(EMinf$binary,EMinf$origine)
ggplot(EMinf, aes(y=inf, x=as.factor(origine)))+geom_boxplot()








