
library(tidyverse)
library(parallel)

data.dir = '/Users/raphaellemomal/these/Data_Oak_remmoa/Oaks-CVacher/'
data.name = 'oaks'
load(paste0(data.dir, data.name, '.RData'))
theme_set(theme_bw())
# Parms
counts = as.matrix(Data$count);
O = Data$offset
covar = as_tibble(Data$covariates)
covar$null=1
n = nrow(counts)
p = ncol(counts)

colnames(counts)<-unlist(lapply(strsplit(colnames(counts),split="_"),function(x){paste0(toupper(x[1]),x[2])}))
colnames(counts)[48]<-"EA"
tested_models=list("null","tree",c("tree","distTOtrunk","distTOground","distTObase"))
models_names=c("Null","Tree","Tree + Distances")

# à refaire oubli offset !!!!!!
compare_output<-ComparEMtree(counts, covar_matrix=covar, models=tested_models, m_names=models_names, 
                             Pt=2/p,  S=100, maxIter=30,O = O,
                             cond.tol=1e-5,cores=3)
compar_graphs(allNets=compare_output,alpha=FALSE,Ft=0.9, nb=6)$G
ggsave(filename = "OaksNets.png",plot = plot,width = 13.6,height = 5.5,
       path = "/Users/raphaellemomal/these/R/images")

# infected samples

infcovar=covar %>% filter(tree==2)
infectedTree=which(covar$tree==2)
infcounts=Y[infectedTree,]   
infO=O[infectedTree,]

vars=c("distTOtrunk","distTOground","distTObase","orientation")
chaine=paste0("~-1+",paste(vars,collapse="+"))
formule=formula(chaine)
matcovar=model.matrix(formule,infcovar)

inftree_covars=ResampleEMtree(infcounts,
               covar_matrix = matcovar,
               O = infO,S =100 , cores=1)
nrow(finalinftree_covs$Pmat)

inftree_covars2=ResampleEMtree(infcounts,
                              covar_matrix = matcovar,
                              O = infO,S1 =101,S2=130 , cores=1)
 
finalinftree_covs=list(Pmat=rbind(inftree_covars$Pmat,inftree_covars2$Pmat),
                       maxIter=c(inftree_covars$maxIter,inftree_covars2$maxIter),
                       times=c(inftree_covars$times,inftree_covars2$times))
saveRDS(finalinftree_covs,file=paste0(data.dir,"infectedTree.rds"))


df<-freq_selec(finalinftree_covs$Pmat,Pt=2/p)

draw_network(df,"infected trees with covariates", layout="circle")



# infected tree without covariates
inftree_null=ResampleEMtree(infcounts,
                              covar_matrix = NULL,
                              O = infO,S =100 , cores=3)
saveRDS(inftree_null,file=paste0(data.dir,"infectedTree_null.rds"))
