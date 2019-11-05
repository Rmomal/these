library(ade4)
library(tidyverse)
library(parallel)
data(baran95)
counts = as.matrix(baran95$fau)
covar = as_tibble(baran95$plan)
covar$null=1
n = nrow(counts)
p = ncol(counts)


tested_models=list("null","date","site",c("date","site"))
models_names=c("Null","Date","Site","Date + Site")
compare_output<-ComparEMtree(counts, covar_matrix=covar, models=tested_models, m_names=models_names, 
                             Pt=2/p,  S=100, maxIter=30,
                             cond.tol=1e-8,cores=3)
plot=compar_graphs(allNets=compare_output,alpha=FALSE,Ft=0.9, nb=4)$G
ggsave(filename = "BaransNets.png",plot = plot,width = 13.6,height = 4.4,
       path = "/Users/raphaellemomal/these/R/images")
