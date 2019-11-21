# lectures données genevieve
Sigma_gen <-   read.table("~/simulations/Data_genevieve/Data/erdos/Sigma_small.txt", quote="\"", comment.char="")
adjmat_gen<-   read.table("~/simulations/Data_genevieve/Data/erdos/adjmat_small.txt", quote="\"", comment.char="")
adjmat_cond <- read.table("~/simulations/Data_genevieve/Data/erdos/adjmat_cond_small.txt", quote="\"", comment.char="")
adjmat_marg <- read.table("~/simulations/Data_genevieve/Data/erdos/adjmat_marg_small.txt", quote="\"", comment.char="")

library(EMtree)
library(tidyverse)
library(ggraph)
library(tidygraph)

p=30
adjmat_cond=matrix(unlist(adjmat_cond),ncol(adjmat_cond),ncol(adjmat_cond))
adjmat_gen=matrix(unlist(adjmat_gen),ncol(adjmat_gen),ncol(adjmat_gen))
adjmat_marg=matrix(unlist(adjmat_marg),ncol(adjmat_marg),ncol(adjmat_marg))
Sigma_gen=matrix(unlist(Sigma_gen),ncol(Sigma_gen),ncol(Sigma_gen))
image(solve(Sigma_gen), label=1:30)
draw_network((adjmat_marg), nodes_label = 1:29, layout="nicely", curv=0.01, size=5, pal="black")


K=solve(Sigma_gen)
Koh=K[]