# test fonctions variogram matthieu
#source('krigeage.R')

# Data @ packages
lapply(c("tidyverse", "lubridate", "rgdal", "geoR", "fields", "mvtnorm", "MASS", "reshape", "WhatIf"), library, character.only = TRUE)
theme_set(theme_minimal(base_size = 14))
load(file="/Users/raphaellemomal/these/Data_Oak_remmoa/REMMOA/output/20180620_PLN_REMMOANC.RData")

colnames(X)[7:8]<-c("Xcart","Ycart")
X[,c(6:10,14:20)]<-apply(X[,c(6:10,14:20)], 2, function(x) as.numeric(x))
null_index=which(rowSums(Y)==0)
# Y=Y[-null_index,] %>% as_tibble()
# X=X[-null_index,] %>% as_tibble()
# N=N[-null_index,] %>% as_tibble()

esw=c(250,285,340,300,310,470,270,310,260,170,430,215,210,170,340,230,240,210,145,165,160,200,200,200,200,200,200,200,200,200,200,200)
species=colnames(Y)
detec=data.frame(esw=as.numeric(esw)/1e3,species)

# tests

# N_C_P2 est une région à courts transects : un segment. Comment estimer les variograms alors ?
reg="N_C_O1" # c("N_C_O1","N_C_O2","N_C_O3","N_C_P1","N_C_P2","N_C_P3")
vec_reg<-c("N_C_O1","N_C_O2","N_C_P1","N_C_P2")
vec_reg="N_C_P1"

x="GREPET"
plots_vario<-lapply(vec_reg,function(reg){
  print(reg)
  emp_vario_intra(N[,x],detec$esw[which(species==x)],l=X$Effort,X[,c(7,8)],projected = TRUE,id_transect = X$Transect.Label,
                  plot = TRUE,pairs.min = 20,n_sim = 500,distmat = NULL, coupes = c(-0.01, 0.01,1:4,seq(5,70,5)),esp=x,
                  region=reg, region.label=X$Region.Label, perTransect = TRUE)$g
})
n <- length(plots_vario)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plots_vario, ncol=nCol))

x=vec_esp[1]
emp_vario_intra(Y[,x],detec$esw[which(species==x)],l=X$Effort,X[,c(7,8)],projected = TRUE,id_transect = X$Region.Label,
                plot = TRUE,pairs.min = 20,n_sim = 500,distmat = NULL, breaks = c(-0.01, 0.01,1:4,seq(5,100,5)),esp=x,
                region=NULL, region.label=X$Region.Label)$g

#####################
# tests fonction originale
library(tidyverse)
library(gridExtra)

vario_reg_species<-function(sp, char){
  enq=enquo(sp)
  vec_reg=c(N %>% as_tibble() %>% mutate(region=X$Region.Label) %>% group_by(region) %>% 
              summarise(sum=sum(!!enq)) %>% filter(sum!=0) %>% dplyr::select(region))$region
  
  x=char
  
  plots_vario<-lapply(vec_reg,function(reg){
    emp_vario_intra(Y[,x],detec$esw[which(species==x)],l=X$Effort,X[,c(7,8)],projected = TRUE,
                    id_transect = X$Transect.Label, breaks = c(-0.01, 0.01,1:4,seq(5,150,10)),plot = TRUE,
                    pairs.min = 100, n_sim = 500, distmat = NULL,region=reg, region.label = X$Region.Label,
                    esp=x)$g
    
  })
  n <- length(plots_vario)
  nCol <- floor(sqrt(n))
  do.call("grid.arrange", c(plots_vario, ncol=nCol))
  
}


vario_reg_species(MARTEAU,"MARTEAU")



x="SMADEL"
y=Y[,x]
esw=detec$esw[which(species==x)]
l=X$Effort
x <- rpois(length(y), mean(y))
model_x <- glm(x ~ 1 + offset(2 * esw * l), family = "poisson", control = list(epsilon = 1e-6, maxit = 10000, trace = FALSE))
exp(mean(rnorm(100, as.numeric(coef(model_x)), as.numeric(sqrt(vcov(model_x))))))
mean(exp(rnorm(100, as.numeric(coef(model_x)), as.numeric(sqrt(vcov(model_x))))))

plot(model_x)


