library(EMtree)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(useful)
library(parallel)
library(ROCR)
library(ggridges)
library(reshape2)#for ggimage
library(gridExtra)
library(ggridges)
library(xtable)
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-V2.R") 
source("/Users/raphaellemomal/these/R/codes/missingActor/VEM_tools.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/modif_pkg.R")
source("/Users/raphaellemomal/these/R/codes/missingActor/fonctions-exactDet.R")
 
pal<-c("#008080","#FA7E1E","#D62976","#756AAB") 
mytheme.dark <-function(legend){list= list(theme_light(), #scale_color_hp(legend,option="Ravenclaw", discrete=TRUE, begin=0.1 ), #scale_fill_hp(legend,option="Ravenclaw", discrete=TRUE, begin=0.1 ),
                                           theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                                                 plot.title = element_text(hjust = 0.5)),
                                           scale_color_manual(legend, values=pal),
                                           scale_fill_manual(legend, values=pal))
return(list)}

# data
H=15
simus=list()
N=403
 
simus=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/simus_V4.rds")
badseeds=which(do.call(rbind, lapply(simus, function(seed){
  length(seed$VEM_1)
}))!=12)
if(length(badseeds)!=0) simus=simus[-badseeds]
N=length(simus)
nH=do.call(rbind,lapply(simus, function(seed){
  G=seed$G
  return(sum(G[,15]))
}))
 
Qual1<-do.call(rbind, lapply(simus, function(seed){
  if(length(seed$VEM_1)==12){
      Pg=seed$VEM_1$Pg
  G=seed$G
  auc=round(auc(pred = Pg, label = G),4)
  H=15
  ppvh=accppvtpr(Pg,G,h=H,seuil=0.5)[5]
  tprh=accppvtpr(Pg,G,h=H,seuil=0.5)[8]
  UH=seed$UH
  MH = seed$VEM_1$M[,H]
  Cor.=abs(cor(UH, MH))
  res=data.frame(auc=auc,ppvh=ppvh,tprh=tprh,Cor.=Cor.)
  }else{
    res=data.frame(auc=NaN,ppvh=NaN,tprh=NaN,Cor.=NaN)
  }
  return(res)
})) %>% as_tibble()
Qual1=Qual1 %>% mutate(nH =nH, influence=unlist(purrr::map(nH, function(x){
  if(x<=5) res="Minor"
  if(x>5 & x<=7) res="Medium"
  if(x>7) res="Major"
  return(res)}))) 
Qual1=Qual1[1:400,]
g=Qual1  %>% rename(TPRH=tprh, AUC=auc, PPVH=ppvh) %>%  gather(key, value, -nH,  -influence) %>% 
  ggplot(aes(value, fct_rev(key), color=key, fill=key))+geom_density_ridges(alpha=0.5)+
  facet_grid(~as.factor(influence))+
  scale_fill_manual(breaks = c("AUC", "PPVH", "TPRH","Cor."),  values=pal)+
  scale_color_manual(breaks =c("AUC", "PPVH", "TPRH","Cor."),  values=pal)+
  labs(x="", y="" )+guides(color=FALSE, fill=FALSE)+
  theme_light()+ theme(strip.background=element_rect(fill="gray50",colour ="gray50"))
ggsave(filename = "simu_densities.png", plot = g,
       path ="/Users/raphaellemomal/these/R/images/", width = 7, height = 4)

#---- table
Qual1 %>% dplyr::select(-nH) %>% group_by(influence) %>% summarize(N=n(),mauc=signif(mean(auc),2), sdauc=signif(sd(auc),2),
                                            mppvh=signif(mean(ppvh, na.rm=TRUE),2), sdppvh=signif(sd(ppvh, na.rm=TRUE),2),
                                            mtprh=signif(mean(tprh),2), sdtprh=signif(sd(tprh),2),
                                            mcor=signif(mean(abs(Cor.)),2), sdcor=signif(sd(abs(Cor.)),2)) %>%
  mutate(auc=paste0(mauc," (",sdauc,")"),ppvh=paste0(mppvh," (",sdppvh,")"),
         tprh=paste0(mtprh," (",sdtprh,")"),Cor.=paste0(mcor," (",sdcor,")")) %>% 
  dplyr::select(influence, N, auc, ppvh, tprh, Cor.) %>% xtable()
