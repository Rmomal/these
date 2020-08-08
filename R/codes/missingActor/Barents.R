load("/Users/raphaellemomal/these/Data_SR/BarentsFish.Rdata")

counts=Data$count
# vec_void=apply(counts, 2, function(col){
#   sum(col==0)/nrow(counts)
# })
# empty_sp=which(vec_void>0.9)
# filtre_counts=counts[,-empty_sp]
X=data.frame(Data$covariates)
Offset=Data$offset#[,-empty_sp]
n=nrow(counts)
p=ncol(counts)

PLNfit<-PLN(counts~1+offset(log(Offset)))
MO<-PLNfit$var_par$M # MiO = ith row of MO
SO<-PLNfit$var_par$S # SiO = diag(ith row of SO)
sigma_obs=PLNfit$model_par$Sigma
#-- normalize the PLN outputs
D=diag(sigma_obs)
matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
MO=MO*matsig
SO=SO*matsig^2

#-------- spca
# cliques_spca<-FitSparsePCA(counts, r=2)$cliques
# complement=lapply(cliques_spca, function(clique){setdiff(1:p,clique)})
# clique=list()
# clique$cliqueList=lapply(c(cliques_spca,complement), function(cl) list(cl))
#-------- spca boot
clique=boot_FitSparsePCA(counts, B=400, r=1,cores = 3)
#-------- blockmodels
# init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0) 
# Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
# resVEM0<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.05,
#                   maxIter=100, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)
# init=initVEM(counts = filtre_counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0) 
# Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
# resVEM0_filtrecounts<- VEMtree(filtre_counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.05,
#                   maxIter=100, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)
# sbm.vem <- BM_bernoulli("SBM_sym",resVEM0_filtrecounts$Pg, plotting="", verbosity=0)
# sbm.vem$estimate()
# clique_filtre=list();c=1
# sapply(2:4, function(k){
#   paramEstimSBMPoisson <- extractParamBM(sbm.vem,k)
#    sapply(1:k, function(z){
#      clique_filtre$cliqueList[[c]]<<-list(which(paramEstimSBMPoisson$Z==z))
#     c<<-c+1
#   })
# })
# clique_filtre$cliqueList=unique(clique_filtre$cliqueList)
#--------- run
# all counts
# spca boot
ListVEM_all_nofiltre_spca200<-List.VEM(cliquesObj =clique, counts, cov2cor(sigma_obs), MO,SO,r=1,alpha=0.05,
                                    eps=1e-2,maxIter=100, cores=3, trackJ=FALSE)
saveRDS(ListVEM_all_nofiltre_spca200,
        file ="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Barents_VEM_200_2.rds" )
ListVEM_all_nofiltre_spca200=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Barents_VEM_200.rds")
vecJ=do.call(rbind, lapply(ListVEM_all_nofiltre_spca200, function(vem){
  if(length(vem)==14){
    J=tail(vem$lowbound$J, 1) 
  }else{
    J=NaN}
}))
VEM_spca_200=ListVEM_all_nofiltre_spca200[[which.max(vecJ)]]
 
mean(do.call(rbind,lapply(ListVEM_all_nofiltre_spca200, function(vem){vem$time})))

init=initVEM(counts = counts,initviasigma=clique$cliqueList[[7]], cov2cor(sigma_obs),MO,r = 1) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
resVEM<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-1, alpha=0.05,
                 maxIter=50, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)

resVEM$features
P1=resVEM$Pg[,p+1]
P2=resVEM$Pg
ggimage(P1)
ggimage(P2)

#----------

MHhat=VEM_spca_200$M[,p+1]
fit=lm(MHhat~X$Temperature)
paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
       "\nCor. Pears. = ",signif(cor(MHhat,X$Temperature),2))


dataMH=data.frame(MH=MHhat, temperature=X$Temperature) 

g=dataMH%>% 
  ggplot(aes(temperature, MH))+geom_point()+theme_light()+
  geom_smooth(method=lm)+labs(x="Temperature", y="Mh")
  annotate("text",x=1,y=max(MHhat),label=stats)

ggsave(filename = "Barents_MH_temp_white.png", plot = g,
       path ="/Users/raphaellemomal/these/R/images/", width = 3, height = 3)

i=i+1
ListVEM_all_nofiltre_spca200[[12]]$lowbound %>% rowid_to_column() %>%  
  ggplot(aes(rowid,J ))+geom_line()+geom_point(aes(color=as.factor(parameter)),size=2, alpha=0.8)+
  labs(x="iteration",y="", title="Lower bound")+theme_light()+ scale_color_manual("",values="#2976d6")+
  guides(color=FALSE)
ListVEM_all_filtreprec_01[[1]]$lowbound  %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>% 
  ggplot(aes(rowid,value, group=key))+geom_line()+geom_point(aes(color=as.factor(parameter)), size=2, alpha=0.8)+
  facet_wrap(~key, scales="free")+  labs(x="sub-iteration",y="", title="Lower bound and components")+mytheme.dark("")

ListVEM_all_nofiltre_spca200[[12]]$features$diffUpsi

#-------------
# les voisins sont corrélés aussi
p=30
 ggimage(VEM_spca_200$Pg)
 lapply(ListVEM_all_nofiltre_spca200, function(x){
   vois<-which(x$Pg[,p+1]>1e-1)
 })
vois<-which(VEM_spca_200$Pg[,p+1]>1e-1)
data_fit=t(apply(VEM_spca_200$M, 2, function(Mk){
  fit=lm(Mk~X$Temperature)
  res=c(adj.r2=summary(fit)$adj.r.squared,c.Pears.=cor(Mk,X$Temperature))
  return(res)
})) %>% as_tibble %>% rowid_to_column() %>% 
  mutate(vois=rowid%in%vois)

data_fit %>% ggplot(aes(vois, abs(c.Pears.), color=vois))+geom_boxplot()+mytheme.dark("")
data_fit %>% group_by(vois) %>% summarise(mean.corP=median(abs(c.Pears.)),
                                          sdcP=sd(abs(c.Pears.))) %>% xtable()

#-------------
g1<-ggimage(VEM_spca_200$Pg)
init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0)
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
resVEM0<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-3, alpha=0.05,
                  maxIter=100, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=FALSE)
g0<-ggimage(resVEM0$Pg)
grid.arrange(g0, g1, ncol=2)

lab0=ifelse(1:30 %in% vois, 1:30, "")
lab1=ifelse(1:31 %in% c(vois,31), c(1:30,"H"), "")


n0=draw_network(resVEM0$Pg, curv=0,nb = 1,pal="#31374f", layout="fr",nodes_label =lab0 ,
             groupes = 1*(1:30 %in% vois))
n0
n1lay=rbind(n0$finallayout[,1:2],c(-8.5,5))
n1lay[30,]=c(-7,2)
n1lay[13,]=c(-10,1)
n1lay[15,]=c(-11.8,4)
n1lay[27,]=c(-13,9)
n1<-draw_network(VEM_spca_200$Pg, curv=0,nb = 1, pal="#31374f",layout=NULL,stored_layout = n1lay, nodes_label =lab1,
                 groupes = c(1*(1:30 %in% vois),2))
n1
g=grid.arrange(n0$G, n1$G, ncol=2)
ggsave(plot=g,filename = "Barents_net_comp3.png", path =  "/Users/raphaellemomal/these/R/images",
       width=7, height=4 )

g1<-ggimage(VEM_spca_200$Pg[c(setdiff(1:30,vois),vois,31),c(setdiff(1:30,vois),vois,31)])+
  labs(title="r = 1")
g0<-ggimage(resVEM0$Pg[c(setdiff(1:30,vois),vois),c(setdiff(1:30,vois),vois)])+labs(title="r = 0")
g=grid.arrange(g0, g1, ncol=2)
ggsave(plot=g,filename = "Barents_mat_comp.png", path =  "/Users/raphaellemomal/these/R/images",
       width=7, height=5 )
