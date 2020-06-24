library(ade4)
data(baran95)
 counts=as.matrix(baran95$fau)
# poor_species<-which(colSums(raw_counts==0)/nrow(raw_counts)>0.9)
# counts=as.matrix(baran95$fau[,-poor_species])

# Lancer VEM avec 2 manquants et représenter les sites sur le plan des 2 vecteurs de moyenne 
# initialization
p=ncol(counts) ; n=nrow(counts) ; B=100 ; O=1:p
PLNfit<-PLN(counts~1)
MO<-PLNfit$var_par$M  
SO<-PLNfit$var_par$S  
sigma_obs=PLNfit$model_par$Sigma
theta=PLNfit$model_par$Theta 
matcovar=matrix(1, n,1)
#-- normalize the PLN outputs
D=diag(sigma_obs)
matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
MO=MO*matsig
SO=SO*matsig^2
# no missing actor
init0=initVEM(counts , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
VEM<-VEMtree(counts, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=1e-3,
                       Wg_init =Wginit,plot = FALSE, maxIter = 100,print.hist = FALSE,
                       alpha=0.05, verbatim=TRUE, trackJ=FALSE )

#- cliques
# cliques_spca<-(FitSparsePCA(counts, r=2)$cliques)
# complement=(lapply(cliques_spca, function(clique){setdiff(1:p,clique)}))
# clique=list()
# clique$cliqueList=c(list(cliques_spca),list(complement))
cliques_spca=boot_FitSparsePCA(counts, B=200,r=2)
# ListVEM_4<-List.VEM(cliquesObj =clique, counts, cov2cor(sigma_obs), MO,SO,r=2,alpha=0.05,
#                   eps=1e-2,maxIter=100, cores=3, trackJ=FALSE)
# ListVEM_boot<-List.VEM(cliquesObj =cliques_spca, counts, cov2cor(sigma_obs), MO,SO,r=2,alpha=0.05,
#                   eps=1e-3,maxIter=100, cores=3, trackJ=FALSE)
ListVEM_boot_200<-List.VEM(cliquesObj =cliques_spca, counts, cov2cor(sigma_obs), MO,SO,r=2,alpha=0.05,
                       eps=1e-3,maxIter=100, cores=3, trackJ=FALSE)
saveRDS(ListVEM_boot_200,
        file ="/Users/raphaellemomal/these/R/codes/missingActor/SimResults/Fatala_VEM_200.rds" )
vec_sdMH=do.call(rbind, lapply(ListVEM_boot_200, function(vem){
  if(length(vem)==14){
   # J=tail(vem$lowbound$J, 1)
 
    sdM=apply(vem$M,2, function(x) sd(x))
    sdMH=tail(sdM,2)
  }else{
   # J=NaN
    sdMH = c(NaN,NaN)}
}))
vec_J=do.call(rbind, lapply(ListVEM_boot_200, function(vem){
  if(length(vem)==14){   J=tail(vem$lowbound$J, 1)
  }else{ J=NaN}
}))

dataFatala=data.frame(vec_sdMH,vec_J)
colnames(dataFatala) = c("sd1","sd2","J")
dataFatala=dataFatala %>% as_tibble() %>% 
  mutate(num=1:nrow(dataFatala), TwoActors=ifelse(log(sd1)>-15 & log(sd2) >-15,"on", "off"))
dataFatala %>% gather(key, value, -J,-TwoActors) %>% 
  ggplot(aes(J, value, color=key))+geom_point()+mytheme.dark("")+facet_wrap(~key)
dataFatala %>% gather(key, value, -J,-TwoActors,-num) %>% 
  ggplot(aes(as.factor(J), log(as.numeric(value)), color=key, group=interaction(J,key)))+
  geom_point(position = position_dodge(width=0.5))+mytheme.dark("") 
dataFatala %>% gather(key, value,-J,-num,-TwoActors) %>% 
  ggplot(aes(log(value), color=key, fill=key))+
  geom_histogram(alpha=0.5,bins = 30)+mytheme.dark("")+facet_wrap(~key)
 
dataFatala %>%  ggplot(aes(log(sd1), log(sd2), color=J))+
  geom_point(position=position_dodge(0.75))+theme_light()  

dataFatala %>%filter(log(sd1)>-20,log(sd2)>-20) %>%  ggplot(aes(log(sd1), log(sd2), color=J))+
  geom_point()+theme_light()  

max(dataFatala$J, na.rm=TRUE)
dataFatala %>%filter(log(sd1)>-25,log(sd2)>-25) %>% 
  filter(J==max(J))


dataFatala  %>% filter(!is.na(J)) %>% 
  ggplot(aes(TwoActors, J))+geom_point()+mytheme.dark("")

dataFatala %>% filter(TwoActors=="on") %>% filter(J==max(J))

VEM_fatala_200=ListVEM_boot_200[[32]]
VEM_fatala_4=ListVEM_4[[which.max(vecJ)]]
VEM_fatala_boot=ListVEM_boot[[which.max(vecJ)]]
ggimage(VEM_fatala_4$Pg)
ListVEM[[4]]$lowbound%>% rowid_to_column() %>%  
  ggplot(aes(rowid,J ))+geom_line()+geom_point(aes(color=as.factor(parameter)),size=2, alpha=0.8)+
  labs(x="iteration",y="", title="Lower bound")+mytheme+ scale_color_manual("",values="#2976d6")+
  guides(color=FALSE)
VEM_fatala$features$diffUpsi
init=initVEM(counts = counts,initviasigma=clique$cliqueList[[1]], cov2cor(sigma_obs),MO,r = 2) 
Wginit= init$Wginit; Winit= init$Winit; upsinit=init$upsinit ; MHinit=init$MHinit
resVEM<- VEMtree(counts,MO,SO,MH=MHinit,upsinit,Winit,Wginit, eps=1e-2, alpha=0.05,
                 maxIter=30, plot=TRUE,print.hist=FALSE, verbatim = TRUE,trackJ=TRUE)
ggimage(resVEM$Pg)
#---------------
r=2
H=(p+1):(p+r)
MH=VEM_fatala_200$M[,H]
data=MH %>% as_tibble() %>% mutate(Site=baran95$plan$site, date=baran95$plan$date,
                              Season=unlist(purrr::map(date, function(x){
                                if(x%in%c("oct93","jun93","aug93")) res="Rainy"
                                if(x%in%c("feb94","dec93","apr93")) res="Dry"
                                return(res)
                              }))) %>% 
  mutate(covar=paste(Site, date), covarnum=as.numeric(as.factor(covar))) 
 
 
 
g1<-data %>% 
  ggplot(aes(V1, V2, color=Site, shape=Site))+geom_point()+ 
  labs(x="Mh1",y="Mh2", title="")+scale_shape_discrete("Site",solid = TRUE)+
  scale_color_brewer("Site",palette="Dark2")+theme_light()

g2<-data%>% 
  ggplot(aes(V1, V2, color=Season, shape=Season))+geom_point()+theme_light()+
  scale_color_brewer("Season",palette="Dark2")+scale_shape_discrete("Season")+
  labs(x="Mh1",y=" Mh2", title="")
# g3<-data %>%
#    ggplot(aes( V1, Site, color=Site, fill=Site))+
#   geom_density_ridges(stat="binline",bins=20,draw_baseline=FALSE, alpha=0.8)+mytheme.dark("")+
#    labs(x="Mh1" ,y="", title="")+guides(color=FALSE, fill=FALSE)
# g4<-data %>%
#    ggplot(aes( V2, Season, color=Season, fill=Season))+
#   geom_density_ridges(stat="binline",bins=20,draw_baseline=FALSE, alpha=0.8)+mytheme.dark("")+
#    labs(x="Mh2",y="", title="")+guides(color=FALSE, fill=FALSE)
 g3<-data %>%
  ggplot(aes( V1, Site, color=Site, fill=Site))+
  geom_density_ridges( alpha=0.5)+mytheme.dark("")+
  labs(x="Mh1" ,y="", title="")+guides(color=FALSE, fill=FALSE)
g4<-data %>%
  ggplot(aes( V2, Season, color=Season, fill=Season))+
  geom_density_ridges( alpha=0.5)+mytheme.dark("")+
  labs(x="Mh2",y="", title="")+guides(color=FALSE, fill=FALSE)
g=grid.arrange(g1, g3,g2, g4, ncol=2)
ggsave(filename = "Fatala_MH2.png", plot = g,
       path ="/Users/raphaellemomal/these/R/images/", width = 7, height = 5)
