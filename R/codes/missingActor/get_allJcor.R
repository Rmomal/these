
AllJcor = data.frame(Jcor=numeric(), diff=numeric(), detEg=numeric(), sumP=numeric(), auc=numeric(), ppvh=numeric(), tprh=numeric(), seed=numeric())
N=400
simus=lapply(1:N, function(seed){
  file= paste0("/Users/raphaellemomal/simulations/15nodes_V2/SF_seed",seed,".rds")
  simu=readRDS(file)
  #G=simu$G
  J=do.call(rbind,lapply(simu$ListVEM, function(vem){
    if(length(vem)==14){
      tail(vem$lowbound$J,1)
  }else{ NaN}}))
  # Jcor_data=data.frame(do.call(rbind,lapply(simu$ListVEM, function(vem){
  #   if(length(vem)==14){
  #     res=getJcor(vem,p=14,1e-14,1e-14)
  #   }else{
  #     res=c(NaN,NaN,NaN,NaN)
  #   }
  #   return(res)
  # })),do.call(rbind,lapply(simu$ListVEM, function(vem){
  #   if(length(vem)==14){
  #     Pg=vem$Pg
  #     sumP=abs(sum(Pg)-28)
  #     auc=round(auc(pred = Pg, label = ome),4)
  #     ppvh=accppvtpr(Pg,ome,h=H,seuil=0.5)[5]
  #     tprh=accppvtpr(Pg,ome,h=H,seuil=0.5)[8]
  #     res=c(sumP=sumP,auc=auc, ppvh=ppvh, tprh=tprh)
  #   }else{
  #     res=c(NaN,NaN,NaN,NaN)
  #   }
  #   return(res)
  # }))) %>% as_tibble() %>% mutate(seed=seed)
  # #AllJcor<<-rbind(AllJcor,Jcor_data)
  # Jcor_data$num = 1:length(simu$ListVEM)
  # #filtre 1 : 0.1 et -35
  # #filtre 2 : 1e-5 et -45
  # num=unlist(Jcor_data %>% filter(sumP<1e-10 ) %>% filter(Jcor==max(Jcor)) %>% 
  #   dplyr::select(num))
  # if(length(num)==0) num=unlist(Jcor_data %>% filter(sumP==min(sumP) )  %>%   dplyr::select(num))
  num=which.max(J)
  simu$ListVEM = simu$ListVEM[num]
  return(simu)
})
for(i in 1:N){
  names(simus[[i]])<-c("G" ,"UH"  , "VEM_1","time_boots", "nbinit"  )
  simus[[i]]$VEM_1=unlist(simus[[i]]$VEM_1, recursive = FALSE)
  }
 
saveRDS(AllJcor, "/Users/raphaellemomal/these/R/codes/missingActor/SimResults/AllJcor_14.rds")
saveRDS(simus, "/Users/raphaellemomal/these/R/codes/missingActor/SimResults/simus_V2_400.rds")
simus=readRDS("/Users/raphaellemomal/these/R/codes/missingActor/SimResults/simus_V2_216.rds")
which(do.call(rbind, lapply(simus, function(seed){
  length(seed)
}))!=5)
AllJcor %>%filter(detEg>-35) %>%  ggplot(aes(y=log(sumP+1e-16), color=(tprh>0.8)))+geom_boxplot()+theme_light()
AllJcor %>% ggplot(aes(y=detEg, color=(tprh>0.8)))+geom_boxplot()+theme_light()

AllJcor %>% group_by(tprh>0.7) %>% summarise(meanP=mean(sumP), maxP=max(sumP),
                                             meandet=mean(detEg) ,minDet=min(detEg) )             

jselec=AllJcor %>% filter(auc>0.5 , tprh>0.8) %>% dplyr::select(Jcor)                                      
which.max(jselec)

AllJcor %>% filter(sumP<1e-8, detEg>-65) %>%  group_by(seed) %>%
  filter(Jcor==max(Jcor))  
setdiff(1:400, seeds)
AllJcor %>% filter(seed%in%c(74, 148))
AllJcor %>% filter(sumP<1e-10,auc<0.5)
 
