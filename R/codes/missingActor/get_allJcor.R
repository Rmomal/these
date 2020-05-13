
AllJcor = data.frame(Jcor=numeric(), diff=numeric(), detEg=numeric(), sumP=numeric(), auc=numeric(), ppvh=numeric(), tprh=numeric(), seed=numeric())
simus=lapply(1:N, function(seed){
  file= paste0("/Users/raphaellemomal/simulations/15nodes_rawdata/SF_seed",seed,".rds")
  simu=readRDS(file)
  ome=simu$omega
  diag(ome)=0
  Jcor_data=data.frame(do.call(rbind,lapply(simu$ListVEM, function(vem){
    if(length(vem)==15){
      getJcor(vem,p=14,1e-8,1e-7)
    }
  })),do.call(rbind,lapply(simu$ListVEM, function(vem){
    if(length(vem)==15){
      Pg=vem$Pg
      sumP=abs(sum(Pg)-28)
      auc=round(auc(pred = Pg, label = ome),4)
      ppvh=accppvtpr(Pg,ome,h=H,seuil=0.5)[5]
      tprh=accppvtpr(Pg,ome,h=H,seuil=0.5)[8]
      return(c(sumP=sumP,auc=auc, ppvh=ppvh, tprh=tprh))
    }
  }))) %>% as_tibble() %>% mutate(seed=seed)
  AllJcor<<-rbind(AllJcor,Jcor_data)
  simu$ListVEM = simu$ListVEM[which.max(Jcor_data$Jcor)]
  return(simu)
})