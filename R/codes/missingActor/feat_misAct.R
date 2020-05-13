# links of missing actor
build_misAct<-function(simus, H){
  linksH<-data.frame(do.call(rbind,lapply(simus, function(seed){seed$omega[H,-H]})))
  cpH<-data.frame(do.call(rbind,lapply(simus, function(seed){vecpartcor=-cov2cor(seed$omega)[H,-H]})))
  colnames(cpH)<-1:ncol(cpH)
  voisvois<-do.call(rbind,lapply(simus, function(seed){
    vois=1*(seed$omega[H,-H]!=0) # vect size 14 1 vois 0 pas vois
    di=round(diag(seed$omega)[-H],0) # degré des 14 premiers noeuds
    di[vois==0]=0 # mise à zéro des degré des non voisins de H
    return(di)
  }))
  voisMaxDeg<-do.call(rbind,lapply(simus, function(seed){
    index.vois=which(seed$omega[H,-H]!=0)
    vecmax<-rep(0,H-1)
    sapply(index.vois, function(vois){
      if(round(diag(seed$omega)[vois],0)>1){
        vv=setdiff(which(seed$omega[vois,-H]!=0),vois) # vecteur des noeuds voisins de vois à part H
        maxDeg=round(max(diag(seed$omega)[vv]),0)
        vecmax[vois]<<-maxDeg
      }
    })
    return(vecmax)
  }))
  voisBetween<-do.call(rbind,lapply(simus, function(seed){
    between=draw_network(seed$omega)$graph_data %>% activate(nodes) %>% as_tibble() %>% dplyr::select(btw)
    between=between$btw#/max(between$btw)
    between[seed$omega[H,-H]==0] = 0
    between=between[-H]
    return(between)
  }))
  Hbtw<-do.call(rbind,lapply(simus, function(seed){
    between=draw_network(seed$omega)$graph_data %>% activate(nodes) %>% as_tibble() %>% dplyr::select(btw)
    betweenH=between$btw[H]
    return(betweenH)
  }))
  
  misAct<-data.frame(nH=rowSums(1*(cpH!=0)), 
                     mean.cpH=apply(cpH, 1,function(vec){mean(vec[vec!=0])}),
                     mean.vv=apply(voisvois, 1,function(vec){mean(vec[vec!=0])}),
                     max.maxDeg=apply(voisMaxDeg, 1,function(vec){max(vec[vec!=0])}), 
                     mean.maxDeg=apply(voisMaxDeg, 1,function(vec){mean(vec[vec!=0])}), 
                     mean.btw=apply(voisBetween, 1,function(vec){mean(vec[vec!=0])}),
                     btwH=Hbtw) %>% as_tibble()
  return(misAct)
}
