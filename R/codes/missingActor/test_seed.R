
#-------
# functions 


get.ListVEM<-function(seed){
  set.seed(seed) ; p=14 ; r=1 ;n=200 # 2 faible influence
  type="scale-free" ; O=1:p ;H=(p+1):(p+r); plot=FALSE 
  # Data
  missing_data<-missing_from_scratch(n,p,r,type,plot)
  counts=missing_data$Y; ZH=missing_data$ZH ; sigmaO= missing_data$Sigma; 
  omega=missing_data$Omega; trueClique=missing_data$TC[[1]]; hidden=missing_data$H
  # Observed parameters
  PLNfit<-PLN(counts~1, control=list(trace=0))
  MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S  ; theta=PLNfit$model_par$Theta ; 
  matcovar=matrix(1, n,1) ; sigma_obs=PLNfit$model_par$Sigma
  #------
  cliques_spca<- boot_FitSparsePCA(scale(MO),B=100,r=1, cores=3)
  ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, sigma_obs, MO,SO,r=1,alpha=0.1,
                    eps=1e-3,maxIter=100, nobeta=FALSE, cores=3,
                    filterDiag = FALSE, filterWg=FALSE,save=FALSE)
}

getDataVEM<-function(ListVEM,seed,n=200){
  set.seed(seed)
  missing_data<-missing_from_scratch(n=n,p=14,r=1,type="scale-free",plot=FALSE)
  omega=missing_data$Omega; 
  hidden=missing_data$H
  sorted_omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden),hidden)]
  diag(sorted_omega)=0
  
  badinit=which(do.call(rbind, lapply(ListVEM, function(x) length(x)))!=15)
  if(length(badinit)!=0) ListVEM=ListVEM[-badinit]
  
  J=do.call(rbind, lapply(ListVEM, function(vem){
    tail(vem$lowbound$J,1)
  }))
  Jcor=do.call(rbind, lapply(ListVEM, function(vem){
    getJcor(vem,14)
  }))
  AUC=do.call(rbind, lapply(ListVEM, function(vem){
    
    Pg=vem$Pg
    AUC=round(auc(pred = Pg, label = sorted_omega),4)
    
    return(AUC)}))
  ppvh=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg
    ppvh=accppvtpr(Pg,sorted_omega,h=15,seuil=0.5)[5] }))  
  penT=do.call(rbind, lapply(ListVEM, function(vem){
    penT=-( 0.5*sum( vem$Pg * log(vem$Wg+(vem$Wg==0)) ) - logSumTree(vem$Wg)$det) 
  }))
  
  sumP=do.call(rbind, lapply(ListVEM, function(vem){
    sum(vem$Pg)
  }))
  projL=do.call(rbind, lapply(ListVEM, function(vem){
    projL=sum(vem$projL)!=0
  }))
  
  LC=do.call(rbind, lapply(ListVEM, function(vem){
    LC= length(vem$clique[[1]])
  }))
  tprh=do.call(rbind, lapply(ListVEM, function(vem){
    Pg=vem$Pg  
    tprh=accppvtpr(Pg,sorted_omega,h=15,seuil=0.5)[8]
  }))
  
  data=data.frame(J, Jcor, AUC, ppvh,tprh,penT, sumP,projL,LC)
  data=data %>% mutate(ICL = Jcor-penT)
  return(data)
}

getPlots<-function(data, seed){
  data=data %>%as_tibble() %>%  mutate( maxJcor = ifelse(is.na(Jcor),FALSE,Jcor==max(Jcor, na.rm=TRUE)),
                                        maxJ=J==max(J),ICL=Jcor-penT)
  
  data %>%dplyr::select(J,Jcor,tprh,AUC,sumP,maxJcor) %>%
    gather(borne, vborne, -tprh,-AUC,-sumP,-maxJcor) %>% 
    gather(mesure, vmesure, -borne,-vborne,-sumP,-maxJcor) %>% 
    ggplot(aes(vmesure,vborne,color=as.factor(round(sumP,0)), size=(maxJcor)))+
    geom_vline(xintercept = 0.5, linetype="dashed",color="orange")+
    geom_point()+facet_grid(borne~mesure,scales = "free")+
    labs(title=paste0("seed ",seed))+
    mytheme.dark("sumP")
  
}

#------------
# processing

ListVEM_1  <-get.ListVEM(1)
ListVEM_7  <-get.ListVEM(7)
ListVEM_19 <-get.ListVEM(19)
ListVEM_2  <-get.ListVEM(2)
ListVEM_3  <-get.ListVEM(3)
ListVEM_6  <-get.ListVEM(6)
ListVEM_10 <-get.ListVEM(10)
ListVEM_111<-get.ListVEM(111)


data1   = getDataVEM(ListVEM_1  ,1  )
data7   = getDataVEM(ListVEM_7  ,7  )
data19  = getDataVEM(ListVEM_19 ,19 )
data2   = getDataVEM(ListVEM_2  ,2  )
data3   = getDataVEM(ListVEM_3  ,3  )
data6   = getDataVEM(ListVEM_6  ,6  )
data10  = getDataVEM(ListVEM_10 ,10 )
data111 = getDataVEM(ListVEM_111,111)

 
getPlots(data1  ,1  )
getPlots(data7  ,7  )
getPlots(data19 ,19 )
getPlots(data2  ,2  )
getPlots(data3  ,3  )
getPlots(data6  ,6  )
getPlots(data10 ,10 )
getPlots(data111,111)
#------
# tests
set.seed(seed)
missing_data<-missing_from_scratch(n=n,p=14,r=1,type="scale-free",plot=FALSE)
omega=missing_data$Omega; 
hidden=missing_data$H
sorted_omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden),hidden)]
diag(sorted_omega)=0
ggimage(ListVEM_7[[20]]$Pg)
