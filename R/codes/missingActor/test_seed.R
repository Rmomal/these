
#-------
# functions 


get.ListVEM<-function(seed,r=1){
  set.seed(seed) ; p=14 ; r=r ;n=200 # 2 faible influence
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
  cliques_spca<- boot_FitSparsePCA(scale(MO),B=100,r=r, cores=3)
  ListVEM<-List.VEM(cliquesObj =cliques_spca, counts, sigma_obs, MO,SO,r=r,alpha=0.1,
                    eps=1e-3,maxIter=100, nobeta=FALSE, cores=3,updateSH = TRUE,
                    filterDiag = FALSE, filterWg=FALSE,save=FALSE)
}

getDataVEM<-function(ListVEM,seed,n=200,r=1,p=14){
  set.seed(seed)
  missing_data<-missing_from_scratch(n=n,p=p,r=r,type="scale-free",plot=FALSE)
  omega=missing_data$Omega; 
  hidden=missing_data$H
  sorted_omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden),hidden)]
  diag(sorted_omega)=0
  H=(p+1):(p+r)
 
  badinit=which(do.call(rbind, lapply(ListVEM, function(x) length(x)))!=15)
  if(length(badinit)!=0) ListVEM=ListVEM[-badinit]
  Jcor=do.call(rbind, lapply(ListVEM, function(vem){
    getJcor(vem,p)
  }))
  qual=do.call(rbind, lapply(ListVEM, function(vem){
    
    Pg=vem$Pg
    J= tail(vem$lowbound$J,1)
    omH=vem$omega[H,H]
    auc=round(auc(pred = Pg, label = sorted_omega),4)
    ppvh=accppvtpr(Pg,sorted_omega,h=H,seuil=0.5)[5]  
    penT=-( 0.5*sum(Pg * log(vem$Wg+(vem$Wg==0)) ) - logSumTree(vem$Wg)$det) 
    penZH= 0.5*sum(apply(matrix(vem$S[,H],n,r),2,function(x){ log(sum(x))}))+r*n*0.5*(1+log(2*pi))
    sumP=abs(sum(vem$Pg)-2*((p+r)-1))
    projL=vem$projL 
    nbinit= vem$nbinit
    tprh=accppvtpr(Pg,sorted_omega,h=H,seuil=0.5)[8]
    return(c(J=J,omH=omH,auc=auc, ppvh=ppvh,tprh=tprh,sumP=sumP,projL=projL, nbinit=nbinit, penT=penT))
  }))
  data=data.frame(Jcor, qual)
  data=data %>% mutate(ICL = Jcor-penT)
  return(data)
}

getPlots<-function(data, seed){
  data=data %>% as_tibble() %>% 
    mutate(maxJcor=ifelse(sumP<1e-10,FALSE,Jcor==max(Jcor, na.rm=TRUE)))
  
  data %>%dplyr::select(J,Jcor,tprh,auc,sumP,maxJcor, omH) %>%
    gather(borne, vborne, -tprh,-auc,-sumP,-maxJcor, -omH) %>% 
    gather(mesure, vmesure, -borne,-vborne,-sumP,-maxJcor, -omH) %>% 
    ggplot(aes(vmesure,vborne,color=omH>25 ))+
    geom_vline(xintercept = 0.5, linetype="dashed",color="orange")+
    geom_point()+facet_grid(borne~mesure,scales = "free")+
    labs(title=paste0("seed ",seed))+#theme_light()
    mytheme.dark("omH>25")
  
}

#------------
# processing

ListVEM_1  <-get.ListVEM(1)
ListVEM_7  <-get.ListVEM(7)
ListVEM_19 <-get.ListVEM(19)
ListVEM_2  <-get.ListVEM(2)
ListVEM_3  <-get.ListVEM(3)
ListVEM_6  <-get.ListVEM(6)
ListVEM_111<-get.ListVEM(111)
ListVEM_10 <-get.ListVEM(10)
ListVEM_27 <-get.ListVEM(27)

data1   = getDataVEM(SF_seed1$ListVEM  ,1  )
data7   = getDataVEM(ListVEM_7  ,7  )
data192  = getDataVEM(SF_seed19$ListVEM ,19 )
data19  = getDataVEM(ListVEM_19,19 )
data2   = getDataVEM(ListVEM_2  ,2  )
data3   = getDataVEM(ListVEM_3  ,3  )
data6   = getDataVEM(ListVEM_6  ,6  )
data10  = getDataVEM(ListVEM_10 ,10 )
data111 = getDataVEM(ListVEM_111,111)
data27 = getDataVEM(ListVEM_27,27)


getPlots(data1  ,1  )
getPlots(data7  ,7  )
getPlots(data19 ,19 )
getPlots(data192 ,19 )
getPlots(data2  ,2  )
getPlots(data3  ,3  )
getPlots(data6  ,6  )
getPlots(data10 ,10 )
getPlots(data111,111)
getPlots(data27,27)
#------
# tests
set.seed(seed)
missing_data<-missing_from_scratch(n=n,p=14,r=1,type="scale-free",plot=FALSE)
omega=missing_data$Omega; 
hidden=missing_data$H
sorted_omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden),hidden)]
diag(sorted_omega)=0
ggimage(ListVEM_7[[20]]$Pg)
