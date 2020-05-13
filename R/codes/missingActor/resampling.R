# codage de l'estimateur GAP, resampling pour choix du nombre de clusters

seed=7 ; p=14 ; n=200
r=1
# données d'origine
set.seed(seed)
# sim data
missing_data<-missing_from_scratch(n,p,r=r,type,plot=FALSE)
counts=missing_data$Y; omega=missing_data$Omega
# Observed parameters
PLNfit<-PLN(counts~1, control=list(trace=0))
MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S ; 
theta=PLNfit$model_par$Theta; sigma_obs=PLNfit$model_par$Sigma
init0=initVEM(counts , initviasigma = NULL,  sigma_obs,r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 

init=initVEM(counts = counts,initviasigma=NULL, sigma_obs,r = 0) #cliques_spca$cliqueList[[4]]
Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit

VEM0=VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=1e-3, 
                alpha=0.1, 
                maxIter=50, plot=TRUE,print.hist=FALSE,filterWg=TRUE,updateSH = TRUE,
                verbatim = TRUE,nobeta = FALSE, filterDiag = FALSE)
Sigma_tree=(as.matrix(solve(nearPD(VEM0$omega*VEM0$Pg, eig.tol = 1e-8, posd.tol = 1e-7,doSym = TRUE)$mat)))
Sigma_tree=Sigma_tree/mean(Sigma_tree) 
det(Sigma_tree)
ggimage(Sigma_tree)
# modèle nul = données sans structure
resamp_nullJcor<-function(sigma_obs,r,alpha=0.1,eps=1e-3){
  p=ncol(sigma_obs)
  null_counts=generator_PLN(Sigma_tree,covariates = NULL, n=200)$Y
  PLN_null<-PLN(null_counts~1, control=list(trace=0))
  MO<-PLN_null$var_par$M  ; SO<-PLN_null$var_par$S  ; sigma_obs=PLN_null$model_par$Sigma
  if(r!=0){
    clique=boot_FitSparsePCA(null_counts,B=10,r=r,cores=3)# FitSparsePCA(null_counts,r=r)$cliques
    
    cat(paste0("computing ", length(clique$cliqueList)," VEM...\n"))
    ListVEM<-List.VEM(cliquesObj =clique, null_counts, sigma_obs, MO,SO,r=r,alpha=0.1,
                      eps=eps,maxIter=100, nobeta=FALSE, cores=3,
                      filterDiag = FALSE,updateSH=TRUE, filterWg=TRUE,save=FALSE)
    cat(paste0("gathering Jcordata...\n"))
    Jcordata=do.call(rbind, lapply(ListVEM, function(vem){
      if(length(vem)==15){
        Jcor=getJcor(vem,p, 1e-8, 1e-7)
      }else{
        Jcor=c(Jcor=NaN,diff=NaN,detEg=NaN)
      } 
    })) %>% as_tibble()
    Jcordata=Jcordata #%>% filter(detEg>-40) %>% filter(Jcor==max(Jcor, na.rm=TRUE))
  }else{
    clique=NULL
  }
  
  return(Jcordata)
}
tic()
echant_J_nullr1<-lapply(1:3,function(x){
  cat(paste0("\n resampling ",x,"\n")) # 7min r=1
  resamp_nullJcor(sigma_obs,r=1)
})
toc()

echant_J_nullr1=do.call(rbind, echant_J_nullr1) 
echant_J_nullr1 %>% as_tibble() %>% filter(!is.na(Jcor)) %>% 
  ggplot(aes(Jcor))+geom_density()+theme_light()
