# codage de l'estimateur GAP, resampling pour choix du nombre de clusters

seed=7 ; p=14 ; n=200
r=1
set.seed(seed)
# sim data
missing_data<-missing_from_scratch(n,p,r=r,type,plot)
counts=missing_data$Y 
# Observed parameters
PLNfit<-PLN(counts~1, control=list(trace=0))
MO<-PLNfit$var_par$M  ; SO<-PLNfit$var_par$S ; 
theta=PLNfit$model_par$Theta; sigma_obs=PLNfit$model_par$Sigma
#-- normalize the PLN outputs
D=diag(sigma_obs)
matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
MO=MO*matsig
SO=SO*matsig^2
#initialize
init0=initVEM(counts , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
# vem with original counts
VEM0<-VEMtree(counts, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=1e-3,
                Wg_init =Wginit,plot = FALSE, maxIter = 100,print.hist = FALSE,
                alpha=0.1, verbatim=FALSE, trackJ=FALSE )

Sigma_tree=cov2cor(as.matrix(solve(
  nearPD(VEM0$Upsilon*VEM0$Pg, eig.tol = 1e-14, posd.tol = 1e-14,doSym = TRUE)$mat)))
#Sigma_tree=Sigma_tree/mean(Sigma_tree) 
det(Sigma_tree)
ggimage(Sigma_tree)
# modèle nul = données sans structure
resamp_nullJcor<-function(sigma_obs,r,B=10,alpha=0.1,eps=1e-3){
  p=ncol(sigma_obs)
  null_counts=generator_PLN(Sigma_tree,covariates = NULL, n=200)$Y
  PLN_null<-PLN(null_counts~1, control=list(trace=0))
  MO<-PLN_null$var_par$M  ; SO<-PLN_null$var_par$S  ; sigma_obs=PLN_null$model_par$Sigma
  #-- normalize the PLN outputs
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  if(r!=0){
    clique=boot_FitSparsePCA(null_counts,B=B,r=r,cores=3) 
    cat(paste0("computing ", length(clique$cliqueList)," VEM...\n"))
    ListVEM<-List.VEM(cliquesObj =clique, null_counts, cov2cor(sigma_obs), MO,SO,r=r,alpha=0.1,
                      eps=eps,maxIter=100, cores=3, trackJ=FALSE)
    cat(paste0("gathering Jcordata...\n"))
    J=max(do.call(rbind, lapply(ListVEM, function(vem){tail(vem$loxbound$J, 1)})))
  }else{
    #initialize
    init0=initVEM(counts , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
    Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
    VEM<-VEMtree(counts, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=1e-3,
                  Wg_init =Wginit,plot = FALSE, maxIter = 100,print.hist = FALSE,
                  alpha=0.1, verbatim=TRUE, trackJ=FALSE )
    J=tail(VEM$lowbound$J,1)
  }
  return(J)
}
tic()
echant_J_nullr1<-lapply(1:3,function(x){
  cat(paste0("\n resampling ",x,"\n")) # 7min r=1
  resamp_nullJcor(sigma_obs,r=0)
})
toc()

echant_J_nullr1=do.call(rbind, echant_J_nullr1) 
echant_J_nullr1 %>% as_tibble() %>% filter(!is.na(Jcor)) %>% 
  ggplot(aes(Jcor))+geom_density()+theme_light()
