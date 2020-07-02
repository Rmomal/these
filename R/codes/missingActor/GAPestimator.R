# codage de l'estimateur GAP, resampling pour choix du nombre de clusters

seed=7 ; p=14 ; n=200
r=0
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
VEM0<-VEMtree(counts, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=1e-5,
              Wg_init =Wginit,plot = TRUE, maxIter = 100,print.hist = FALSE,
              alpha=0.1, verbatim=TRUE, trackJ=TRUE )
tail(VEM0$lowbound$J, 1)
Sigma_tree=as.matrix(solve( VEM0$Upsilon*(VEM0$Pg+diag(p))))
ggimage(VEM0$Pg)
ggimage(missing_data$G)
# modèle nul = données sans structure
resamp_nullJcor<-function(sigma_obs,r,B=10,alpha=0.1,eps=1e-3){
  p=ncol(sigma_obs)
  null_counts=generator_PLN(sigma_obs,covariates = NULL, n=200)$Y
  PLN_null<-PLN(null_counts~1, control=list(trace=0))
  MO<-PLN_null$var_par$M  ; SO<-PLN_null$var_par$S  ; sigma_obs=PLN_null$model_par$Sigma
  #-- normalize the PLN outputs
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  
  #initialize
  init0=initVEM(null_counts , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
  Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
  VEM<-VEMtree(null_counts, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=1e-3,
               Wg_init =Wginit,plot = FALSE, maxIter = 100,print.hist = FALSE,
               alpha=0.1, verbatim=TRUE, trackJ=FALSE )
  J=tail(VEM$lowbound$J,1)
  
  return(J)
}
# apprentissage null model r=0
set.seed(7)
tic()
echant_J_nullr0_0<-mclapply(1:100,function(x){
  cat(paste0("\n resampling ",x,"\n")) # 63s 3 cores r=0
  resamp_nullJcor(Sigma_tree,r=0)
}, mc.cores=3)
toc()
badinit=which(lapply(echant_J_nullr0_0, typeof)=="character")
if(length(badinit!=0)) echant_J_nullr0=echant_J_nullr0[-badinit]
echant_J_nullr0_0=data.frame(J=do.call(rbind, echant_J_nullr0_0) )
echant_J_nullr0_0 %>% as_tibble() %>%  
  ggplot(aes(J))+geom_density()+theme_light()
echant_J_nullr0 %>% as_tibble() %>%  
  ggplot(aes(J))+geom_boxplot()+theme_light()
summary(echant_J_nullr0)

# means
m0 = mean(echant_J_nullr0_0$J)

##########
# tests
# trueR=0
# test with r 0
init0=initVEM(counts , initviasigma = NULL,  cov2cor(sigma_obs),r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; upsinit=init0$upsinit 
# vem with original counts
VEM0<-VEMtree(counts, MO, SO, MH=NULL,upsinit,W_init =Winit,eps=1e-3,
              Wg_init =Wginit,plot = FALSE, maxIter = 100,print.hist = FALSE,
              alpha=0.1, verbatim=FALSE, trackJ=FALSE )
Jtest_0=tail(VEM0$lowbound$J,1)
echant_J_nullr0_0 %>% as_tibble() %>% ggplot(aes(J))+
  geom_density(color="cornflowerblue",fill="cornflowerblue", alpha=0.5)+theme_light()+
  geom_vline(xintercept = Jtest_0, color="black", linetype="dashed")
# gap estimator : comparer les distancecs et choisir le modèle qui a 
# la plus petite distance aux prédictions du modèle nul
plot(ecdf(echant_J_nullr0_0$J))
# plot
echant_J_nullr1=echant_J_nullr1[is.finite(unlist(echant_J_nullr1)),1]
data=rbind(cbind(J=echant_J_nullr0$J,r=0),
           cbind(J=echant_J_nullr1, r=1))
data %>%as_tibble() %>%
  ggplot(aes(J,fill=as.factor(r),color=as.factor(r), group=as.factor(r)))+
  geom_density(alpha=0.5)+theme_light()+guides(color=FALSE, fill=FALSE)+
  #geom_vline(xintercept = Jtest, color="black", linetype="dashed")+
  geom_vline(xintercept = Jtest_0, color="black", linetype="dashed")+
  #geom_vline(xintercept = m1, color="darkturquoise")+
  geom_vline(xintercept = m0, color="coral")

# diff in J
(Jtest-m0)/sd(echant_J_nullr0$J)^2
(Jtest-m1)/sd(echant_J_nullr1 )^2
(Jtest_0-m0)/sd(echant_J_nullr0$J)^2
(Jtest_0-m1)/sd(echant_J_nullr1 )^2
