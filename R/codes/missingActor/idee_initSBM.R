
######
# r=1
# data
#run EM
# SBM cliques

######
# r=0
r=0
p=14
set.seed(7)
missing_data0<-missing_from_scratch(n,p,r,type,plot=TRUE)
counts0=missing_data0$Y;  omega=missing_data0$Omega; 
hub=missing_data0$H ; clique=missing_data0$TC
# regarder si spce veut les voisins de 2 comme clique manquante. Si oui peut-être la différence
# entre corrélations et corrélations partielles et mettre au point une méthode
# d'initialisation à partir de la connectivité intra des communautés construites avec les modèles graphiques

FitSparsePCA(scale(MO), r=1)$cliques
# Observed parameters
PLNfit0<-PLN(counts0~1, control=list(trace=0))
MO<-PLNfit0$var_par$M  ; SO<-PLNfit0$var_par$S  ; theta=PLNfit0$model_par$Theta ; 
matcovar=matrix(1, n,1) ; sigma_obs=PLNfit0$model_par$Sigma

ggimage(cov2cor(solve(sigma_obs)))
if(r!=0){
  sorted_omega=omega[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
}else{
  sorted_omega=omega
}
#VEM with no missing actor
init0=initVEM(counts = counts, initviasigma = NULL,  sigma_obs,r = 0)
Wginit= init0$Wginit; Winit= init0$Winit; omegainit=init0$omegainit 
alpha<-tryCatch(expr={computeAlpha(omegainit,default =0.3, MO, SO, plot=plot)},
                error=function(e){message("sythetic alpha")
                  return(0.3)})

# data
#run EM
# SBM cliques
#r=1
r=1
missing_data<-missing_from_scratch(n,p,r,type,plot=TRUE)
counts=missing_data$Y;  omega=missing_data$Omega; 
sigma0=missing_data$Sigma ; 
h=missing_data$H ; trueClique=missing_data$TC
sorted_omega=omega[c(setdiff(1:(p+r), h), h),c(setdiff(1:(p+r), h), h)]
PLNfit<-PLN(counts~1, control=list(trace=0))
probhat=EMtree(PLNfit)$edges_prob
G=draw_network(probhat, 
               layout="kk",curv=0,nb=2,pal="black",nodes_label =1:(p))$G
print(G)
