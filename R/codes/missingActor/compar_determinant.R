#test formule lauritzen
# arbre à 3 noeuds
arbre=matrix(c(10,1,0,1,10,1,0,1,1),3, 3, byrow = TRUE)
Iphi=1*(arbre!=0)*(CorOmegaMatrix(arbre))
(det(arbre)-exp(sum(log(diag(arbre)))+0.5*sum(log(Iphi+(Iphi==0)))))/det(arbre)



######
B=30
diffs=lapply(1:B, function(i){
  set.seed(i)
  n=200 ;p=14;r=0;type="scale-free" 
  missing_data<-missing_from_scratch(n,p,r,type,plot=FALSE)
  arbre=missing_data$Omega
  sigma=missing_data$Sigma
  det(arbre)
  
  Iphi=1*(arbre!=0)*(CorOmegaMatrix(sigma))
  diffSig=(det(sigma)-exp(sum(log(diag(sigma)))+sum(log(Iphi+(Iphi==0))))) 
  Iphi=1*(arbre!=0)*(CorOmegaMatrix(arbre))
  diffOm=(det(arbre)-exp(sum(log(diag(arbre)))+sum(log(Iphi+(Iphi==0)))))
  res=data.frame(diffSig=diffSig, diffOm=diffOm)
  return(res)
})
datadiff=do.call(rbind, diffs) %>% as_tibble()
summary(datadiff)
Iphi=1*(arbre!=0)*(CorOmegaMatrix(sigma))
det(Iphi)
det(arbre%*%CorOmegaMatrix(sigma))

# comparer corrélations et corrélations partielles
set.seed(3)
n=1e+6 ;p=14;r=0;type="scale-free" 
missing_data<-missing_from_scratch(n,p,r,type,plot=FALSE)
arbre=missing_data$Omega
diag(arbre)=2
sigma=solve(arbre)
corr=F_Sym2Vec(1*(arbre!=0)*(1-cov2cor(sigma)^2))
corrp=F_Sym2Vec(1*(arbre!=0)*CorOmegaMatrix(solve(sigma)))
plot(corr, corrp)

#comparer matrice de pécision et minverse de la matrice de corrélation
# vérifier la formule pour sigmaTilde
set.seed(1)
n=200 ;p=14;r=0
missing_data<-missing_from_scratch(n,p,r,type,plot=TRUE)
counts=missing_data$Y
PLNfit=PLN(counts~1)
Sigma=PLNfit$model_par$Sigma
Cor = cov2cor(Sigma)
D=sqrt(diag(1/diag(Sigma)))
Cor_t=D%*%Sigma%*%D
plot(Cor, Cor_t)

omega=solve(Sigma)
corom=cov2cor(solve(Sigma))
corom_t=solve(cov2cor(Sigma))
plot(corom, corom_t)
abline(0,1)


plot(diag(solve(Cor)),diag(solve(Sigma)))
abline(0,1)
# en sortie de PLN, on normalise les Z

# corrélation en fonction des corrélations partielles obtenues
# par l'inverse de la matrice de corrélation
arbre=1*(missing_data$Omega!=0)
Sigma=missing_data$Sigma
cor_t=cov2cor(Sigma)
gamma= solve(cor_t)
plot((arbre*cor_t)^2, -(arbre*cov2cor(gamma)))
abline(0,1)

# corrélations en fonction des corrélations partielles obtenues par l'inverse
# de la matrice de variance
set.seed(1)
n=200 ;p=14;r=0
missing_data<-missing_from_scratch(n,p,r,type,plot=TRUE)
counts=missing_data$Y
PLNfit=PLN(counts~1)
Sigma=PLNfit$model_par$Sigma
omega=missing_data$Omega
arbre=1*(missing_data$Omega!=0)
plot(F_Sym2Vec((1-arbre*cov2cor(Sigma)^2)^2), F_Sym2Vec(arbre*cov2cor(omega))^2)
abline(0,1)
omega[1,2]^2/(omega[1,1]*omega[2,2])
cov2cor(omega)[1,2]^2

######
# verif formule 5.12
omega=missing_data$Omega
arbre=1*(omega!=0)
corp=cov2cor(omega)
sigma=solve(omega)
detS=det(sigma)
res=matrix(NA, 14, 14)
for(mu in 1:14){
  for(gam in mu:14){
    detS_mu = det(sigma[-mu,-mu])
    detS_gam=det(sigma[-gam,-gam])
    detS_mugam=det(sigma[-c(mu,gam),-c(mu,gam)])
   res[mu,gam]=arbre[mu,gam]*(corp[mu,gam]^2-(1-((detS*detS_mugam)/(detS_mu*detS_gam))))
  }
}
ggimage(res)
diag(arbre)=0

det.fractional(sigma, log=TRUE)
res=matrix(NA, 14, 14)
for(mu in 1:14){
  for(gam in mu:14){
    detS_mu = 0.5*sum(log((1-cov2cor(sigma[-mu,-mu])[arbre[-mu,-mu]!=0]^2)))+
      sum(log(diag(sigma[-mu,-mu])))
    detS_gam=0.5*sum(log((1-cov2cor(sigma[-gam,-gam])[arbre[-gam,-gam]!=0]^2)))+
      sum(log(diag(sigma[-gam,-gam])))
    detS_mugam=0.5*sum(log((1-cov2cor(sigma[-c(mu,gam),-c(mu,gam)])[arbre[-c(mu,gam),-c(mu,gam)]!=0]^2)))+
      sum(log(diag(sigma[-c(mu,gam),-c(mu,gam)])))
    res[mu,gam]=arbre[mu,gam]*(corp[mu,gam]^2-(1-((detS*detS_mugam)/(detS_mu*detS_gam))))
  }
}
ggimage(res)

plot(CorOmegaMatrix(sigma), 1-cov2cor(sigma)^2)
prod(1-cov2cor(sigma)^2*arbre)
prod(arbre*CorOmegaMatrix(sigma) + (arbre==0))

mu=1
detS_1=log(det(sigma[-1,-1]))
detS_mu = 0.5*sum(log((1-cov2cor(sigma[-mu,-mu])[arbre[-mu,-mu]!=0]^2)))+
  sum(log(diag(sigma[-mu,-mu])))

log(det(sigma)) - 0.5*sum(log((1-cov2cor(sigma)[arbre!=0]^2)))-sum(log(diag(sigma)))
