setwd("/home/momal/Git/these/Pgm_SR")
source('FunctionsMatVec.R')
source('FunctionsTree.R')
source('FunctionsInference.R')
source('fonctions.R')
library(mvtnorm)
library(PLNmodels)

# Sim parms
p = 20; n =100; d = .1
par(pch=20)
B = 1e3

# Donne l'estimateur du nombre d'abdence d'arêtes
EstimM <- function(Prob, G){
   p = ncol(Prob); 
   M = 2*sum(Prob[upper.tri(Prob)]>.5); 
   # hist(as.vector(Prob), breaks=p, main=paste(M, '/', sum(G[upper.tri(G)]==0)))
   return(M)
}

# Sim loop
Res = matrix(0, B, 7)
for (b in 1:B){
   if (b %% round(sqrt(B))==0){cat(b, '')}
   # Graph
   #G = erdos(p, d)  
  G = SimCluster(p, 3, d, 10)
   
   # Model parms
   lambda = 1; Omega = diag(rep(lambda, p)) + G; 
   while (min(eigen(Omega)$values) < 0){lambda = 1.1*lambda; Omega = diag(rep(lambda, p)) + G}
   Sigma = solve(Omega)
   O = matrix(rnorm(n*p), n, p)
   
   # Data
   Z = rmvnorm(n, sigma=Sigma); 
   Y = matrix(rpois(n*p, exp(O+Z)), n, p)
   PLN = PLN(Y ~ -1 + offset(O))

   # Inference based on Z : regression
   ZpZ = t(Z)%*%Z; invZpZ = solve(ZpZ)
   StatZ1 = PvalZ1 = matrix(0, p, p)
   for (j in 1:p){
      # LM = lm(Z[, j] ~ -1 + Z[, -j])
      # StatZ1[j, -j] = summary(LM)$coef[, 3]
      # PvalZ1[j, -j] = summary(LM)$coef[, 4]
      invXpX = solve(ZpZ[-j, -j])
      Beta = invXpX%*%ZpZ[-j, j]
      SdBeta = sqrt(diag(invXpX)) / sqrt(invZpZ[j, j]) / sqrt(n-p+1)
      StatZ1[j, -j] = Beta/SdBeta
      PvalZ1[j, -j] = 2*pt(abs(StatZ1[j, -j]), df=(n-p+1), lower.tail=F)
   }
   
   # Inference based on Y : regression
   ZpZ.hat = PLN$model_par$Sigma; invZpZ.hat = solve(ZpZ.hat)
   StatY1 = PvalY1 = matrix(0, p, p)
   for (j in 1:p){
      X = log(1+Y[, -j]/exp(O[, -j]))
      # GLM = glm(Y[, j] ~ X + offset(O[, j]), family=poisson)
      # StatY1[j, -j] = summary(GLM)$coef[-1, 3]
      # PvalY1[j, -j] = summary(GLM)$coef[-1, 4]
      invXpX = solve(ZpZ.hat[-j, -j])
      Beta = invXpX%*%ZpZ.hat[-j, j]
      SdBeta = sqrt(diag(invXpX)) / sqrt(invZpZ.hat[j, j]) / sqrt(n-p+1)
      StatY1[j, -j] = Beta/SdBeta
      PvalY1[j, -j] = 2*pt(abs(StatY1[j, -j]), df=(n-p+1), lower.tail=F)
   }
   
   # par(mfrow=c(2, 2))
   # plot(as.vector(cor(Z)), as.vector(cov2cor(PLN$model_par$Sigma))); abline(0, 1)
   # plot(as.vector(cov(Z)), as.vector(PLN$model_par$Sigma)); abline(0, 1)
   # plot(as.vector(StatZ1), as.vector(StatY1)); abline(0, 1)
   # plot(as.vector(PvalZ1), as.vector(PvalY1), log='xy'); abline(0, 1)
   
   # Inference based on Z : correlation
   RZ = cor(Z)
   StatZ2 = RZ * sqrt((n-2)/(1-RZ^2))
   PvalZ2 = matrix(2*pt(abs(StatZ2), lower.tail=F, df=n-2), p, p)
   
   # Inference based on Z : partial correlation
   OmegaZ = solve(cov(Z))
   RpartZ = -diag(1/sqrt(diag(OmegaZ)))%*%OmegaZ%*%diag(1/sqrt(diag(OmegaZ)))
   StatZ3 = RpartZ * sqrt((n-2)/(1-RpartZ^2))
   PvalZ3 =  matrix(2*pt(abs(StatZ3), lower.tail=F, df=n-p-2), p, p)
   
   # Inference based on Y : correlation
   RY = cov2cor(PLN$model_par$Sigma)
   StatY2 = RY * sqrt((n-2)/(1-RY^2))
   PvalY2 = matrix(2*pt(abs(StatY2), lower.tail=F, df=n-2), p, p)
   #PvalY2 = matrix(2*pnorm(abs(StatY2), lower.tail=F), p, p)
   
   # Inference based on Y : partial correlation
   OmegaY = solve(PLN$model_par$Sigma)
   RpartY = -diag(1/sqrt(diag(OmegaY)))%*%OmegaY%*%diag(1/sqrt(diag(OmegaY)))
   StatY3 = RpartY * sqrt((n-2)/(1-RpartY^2))
   PvalY3 =  matrix(2*pt(abs(StatY3), lower.tail=F, df=n-p-2), p, p)
   #PvalY3 =  matrix(2*pnorm(abs(StatY3), lower.tail=F), p, p)
   
   # # Comparison
   # Cor = as.data.frame(cbind(F_Sym2Vec(cor(Z)), F_Sym2Vec(cov2cor(PLN$model_par$Sigma)), 
   #                           F_Sym2Vec(cov2cor(OmegaZ)), F_Sym2Vec(cov2cor(OmegaY))))
   # names(Cor) = c('corZ', 'corY', 'partZ', 'partY')
   # plot(Cor)
   
   # Stat = as.data.frame(cbind(F_Sym2Vec(StatZ1), F_Sym2Vec(StatY1), F_Sym2Vec(StatZ2), F_Sym2Vec(StatY2), F_Sym2Vec(StatZ3), F_Sym2Vec(StatY3)))
   # Pval = as.data.frame(cbind(F_Sym2Vec(PvalZ1), F_Sym2Vec(PvalY1), F_Sym2Vec(PvalZ2), F_Sym2Vec(PvalY2), F_Sym2Vec(PvalZ3), F_Sym2Vec(PvalY3)))
   # names(Stat) = names(Pval) = c('Z1', 'Y1', 'Z2', 'Y2', 'Z3', 'Y3')
   # cor(Stat); plot(Stat); plot(Pval)
  
   MZ1 = EstimM(PvalZ1, G); MZ2 = EstimM(PvalZ2, G); MZ3 = EstimM(PvalZ3, G)
   MY1 = EstimM(PvalY1, G); MY2 = EstimM(PvalY2, G); MY3 = EstimM(PvalY3, G)
   Res[b, ] = c(sum(G[upper.tri(G)]==0), MZ1, MY1, MZ2, MY2, MZ3, MY3)
   colnames(Res)<-c("ref","regZ","regY","corZ","corY","partCorZ","partCorY")
}

Diff = Res[, 2:7] - Res[, 1]%o%rep(1, 6)
graphics.off()
#boxplot(Diff, col=rep(c(2, 4), 3)); abline(h=0)
dataDiff<-Diff %>%  data.frame() %>%  gather(col,val)
dataDiff$layer<-ifelse(grepl("Z", dataDiff$col),"Z","Y")
ggplot(dataDiff,aes(col,val,fill=layer))+ 
  geom_hline(yintercept=0,col="darkcyan",linetype="dashed",size=1.2)+
  geom_violin(position=position_dodge(1),alpha=0.5)+
  geom_boxplot(width=0.15,outlier.colour="black", outlier.shape=16,
             outlier.size=1.5,fill="white",alpha=0.5, notch=FALSE)+
 
  scale_fill_brewer(palette="Dark2") + 
  labs(x=" ",y="Difference")+
  theme_minimal()+
  theme(axis.text=element_text(size=11))
  
  #geom_dotplot(binaxis='y', stackdir='center',binwidth=.3,   position=position_dodge(1))+ 



