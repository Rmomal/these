# try to understand numerical issues again.
# big alpha 
library(nestor)
data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
PLNfit<-norm_PLN(data$Y)
MO<-PLNfit$MO
SO<-PLNfit$SO
sigma_O=PLNfit$sigma_O
#-- initialize with true clique for example
initClique=data$TC
#-- initialize the VEM
initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=1 )
#-- run core function nestorFit
fit=nestorFit(data$Y, MO,SO, initList=initList, maxIter=5,verbatim=1,alpha = 1)
str(fit)
