# inférer le réseau avec 1 variable manquante et voir si les moyennes correspondent à celles de la température
library(PLNmodels)
library(EMtree)
# DATA
data.dir = "/Users/raphaellemomal/these/R/codes/Spatial/"
data.name = 'BarentsFish'
load(paste0(data.dir, data.name, '.Rdata'))
n=nrow(Data$count)
#PLN, EMtree
fitPLN = PLN(Data$count~1)
sigma=fitPLN$model_par$Sigma
fitEM=EMtree(fitPLN,maxIter = 50,cond.tol = 1e-6,plot = TRUE,verbatim = TRUE)

# LiTree
cliques=findCliques(sigma,1)
initial.param<-initEM(sigma,n=n,cliquelist = cliques,pca=TRUE) # ajout de trois variables manquantes
K0=initial.param$K0
Sigma0=initial.param$Sigma0

XO=rmvnorm(n=n,sigma=sigma)
infEM=EMtreeMissing(S = sigma,k = 1,K0 =K0 ,Sigma0 = Sigma0,n=n, XO=XO, condMeans=TRUE)

vecmu=infEM$mu
Data$covariates %>% as_tibble() %>% mutate(mu=vecmu) %>% gather(covar,value,-mu) %>% 
  ggplot(aes(log(abs(mu)),value))+geom_point()+facet_wrap(~covar, scales = "free")+theme_minimal()
