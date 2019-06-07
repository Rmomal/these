
##########
# Settings
##########
path<-"/Users/raphaellemomal/simulations/experiment/"
type<-"erdos"
p<-30
dens<-4/p
prob=dens
n=50
B=100
# covar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),
#                   X3=(runif(n)))
# saveRDS(covar,paste0(path,"covar_n",n,".rds"))
dat <- data_from_scratch(type, p=p,r=10, covar)
edgesOrigin<-ifelse(abs(F_Sym2Vec(dat[[2]]))<1e-16,0,1)
Y<-dat[[1]]
saveRDS(Y,paste0(path,"Ymint_",type,"_",difficulty,nbgraph,".rds"))

##########
# functions
##########


m<- model.matrix(~X1+X2+X3,covar)# on choisi de mettre la constante dans PLN 
T1<-Sys.time()
pmat<-F_ResampleTreePLN(Y, X=m, O=matrix(0,nrow=nrow(covar),ncol=p), v=0.8, B=B, maxIter=20,
                        cond.tol=1e-12,cores=3)
T2<-Sys.time()
difftime(T2,T1)

pMoy<-colSums(pmat)/B
pSum<-colSums(pmat)
##################
# post treatment  
##################
# tables
doubleThresh<-1*(ifelse(colSums(ifelse(pmat<2/p,0,1))/B >0.8,1,0))
pMoyThresh<-1*(ifelse(pMoy>2/p,1,0))
pSumThresh<-1*(ifelse(pSum>2*B*0.8/p,1,0))

table(doubleThresh,edgesOrigin)
table(pMoyThresh,edgesOrigin)
table(pSumThresh,edgesOrigin)

table(pMoyThresh,doubleThresh)
table(pSumThresh,doubleThresh,edgesOrigin)
table(pSumThresh,pMoyThresh,edgesOrigin)

# graph densities
pMoy %>%as_tibble() %>%  gather(key,value) %>% rowid_to_column() %>% mutate(origin=edgesOrigin, double=doubleThresh) %>% 
  mutate(double=ifelse(double==1 ,double*value,NA),origin=ifelse(origin==1,origin*value,NA),
         selected=ifelse(value<2/p,NA,value)) %>% 
  gather(type,value,-rowid) %>% 
  ggplot(aes(value, color=type))+theme_minimal()+
  geom_vline(xintercept = 2/p)+
  geom_density()+geom_rug()+ coord_cartesian(ylim=c(0,3))+
  scale_color_manual("Mean \nprobabilities of:",values=c("#7CAE00","#F8766D","#C77CFF","#00BFC4"),
                     label=c("Double \nTreshold","Origin","Selected","Mean \nprobability"))





