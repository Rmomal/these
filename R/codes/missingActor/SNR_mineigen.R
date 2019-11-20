source('/Users/raphaellemomal/these/R/codes/missingActor/LiTree/functions-simulation.R')
compute_SNR<-function(K, indexmissing){
  H=indexmissing ; p=ncol(K)
  O=(1:p)[-H]
  num=norm(K[O,H]%*%K[H,O]/K[H,H],type='F')^2
  denom=norm(K[O,O],type='F')^2
  return(num/denom)
}


Data=gener.data(counts = FALSE,p=p,n=n,type="scale-free")
omega=Data$omega
adjmat=1*(omega!=0)

p=30
star.graph <- graphModel$new(type = "starerdos",size=p, p.or.m = 10/p)
star.model <- GGMmodel$new(graph=star.graph,nb.missing.var= 1)
#star.model$randomSample(n=40)
star.model$missing.var.list
adjmat=star.model$getAdjmat()

missing=chooseMissingVar(adjmat,missing.var.number = 1)
B=100

L=lapply(1:B, function(x){
  sapply(list("default","gen","chris", "raph"), function(method){
    sapply(seq(0,1,0.1), function(prop){
      K=covarianceFromGraph(adjmat,missing.var.list = missing,
                            method=method,prop.positive.cor = prop ,
                            alpha.hidden = 1, alpha.observed = 1.5)
      return(compute_SNR(K,missing))
    })
  })
})

L2=lapply(L,function(x){
  data.frame(x) %>% as_tibble() %>% mutate(prop=seq(0,1,0.1))
})

L2=do.call(rbind,(L2))
colnames(L2)=c("default","gen","chris","raph","prop.pos")
L2 %>% as_tibble() %>% gather(key,value,-prop.pos) %>% 
  ggplot(aes(x=key,y=value, color=as.factor(prop.pos)))+geom_boxplot()+theme_light()+labs(x="",y="SNR")


(3^2/1)^2
##########################


Q=eigen(K)$vectors
lambda=eigen(K)$values
min(lambda)
Lambda_corrected =  diag( pmax(lambda, max(lambda)/5)) 
test=Q %*% Lambda_corrected %*% t(Q)

min(eigen(test)$values)
compute_SNR(test,30)

