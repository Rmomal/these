
sapply(1:403, function(seed){
 
  MO=tryCatch({readRDS(paste0("/Users/raphaellemomal/simulations/15nodes_V4_oracle/SF_seed",seed,".rds"))$ListVEM[[1]]$M[,-15] 
  n=nrow(MO) 
  M=t(MO)%*%MO 
  res=max(M/n) 
  print(res)},
              error=function(e){res=NULL}, finally={})
  return(res)
})
f<-function(x){
  -0.5*log(1-x^2)+x^2/(1-x^2)
}

borne_alpha<-function(x,q=15,delta=.Machine$double.xmax,n=200){
  C=f(x)
  res=((1/(q-1))*log(delta)-log(q-1))/(C*n)
  return((res))
}
curve(f, 0,0.8)
curve(borne_alpha,0.5,0.8)
borne_alpha(x=0.8)


MO=readRDS(paste0("/Users/raphaellemomal/simulations/15nodes_V4_oracle/SF_seed1.rds"))$ListVEM[[1]]$M[,-15] 
M=t(MO)%*%MO /n
res=max(abs(M/n) )
sapply(F_Sym2Vec(M), f)


plot(density(M/n))  
hist(M/n, breaks=50)
