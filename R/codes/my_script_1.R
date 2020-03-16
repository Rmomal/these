
x=c()
for(i in 1:50){
  x=c(x,mean(rnorm(10,0,2)))
}
saveRDS(x,"results.rds")