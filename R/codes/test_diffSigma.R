# effondrement sigma

dev.off()
effondrementSigma<-function(data,dataShuffle, var){
  sigma<-PLN(data ~ var)$model_par$Sigma
  sigma_shuffle<-PLN(dataShuffle ~ var)$model_par$Sigma
  diffSigma=c(sigma)-c(sigma_shuffle)
  par(mfrow=c(1,2), oma=c(0,0,2,0))
  hist(diffSigma,prob=TRUE,main=" ")
  curve(dnorm(x,mean(diffSigma),sd(diffSigma)),col="blue", add=TRUE)
  qqnorm(diffSigma, cex=0.3, col="blue", main=" "); abline(mean(diffSigma),sd(diffSigma))
  title(paste0("q1=",round(quantile(diffSigma,0.25),2), ", q2=",round(median(diffSigma),2), ", q3=",
               round(quantile(diffSigma,0.75),2), ",\n mean=",round(mean(diffSigma),2),", sd=", round(sd(diffSigma,0.25),2)), outer=TRUE)
}
effondrementSigma(counts,data_shuffleSite, covar$site)
effondrementSigma(Y,data_shuffleTree, X$tree)

############
shuffle<-function(counts,var){
  levels<-levels(var)
  counts_shuffle=0*counts
  for (bloc in levels){
    indices<-which(var==bloc)
    counts_shuffle[indices,]<-apply(counts[indices,], 2, function(x){
      sample(x)
    })
  }
  return(counts_shuffle)
}
# shuffled data
data_shuffleSite<-shuffle(counts,covar$site)
data_shuffleTree<-shuffle(Y,X$tree)
