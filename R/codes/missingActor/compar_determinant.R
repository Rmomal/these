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
