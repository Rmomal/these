################
# diagnostic functions
auc<-function(pred,label){ #require(ROCR)
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}
accppvtpr<-function(probs,omega,h, seuil=0.5){
  PPV=round(sum((omega!=0)*(probs>seuil))/(sum((omega!=0)*(probs>seuil))+ sum((omega==0)*(probs>seuil))),2)#TP/(TP+FP)
  PPVH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/(sum((omega[h,]!=0)*(probs[h,]>seuil))+ sum((omega[h,]==0)*(probs[h,]>seuil))),2)
  PPVO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))+  sum((omega[-h,-h]==0)*(probs[-h,-h]>seuil))),2)
  
  TPR=round(sum((omega!=0)*(probs>seuil))/sum(omega!=0), 2)
  TPRH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/sum(omega[h,]!=0), 2)
  TPRO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/sum(omega[-h,-h]!=0), 2)
  P=(omega[h,]!=0) ; N=(omega[h,]==0)
  FP=sum((probs[h,]>seuil)*(omega[h,]==0))
  TN=sum((probs[h,]<seuil)*(omega[h,]==0))
  FPRH=FP/(FP+TN)
  FNRH=1-TPRH
  return(c(PPV,PPVH,PPVO,TPR,TPRH,TPRO,FPRH,FNRH))
}
get_PPV_TPR<-function(probs, omega, seuil){
  PPV=round(sum((omega!=0)*(probs>seuil))/(sum((omega!=0)*(probs>seuil))+ sum((omega==0)*(probs>seuil))),2)#TP/(TP+FP)
  TPR=round(sum((omega!=0)*(probs>seuil))/sum(omega!=0), 2)
  return(c(PPV=PPV, TPR=TPR))
}
seq_PPV_TPR<-function(probs,omega,seq_seuil=seq(0,max(probs),max(probs)/50)){
  tmp=sapply(seq_seuil,function(x)  get_PPV_TPR(seuil=x,probs=probs,omega=omega))
  res=data.frame(cbind(t(tmp),seq_seuil))
  colnames(res)=c("PPV","TPR","seuil")
  return(as_tibble(res))
}
courbes_seuil<-function(probs,omega,h,seq_seuil){
  tmp=sapply(seq_seuil,function(x)  accppvtpr(seuil=x,probs=probs,omega=omega,h=h))
  res=data.frame(cbind(t(tmp),seq_seuil))
  colnames(res)=c("Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","FPRH","FNRH","seuil")
  return(as_tibble(res))
}
compute_nSNR<-function(K, indexmissing){
  H=indexmissing ; p=ncol(K)
  O=(1:p)[-H]
  num=norm(K[O,H]%*%solve(K[O,O])%*%K[H,O],type='F')^2
  denom=norm(K[H,H], type='F')^2
  return(num/denom)
}
computeFPN<-function(cliqueList, trueClique,p){
  N=setdiff(1:p,trueClique)
  FP=unlist(lapply(cliqueList, function(init){init=init[[1]]
  sum(init%in%N)/length(N)}))
  FN=unlist(lapply(cliqueList, function(init){init=init[[1]]
    sum(setdiff(1:p,init)%in%trueClique)/length(trueClique)}))
  return(data.frame(FP=FP, FN=FN))
}
compute_TPPV<-function(cliqueList, trueClique,p){
  N=setdiff(1:p,trueClique)
  TP= sum(init%in%trueClique)
  FN= sum(setdiff(1:p,init)%in%trueClique)
  FP= sum(init%in%N)
  TPR = TP/(TP+FN)
  PPV=TP/(TP+FP)
  return(data.frame(TPRH=TPR, PPVH=PPV))
}
################
# graphics
ggimage<-function(data){
  melted_data <- melt(data)
  ggplot(melted_data, aes(x=Var1, y=Var2, fill=value)) + theme_light()+labs(x="",y="")+
    geom_tile() +guides(fill=FALSE)+ theme(plot.title = element_text(size=10, hjust=0.5))+ coord_fixed()
}
plotVEM<-function(probs,omega,r,seuil){
  # plots heatmaps for the chosen threshold and print verdicts as title
  q=ncol(omega)
  h=(q-r):q
  performance=accppvtpr(probs,omega,h,seuil)
  PPV=performance[1] ;PPVH=performance[2] ; PPVO=performance[3]
  TPR=performance[4] ;TPRH=performance[5] ;TPRO=performance[6] 
  p1<-ggimage(probs)+labs(title=paste0("G hat"))
  p2<-ggimage(omega)+labs(title="True G")
  auc<-round(auc(pred = probs, label = omega),3)
  grid.arrange(p1,p2,ncol=2, top=paste0("Recall=",TPR," (Obs=",TPRO," , Hid=",TPRH,
                                        ")\n Precision=",PPV," (Obs=",PPVO," , Hid=",PPVH,")",
                                        "\n AUC=",auc))
}


plotVerdict<-function(values,colonne){
  # plots verdict curves along threshold from VEM result
  colonne<-enquo(colonne)
  values %>%as_tibble() %>% gather(key,value, - !!colonne)%>%
    mutate(status=ifelse(substr(key,4,4)=="","all",substr(key,4,4))) %>% 
    mutate(key=substr(key,1,3)) %>% spread(key,value)%>% 
    mutate(PPV=ifelse(PPV<0,NA,PPV)) %>% 
    ggplot(aes(TPR,PPV,color=status))+
    geom_rect(aes(xmin=0.5, xmax=1, ymin=0.5, ymax=1), fill="gray90", color="gray90",alpha=0.1)+
    geom_point()+  geom_line()+
    facet_wrap(~status)+mytheme.dark("")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
mytheme.light <- list(theme_light(), scale_color_brewer("",palette="Set3"),guides(color=FALSE),
                      theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                            plot.title = element_text(hjust = 0.5)))

mytheme.dark <-function(legend){
  list= list(theme_light(), scale_color_brewer(legend,palette="Dark2"), 
             scale_fill_brewer(legend,palette="Dark2"), 
             theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                   plot.title = element_text(hjust = 0.5)))
  return(list)}
mytheme <- list(theme_light(), scale_fill_brewer("",palette="Dark2"),#scale_colour_hp_d(option = "RavenClaw", name = ""), #scale_color_manual("",values=mypal),
                theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                      plot.title = element_text(hjust = 0.5)))

############
# selection criteria
VBIC<-function(TrueJ,p,r,d,n){
  q=p+r
  nbparam<-p*(d) + (p*(p+1)/2 +r*p)+(q*(q-1)/2 - 1) #d comprends l'intercept
  vBIC=TrueJ - nbparam*log(n)/2
  return(vBIC)
}

ICL_T<-function(TrueJ, Wg,Pg,n,r){
  pen_T=-( sum( Pg * log(Wg) ) - logSumTree(Wg)$det )
  # if(pen_T<0) browser()
  ICL=TrueJ - pen_T 
  return(ICL)
}
ICL_ZH<-function(TrueJ,S,n,r){
  q=ncol(S) ; p=q-r
  H=(p+1):q
  if(r!=0){
    pen_ZH= 0.5*sum(log(S[,H])) + n*r*0.5*(1+log(2*pi))
  }else{
    pen_ZH= 0
  }
  ICL=TrueJ - pen_ZH
  return(ICL)
}

ICL<-function(TrueJ, Pg,Wg ,S,n,r,d, omega){
  q=ncol(S) ; p=q-r
  H=(p+1):q
  if(r!=0){
    pen_ZH= 0.5*sum(log(S[,H])) + n*r*0.5*(1+log(2*pi))
  }else{
    pen_ZH= 0
  }
  Pg=Kirshner(Wg)
  pen_T=-( sum( Pg * log(Wg+(Wg==0)) ) - logSumTree(Wg)$det) 
  pen_r<-p*(d) + (p*(p+1)/2 +r*p+r)+(q*(q-1)/2 - 1) #d comprends l'intercept
  # norm=n*q*(q-1)/2
  
  
  ICL=TrueJ   - (pen_T + pen_ZH + pen_r*log(n)/2)
  return(ICL)
}

criteria<-function(List.vem,counts,theta, matcovar,r){
  n=nrow(counts);p= ncol(counts)
  data<-lapply(List.vem,function(vem){
    J<-True_lowBound(counts,vem$M,vem$S, theta, matcovar,vem$W, vem$Wg, 
                     vem$Pg, vem$omega )
    vBIC<-VBIC(J,p,r=r, d=ncol(matcovar), n=n)
    # ICLT<-ICL_T(J, vem$Wg,vem$Pg,n,r)
    # ICLZH<-ICL_ZH(J,vem$S, n,r)
    ICL<-ICL(J, vem$Pg,vem$Wg,vem$S, n,r,d=ncol(matcovar), omega=vem$omega)
    res=data.frame(vBIC=vBIC, ICL=ICL,J)
    return(res)
  })
  res=do.call(rbind,data)
  res$r=r
  return(res)
}
cond_lap<-function(W,index){
  L=Laplacian(W)[-index,-index]
  lambda=svd(L)$d
  cond=min(abs(lambda))/max(abs(lambda))
  return(cond)
}
init.mclust<-function(S,nb.missing=1, n.noise=50,plot=TRUE, title="",trueClique=NULL){
  Scomp=prcomp(S,scale. = TRUE)
  data=data.frame(Scomp$rotation[,1:2]%*%diag(Scomp$sdev[1:2]))
  datapolar=cart2pol(x=data[,1],y=data[,2])[,1:2] # recup r et theta
  datapolar_half=datapolar %>% mutate(theta2=ifelse(theta>pi,theta-pi,theta)) %>%
    dplyr::select(r,theta2)
  
  # noise in polar coords
  r <- sqrt(runif(n.noise))
  theta2 <- runif(n.noise, 0, pi)
  datapolarall=rbind(datapolar_half,cbind(r,theta2)) 
  newdata=pol2cart(datapolarall$r,datapolarall$theta2)[,1:2]#plot(newdata,col=col, xlim=c(-1,1),ylim=c(-1,1), pch=20)
  noiseInit<-sample(c(T,F), size=ncol(S), replace=T, prob=c(3, 1)) # noiseInit<-c(rep(F,ncol(S)),rep(T,n.noise))
  
  clust= tryCatch({
    Mclust(data=newdata, initialization = list(noise=noiseInit),  G=nb.missing, verbose = FALSE)
  }, error = function(e) {
    message("new noise")
    r <- sqrt(runif(n.noise))
    theta2 <- runif(n.noise, 0, pi)
    datapolarall=rbind(datapolar_half,cbind(r,theta2)) 
    newdata=pol2cart(datapolarall$r,datapolarall$theta2)[,1:2]
    Mclust(data=newdata, initialization = list(noise=noiseInit),  G=nb.missing, verbose = FALSE)
  }, finally = { })
  groups<-mclust::map(clust$z)[1:nrow(S)]
  res<-lapply(1:nb.missing, function(c){
    which(groups==c)
  })
  
  return(list(init=res, data=data))
}
plotInitMclust<-function(res, title){
  tmp=sapply(res$init, function(c){
    res=1*(1:p) %in% c
  })
  tmp=tmp %>% as_tibble() %>% mutate(null=ifelse(rowSums(tmp)==0,1,0))
  colors= as.matrix((tmp))%*%matrix(1:ncol(tmp), ncol(tmp), 1)
  g= ggplot(res$data,aes(X1,X2, label=1:p,color=as.factor(colors)))+geom_point(size=0.1)+
    theme_light()+geom_text()+labs(x="eig vect 1",y="eig vect 2", title=title)+
    guides(color=FALSE)+scale_color_brewer(palette="Dark2")+
    geom_hline(yintercept=0, color="gray50")+geom_vline(xintercept=0, color="gray50")
  print(g)
  
}

# Variable clustering from hierachical PCA

###############################################################################
WhichMinMat <- function(A){
  jmin = apply(A, 1, which.min)
  imin = which.min(sapply(1:nrow(A), function(i){A[i, jmin[i]]}))
  jmin = jmin[imin]
  return(c(imin, jmin))
}

###############################################################################
MinMat <- function(A){ijmin = WhichMinMat(A); return(A[ijmin[1], ijmin[2]])}

###############################################################################
F_VarClustPCA <- function(S,traceS=FALSE){
  p = ncol(S)
  
  # Cost matrix
  C = matrix(Inf, p, p); PC = array(dim=c(n, p, p))
  sapply(1:(p-1), function(j){sapply((j+1):p, function(k){
    C[j, k] <<- eigen(S[c(j, k), c(j, k)])$values[2]
  })})
  # image(1:p, 1:p, C)
  
  # Hierarchical clustering
  if(traceS){listS = list()}
  Stmp = S; Ctmp = C; 
  clustPath = matrix(0, (p-1), 8); colnames(clustPath) = c('j', 'k', 'j.num', 'k.num', 'clust.num', 'coef.j', 'coef.k', 'cost')
  clustContent = as.list(1:p)
  varNum = 1:p; step = 0; 
  while(step < p-2){
    step = step + 1; 
    # Pair to merge
    jkmin = sort(WhichMinMat(Ctmp)); 
    clustPath[step, 1:2] = varNum[jkmin]; 
    clustPath[step, 3:4] = jkmin; 
    clustPath[step, 5] = p+step; 
    clustContent[[p+step]] = sort(c(clustContent[[varNum[jkmin][1]]], clustContent[[varNum[jkmin][2]]]))
    # Update covariance
    eigenSjk = eigen(Stmp[jkmin, jkmin]); 
    clustPath[step, 6:7] = eigenSjk$vectors[, 1]
    Spp = Stmp[, jkmin]%*%eigenSjk$vectors[, 1]
    Stmp = rbind(cbind(Stmp, Spp), cbind(t(Spp), eigenSjk$values[1])); 
    # Update distances
    Cpp = sapply(1:(ncol(Stmp)-1), function(j){
      eigen(Stmp[c(j, ncol(Stmp)), c(j, ncol(Stmp))])$values[2]
    }); 
    Ctmp = rbind(cbind(Ctmp, Cpp), rep(Inf, (ncol(Ctmp)+1))); 
    clustPath[step, 8] = eigenSjk$values[2]
    # Removing merges rows and columns
    Stmp = Stmp[-jkmin[2], ]; Stmp = Stmp[, -jkmin[2]]; 
    Stmp = Stmp[-jkmin[1], ]; Stmp = Stmp[, -jkmin[1]]; 
    if(traceS){listS[[step]]=Stmp}
    Ctmp = as.matrix(Ctmp)[-jkmin[2], ]; Ctmp = as.matrix(Ctmp)[, -jkmin[2]]; 
    Ctmp = as.matrix(Ctmp)[-jkmin[1], ]; Ctmp = as.matrix(Ctmp)[, -jkmin[1]]; 
    varNum = varNum[-jkmin]; varNum[p-step] = p+step
    #cat(step, ':', clustPath[step, ], '\n')
    # image(Ctmp, main=paste(step, dim(Ctmp)[1]), xlab='', ylab='')
  }
  # Last step
  eigenSjk = eigen(Stmp); 
  clustPath[p-1, ] = c(varNum, c(1, 2), 2*p-1, eigenSjk$vectors[, 1], eigenSjk$values[2])
  clustContent[[2*p-1]] = sort(c(clustContent[[varNum[jkmin][1]]], clustContent[[varNum[jkmin][2]]]))
  clustPath = as.data.frame(clustPath); clustPath$height = cumsum(clustPath$cost)
  
  # Clustering matrix
  clustMatrix = (1:p) %o% rep(1, (p-1))
  sapply(1:(p-1), function(h){
    if(h>1){clustMatrix[, h] <<- clustMatrix[, (h-1)]}
    clustMatrix[clustContent[[p+h]], h] <<- p+h
  })
  
  # Merging path 
  clustMerge = clustPath[, 1:2]
  sapply(1:2, function(c){
    clustMerge[which(clustMerge[, c] <= p), c] <<- -clustMerge[which(clustMerge[, c] <= p), c]
    clustMerge[which(clustMerge[, c] > p), c] <<- clustMerge[which(clustMerge[, c] > p), c] - p
  })
  
  
  return(list(clustPath=clustPath, clustContent=clustContent, clustMatrix=clustMatrix, 
              lastCost=eigenSjk$values[1], clustMerge=clustMerge))
}

rSpanTreeV1 <- function(beta, prob){
  # Approximate sampling of a spanning according to edge probabilities
  # Rejection sampling approach
  # !!! Beta's must be normalized !!! (constant B = 1)
  # Max ratio between proposal and target heuristicaly set to p^2
  p <- nrow(prob); P <- p*(p-1)/2
  prob <- prob / max(prob) # To enforce connectivity
  probVec <- F_Sym2Vec(prob); betaVec <- F_Sym2Vec(beta)
  # # Optimal constant
  # w <- log(prob/beta); optTreeVec <- F_Sym2Vec(mst(w)); 
  # M <- p^(p-2) / exp(sum(F_Sym2Vec(w)[which(optTreeVec==1)]))
  M <- p
  OK <- FALSE; tries <- 0; pTree <- qTree <- rep(0, 1e4)
  while(!OK){
    tries <- tries + 1
    
    # graph = heterogeneous ER connected graph
    graph <- F_Vec2Sym(matrix(rbinom(P, 1, probVec)))
    while(!is.connected(graph)){graph <- F_Vec2Sym(matrix(rbinom(P, 1, probVec)))}
    # gplot(graph, gmode='graph', main=round(SumTree(graph)))
    graphVec <- F_Sym2Vec(graph)
    # qGraph <- prod(dbinom(graphVec, 1, probVec))
    
    # tree = uniform spanning tree
    w <- F_Vec2Sym(F_Sym2Vec(matrix(runif(p^2), p, p)))
    w[which(graph==0)] <- Inf
    tree <- mst(w)
    treeVec <- F_Sym2Vec(tree)
    
    # Monte-Carlo sample of graphs containing the tree
    biasedProbVec <- probVec; biasedProbVec[which(treeVec==1)] <- 1
    invSpanTreeMean <- 0; G <- 1e3
    for(g in 1:G){invSpanTreeMean <- invSpanTreeMean + 1/SumTree(F_Vec2Sym(matrix(rbinom(P, 1, biasedProbVec))))}
    invSpanTreeMean <- invSpanTreeMean/G
    qTree[tries] <- prod(probVec[which(treeVec==1)]) * invSpanTreeMean
    
    # qTree[tries] <- qGraph / SumTree(graph)
    # Actually qTree should be qGraph / SumTree(graph)
    # qTree[tries] <- prod(probVec[which(treeVec==1)]) / SumTree(graph)
    pTree[tries] <- prod(betaVec[which(treeVec==1)])
    
    if(pTree[tries] > (M * qTree[tries])){cat('[*] ')}
    if(runif(1) < pTree[tries] / (M * qTree[tries])){OK <- TRUE}
  }
  pTree <- pTree[1:tries]; qTree <- qTree[1:tries]
  cat(tries, '')
  
  return(list(tree=tree, qTree=qTree, pTree=pTree, qTree=qTree, M=M))
}

##################################################################
F_Vec2Sym <- function(A.vec){
  # Makes a symmetric matrix from the vector made of its lower tirangular part
  n = (1+sqrt(1+8*length(A.vec)))/2
  A.mat = matrix(0, n, n)
  A.mat[lower.tri(A.mat)] = A.vec
  A.mat = A.mat + t(A.mat)
  return(A.mat)
}

##################################################################
F_Sym2Vec <- function(A.mat){
  # Makes a vector from the lower triangular par of a symmetric matrix
  return(A.mat[lower.tri(A.mat)])
}
