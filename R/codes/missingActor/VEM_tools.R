################
# diagnostic functions
auc<-function(pred,label){ #require(ROCR)
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  ROC_auc <- round(performance(prediction,"auc")@y.values[[1]],digits=2)
  return(ROC_auc)
}
accppvtpr<-function(probs,omega,h, seuil=0.5){
  Acc=round(mean(1*(probs>seuil)==omega),2) #(TP+TN)/(TP+FP+TN+FN)
  AccH=round(mean(1*(probs[h,]>seuil)==omega[h,]),2)
  AccO=round(mean(1*(probs[-h,-h]>seuil)==omega[-h,-h]),2)
  PPV=round(sum((omega!=0)*(probs>seuil))/(sum((omega!=0)*(probs>seuil))+ sum((omega==0)*(probs>seuil))),2)#TP/(TP+FP)
  PPVH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/(sum((omega[h,]!=0)*(probs[h,]>seuil))+ sum((omega[h,]==0)*(probs[h,]>seuil))),2)
  PPVO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))+  sum((omega[-h,-h]==0)*(probs[-h,-h]>seuil))),2)
  
  TPR=round(sum((omega!=0)*(probs>seuil))/sum(omega!=0), 2)
  TPRH=round(sum((omega[h,]!=0)*(probs[h,]>seuil))/sum(omega[h,]!=0), 2)
  TPRO=round(sum((omega[-h,-h]!=0)*(probs[-h,-h]>seuil))/sum(omega[-h,-h]!=0), 2)
  
  return(c(Acc,AccH,AccO,PPV,PPVH,PPVO,TPR,TPRH,TPRO))
}
courbes_seuil<-function(probs,omega,h,seq_seuil){
  tmp=sapply(seq_seuil,function(x)  accppvtpr(seuil=x,probs=probs,omega=omega,h=h))
  res=data.frame(cbind(t(tmp),seq_seuil))
  colnames(res)=c("Acc","AccH","AccO","PPV","PPVH","PPVO","TPR","TPRH","TPRO","seuil")
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
  FP=unlist(lapply(cliqueList, function(init){sum(init%in%N)/length(N)}))
  FN=unlist(lapply(cliqueList, function(init){
    sum(setdiff(1:p,init)%in%trueClique)/length(trueClique)}))
  return(data.frame(FP=FP, FN=FN))
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
  Acc=performance[1] ;AccH=performance[2] ;AccO=performance[3] 
  PPV=performance[4] ;PPVH=performance[5] ; PPVO=performance[6]
  TPR=performance[7] ;TPRH=performance[8] ;TPRO=performance[9] 
  p1<-ggimage(probs)+labs(title=paste0("G hat"))
  p2<-ggimage(omega)+labs(title="True G")
  auc<-round(auc(pred = probs, label = omega),3)
  grid.arrange(p1,p2,ncol=2, top=paste0("Tpr=",TPR," (TprO=",TPRO," , TprH=",TPRH,
                                        ")\n Ppv=",PPV," (PpvO=",PPVO," , PpvH=",PPVH,")",
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