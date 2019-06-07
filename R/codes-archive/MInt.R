
# require Kit, ROCR, MInt
library(parallel)
library(ROCR)
library(grid)
library(gridExtra)

parameters<-list(d=c(seq(10,30,2)),n=c(seq(20,120,10)), prob=c(seq(1,5,0.5)/20), 
                 r=c(seq(1,50,5)), dens=c(seq(1,5,0.5)/20))
B.resample<-500
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"

eval_store_mint<-function(Y,covar,path){
  #  browser()
  data<-data_for_MInt(Y,covar,path)
  x <- data[["x"]]
  y <- data[["y"]]
  T1<-Sys.time()
  m <- mint(y,x,fmla = ~feature1 + feature2+ feature3)
  m <- estimate(m)
  T2<-Sys.time()
  temps<-difftime(T2,T1)
  pred<-m$param$P
  #obs<-readRDS(paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,".rds"))$omega
  #roc_curve(pred,obs),
  return(list(pred=pred,temps= temps))
}

######
# Resampling
######
resampling_mint<-function(Y,X,v,B,path,cores){
  n = nrow(Y); p = ncol(Y);  V = round(v*n)
  
  obj<-mclapply(1:B,function(b){
    cat('\n', b, '')
    set.seed(b)
    #browser()
    sample = sample(1:n, V, replace = F)
    Y.sample = Y[sample,]
    X.sample = X[sample,]
    inf<-eval_store_mint(Y.sample,X.sample,path)
    return(inf)
  },mc.cores=cores)
  Pmat<-do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$pred)}))
  ListTmp = do.call(c,lapply(obj,function(x){x$temps}))
  return(list(Pmat=Pmat,temps=ListTmp))
}
resampling_gcoda_spiec<-function(Y,v,B,path,cores,covar, gcoda=TRUE){
  n = nrow(Y); p = ncol(Y);  V = round(v*n)
  
  obj<-mclapply(1:B,function(b){
    #browser()
    cat('\n', b, '')
    set.seed(b)
    #browser()
    sample = sample(1:n, V, replace = F)
    Y = Y[sample,]
    covar = covar[sample,]
    if(gcoda){
      out_gcodaResid<-gcoda(Y, counts=T, covar=covar)
     # K.score <- Reduce("+",out_gcodaResid$path)
      scores<-out_gcodaResid$opt.icov
    }else{
      U<-t(clr.matrix(Y,mar=1))
      m<- model.matrix(~X1+X2+X3,covar)
      model<-lm(U~m)
      inf<- spiec.easi(model$residuals, icov.select = FALSE, nlambda = 50, verbose = FALSE)
      K.score <- Reduce("+",inf$est$path)
    } 
   # scores<- K.score / max(K.score)
    return(F_Sym2Vec( scores))
  },mc.cores=cores)
  
  Pmat<-do.call(rbind,obj)
  return(Pmat)
}
type<-"cluster"
variable<-"dens"
nbgraph<-20
cval<-c(0.25)
n=100
covar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),
                  X3=(runif(n)))
saveRDS(covar,paste0(path,"mint_covar.rds"))
T1<-Sys.time()
for(valeur in cval){
  Y<-data_from_stored_graphs(type, variable, nbgraph, valeur, covar=covar,path,fixe=FALSE)
  saveRDS(Y,paste0(path,"Ymint_",variable,"_",valeur,".rds"))
  Y<-readRDS(paste0(path,"Ymint_",variable,"_",valeur,".rds"))[[1]]
  
  T1<-Sys.time()
  mint<-resampling_mint(Y,covar,v=0.8,B=B.resample,path=path, cores=1)
  T2<-Sys.time()
  mintime<-difftime(T2,T1)
  saveRDS(mint,paste0(path,valeur,"_Pmat_time_mint.rds"))
  m<- model.matrix(~X1+X2+X3,covar)[,-1]
  resEM<-F_ResampleTreePLN(Y, X=m, O=matrix(0,nrow=100,ncol=20), v=0.8, B=500, maxIter=10, cond.tol=1e-14)
  saveRDS(resEM,paste0(path,valeur,"_Pmat_time_EM.rds"))
}
T2<-Sys.time()
difftime(T2,T1)

U<-t(clr.matrix(Y,mar=1))
m<- model.matrix(~X1+X2+X3,covar)
model<-lm(U~m)
inf<- spiec.easi(model$residuals, method='mb', lambda.min.ratio=1e-8, nlambda=50,
                 icov.select.params = list(rep.num=500))
huge::huge.roc(inf$est$path, omegaOriginal, verbose=FALSE)
stars.pr(getOptMerge(inf), omegaOriginal, verbose=FALSE)
spiec_scores<-getOptMerge(inf)
spiec_scores<-resampling_gcoda_spiec(Y,v=0.8,B=500,cores=3,covar=covar,gcoda=FALSE)
gcoda_scores<-resampling_gcoda_spiec(Y,v=0.8,B=500,cores=3,covar=covar)

##############
##############

resMInt<- readRDS(paste0(path,valeur,"_Pmat_time_mint.rds"))
resEM<- readRDS(paste0(path,valeur,"_Pmat_time_EM.rds"))
omegaOriginal <- readRDS("/Users/raphaellemomal/simulations/Simu/PLN.2.0/cluster/dens/Sets_param/Graph20_0.25.rds")$omega
edgesOrigin<-ifelse(abs(F_Sym2Vec(omegaOriginal))<1e-16,0,1)
# plot(density(as.numeric(resMInt$temps)))
# mean(c(resMInt$temps))

threshMInt<-ifelse(abs(resMInt$Pmat)<1e-16,0,1) #0 si nul
threshEM<-ifelse(resEM$Pmat<0.1,0,1) #0 si <2/p
threshgCoda<-ifelse(abs(gcoda_scores)<1e-16,0,1)
fregCoda<-colSums(threshgCoda)
freEM<-colSums(threshEM)
freMInt<-colSums(threshMInt)
datafreq<-data.frame(EM=freEM,MInt=freMInt,color=as.factor(edgesOrigin+1))

visufreqs<-function(datafreq){
  pal<-c("steelblue4","orangered2")
  p1<-ggplot(datafreq,aes(freEM,freMInt,color=color))+geom_point()+
    scale_color_manual("True edge",values=pal,breaks=c(2,1),labels=c("yes","no"))+
    theme_minimal()+ geom_abline()+geom_abline(slope=0,intercept = 400,linetype="dashed")+
    geom_vline(xintercept=400,linetype="dashed")+
    labs(x="EM frequences",y="MInt frequences")
  
  p2<-datafreq %>% 
    rowid_to_column() %>% 
    ggplot(aes(rowid,EM,color=color))+geom_point(show.legend=FALSE)+
    scale_color_manual("True edge",values=pal,breaks=c(2,1),labels=c("yes","no"))+
    theme_minimal()+ geom_hline(yintercept = 400,linetype="dashed")+
    labs(x="Edge index")
  p3<-datafreq %>% 
    rowid_to_column() %>% 
    ggplot(aes(rowid,MInt,color=color))+geom_point(show.legend=FALSE)+
    scale_color_manual("True edge",values=pal,breaks=c(2,1),labels=c("yes","no"))+
    theme_minimal()+ geom_hline(yintercept = 400,linetype="dashed")+
    labs(x="Edge index")
  p1/(p2+p3)
}
visufreqs(datafreq)
##############
##############

##############
# essai ROC
#


prediction_MInt<-prediction(datafreq$MInt/500,as.numeric(datafreq$color)-1)
prediction_EM<-prediction(datafreq$EM/500,as.numeric(datafreq$color)-1)
prediction_spiec<-prediction(F_Sym2Vec( spiec_scores)/max(spiec_scores),as.numeric(datafreq$color)-1)
prediction_gcoda<-prediction(fregCoda/max(fregCoda),as.numeric(datafreq$color)-1)

# courbes fdr
fdrmint <- data.frame(x=performance(prediction_MInt,"pcfall")@x.values[[1]],y=performance(prediction_MInt,"pcfall")@y.values[[1]])
fdrem <- data.frame(x=performance(prediction_EM,"pcfall")@x.values[[1]],y=performance(prediction_EM,"pcfall")@y.values[[1]])
fdrgcoda <- data.frame(x=performance(prediction_gcoda,"pcfall")@x.values[[1]],y=performance(prediction_gcoda,"pcfall")@y.values[[1]])

fdr<-rbind(cbind(fdrmint,method="MInt"),cbind(fdrem,method="EM"))#,cbind(fdrgcoda,method="gCoda"))
#courbes precrec
precrec_mint <- data.frame(x=performance(prediction_MInt,"rec")@y.values[[1]],y=performance(prediction_MInt,"prec")@y.values[[1]])
precrec_em <- data.frame(x=performance(prediction_EM,"rec")@y.values[[1]],y=performance(prediction_EM,"prec")@y.values[[1]])
precrec_spiec<- data.frame(x=performance(prediction_spiec,"rec")@y.values[[1]],y=performance(prediction_spiec,"prec")@y.values[[1]])
precrec_gcoda<- data.frame(x=performance(prediction_gcoda,"rec")@y.values[[1]],y=performance(prediction_gcoda,"prec")@y.values[[1]])

precrec<-rbind(cbind(precrec_mint,method="MInt"),cbind(precrec_em,method="EM"))#, cbind(precrec_gcoda,method="gcoda"))
#courbes roc
roc_mint <- data.frame(y=performance(prediction_MInt,"rec")@y.values[[1]],x=performance(prediction_MInt,"fpr")@y.values[[1]])
roc_em <- data.frame(y=performance(prediction_EM,"rec")@y.values[[1]],x=performance(prediction_EM,"fpr")@y.values[[1]])
roc_spiec <- data.frame(y=performance(prediction_spiec,"rec")@y.values[[1]],x=performance(prediction_spiec,"fpr")@y.values[[1]])
roc_gcoda <- data.frame(y=performance(prediction_gcoda,"rec")@y.values[[1]],x=performance(prediction_gcoda,"fpr")@y.values[[1]])

roc<-rbind(cbind(roc_mint,method="MInt"),cbind(roc_em,method="EM"))#,cbind(roc_gcoda,method="gcoda"))

pal<-c("#6abd35","steelblue4","indianred")
p1<-ggplot(fdr,aes(x,y,color=method))+
  geom_point()+geom_line()+theme_minimal()+
  scale_color_manual("Method:",values = pal)+
  labs(x="Edge frequence threshold (%)",y="FDR (%)")+
  theme(legend.text = element_text(size=13))

p2<-ggplot(precrec,aes(x,y,color=method))+
  geom_point()+geom_line()+theme_minimal()+
  scale_color_manual("Method:",values=pal)+labs(x="Recall",y="Precision")+
  theme(legend.text = element_text(size=13))

p3<-ggplot(roc,aes(x,y,color=method))+
  geom_point()+geom_line()+theme_minimal()+
  scale_color_manual("Method:",values=pal)+labs(x="FPR",y="Recall")+
  theme(legend.text = element_text(size=13))


grid_arrange_shared_legend(p1,p2,p3,nrow=1,ncol=3)


