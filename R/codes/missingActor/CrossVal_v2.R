# Model selection through cross validation using pairwise composite likelihood


library(nestor)
library(gridExtra)
library(harrypotter)
library(tictoc)
library(sna)
library(ape)
library(poilog)
library(parallel)
library(EMtree)
library(Matrix)
library(tidyverse)
F_Vec2Sym <- function(A.vec){
  # Makes a symmetric matrix from the vector made of its lower tirangular part
  n = (1+sqrt(1+8*length(A.vec)))/2
  A.mat = matrix(0, n, n)
  A.mat[lower.tri(A.mat)] = A.vec
  A.mat = A.mat + t(A.mat)
  return(A.mat)
}
F_Sym2Vec <- function(A.mat){
  # Makes a vector from the lower triangular par of a symmetric matrix
  return(A.mat[lower.tri(A.mat)])
}
mytheme.dark <-function(legend){
  list= list(theme_light(), scale_color_brewer(legend,palette="Dark2"), 
             scale_fill_brewer(legend,palette="Dark2"), 
             theme(strip.background=element_rect(fill="gray50",colour ="gray50"),
                   plot.title = element_text(hjust = 0.5)))}
rSpanTreeV1 <- function(beta, prob,maxTries=30){
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
    c=0
    while(!is.connected(graph)){
      c=c+1
      #cat(paste0("\n---- c=",c," ----"))
      if(c<maxTries){
        graph <- F_Vec2Sym(matrix(rbinom(P, 1, probVec)))
      }else{
        probVec2=probVec; probVec2[probVec2<0.05]=0.05
        graph <- F_Vec2Sym(matrix(rbinom(P, 1, probVec2)))
      }
    }
    # gplot(graph, gmode='graph', main=round(SumTree(graph)))
    graphVec <- F_Sym2Vec(graph)
    # qGraph <- prod(dbinom(graphVec, 1, probVec))
    # tree = uniform spanning tree
    w <- F_Vec2Sym(F_Sym2Vec(matrix(runif(p^2,0,beta), p, p)))
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
    #  if(pTree[tries] > (M * qTree[tries])){cat('[*] ')}
    if(runif(1) < pTree[tries] / (M * qTree[tries])){OK <- TRUE}
  }
  pTree <- pTree[1:tries]; qTree <- qTree[1:tries]
  # cat(tries, '')
  return(list(tree=tree, qTree=qTree, pTree=pTree, qTree=qTree, M=M))
}
split_fit<-function(Y,X=NULL,v=0.8,r=1,alpha=0.1, bloc=1,eps=1e-3){#shuffle data prior to this
  p=ncol(Y) ;n=nrow(Y); O=1:p ; H=(p+1):(p+r)
  # split data
  T1<-Sys.time()
  ntrain=round(n*v) ; ntest=n-ntrain ; nblocs=round(n/ntest)
  vec_blocs=rep(0,n)
  vec_blocs[1:((nblocs-1)*ntest)]=rep(1:(nblocs-1), each=ntest)
  vec_blocs[vec_blocs==0]=nblocs
  samp=which(vec_blocs==bloc)
 
  #-- normalized PLN outputs
  Y.train<-Y[-samp,] ; Y.test=Y[samp,]
  X.train<-X[-samp,]; X.test=X[samp,]
  #normPLNfit.train<-norm_PLN(Y.train,X = X.train)
  normPLNfit.train<-norm_PLN(Y,X = X)
  MO.train<-normPLNfit.train$MO[-samp,]  ;SO.train<-normPLNfit.train$SO[-samp,]  
  sigma_obs=t(MO.train)%*%MO.train+diag(colSums(SO.train))
  ntrain=nrow(MO.train)
  
  if(r!=0){
    if(r==1){minV=2
    }else{minV=0}
    #--initial cliques and run
    # if(r==1){
    #   clique=complement_spca(MO.train, 2)
    # }else{
      # normPLNfit.all<-norm_PLN(Y,X = X)
      # MO<-normPLNfit.all$MO
      #   clique1=boot_FitSparsePCA(MO, B=200,r=r,min.size=1, cores=3)
      clique=boot_FitSparsePCA(MO.train, B=200,r=r,min.size=minV, cores=3)$cliqueList
      # clique=unique(c(clique1$cliqueList,clique2$cliqueList))
    #}
    # if(length(clique$cliqueList)>50){
    #   clique_spca<-FitSparsePCA(counts, r=2)$cliques
    #   complement=lapply(clique_spca, function(clique){setdiff(1:p,clique)})
    #   clique=list()
    #   clique$cliqueList=lapply(c(clique_spca,complement), function(cl) list(cl))
    # }
    #VEM
    cat(paste0("\n   Fitting ",length(clique)," VEM..."))
    ListVEM<-List.nestorFit(cliqueList =clique, sigma_obs, MO.train,SO.train,r=r,alpha=alpha,
                            eps=eps,maxIter=100, cores=3, trackJ=FALSE)
    #-- select best J among models which gives r Mh with reasonable sd 
    if(r>1){
      vec_sdMH_J=do.call(rbind, lapply(ListVEM, function(vem){
        if(length(vem)>5){
          J=tail(vem$lowbound$J, 1)
          sdM=apply(vem$M,2, function(x) sd(x))
          sdMH=tail(sdM,r)
          res=c(sdMH,J=J)
        }else{# if degenerates: NaN
          J=NaN ;sdMH =rep(NaN,r)
          res=c(sdMH,J=J)}
        return(res)
      }))
      data_sd=data.frame(vec_sdMH_J) %>% dplyr::select(-J)
      AllActors=apply(data_sd, 1, function(x){
        sum(log(x)>-20)==r # min threshold for log(sd)
      })
      if(sum(AllActors, na.rm=TRUE)!=0){
        index=unlist(data.frame(vec_sdMH_J,AllActors=AllActors,num=1:nrow(data_sd)) %>% as_tibble() %>% 
                       filter(AllActors) %>% filter(J==max(J,na.rm=TRUE)) %>% dplyr::select(num))
      }else{
        if(sum(!is.na(AllActors))!=0){
          index=which.max(data.frame(vec_sdMH_J)$J)
        }else{
          index=NA
        }
      }
      
    }else{
      vec_J=do.call(rbind, lapply(ListVEM, function(vem){
        if(length(vem)>5){
          J=tail(vem$lowbound$J, 1)
        }else{
          J=NaN }
      }))
      index=which.max(vec_J)
    }
    if(length(index)==0) index=NA
    VEM=ListVEM[[index]]
  }else{
    cat(paste0("\n   Fitting 1 VEM with r=0..."))
    init0=initVEM(cliqueList=initClique,sigma_obs, MO.train,r=0)
    VEM<-tryCatch({nestorFit(MO.train, SO.train,  initList=init0, eps=eps,maxIter = 100,
                             print.hist = FALSE, alpha=alpha, verbatim=FALSE, trackJ=FALSE )},
                  error=function(e){e}, finally={})
  }
  T2<-Sys.time()
  runtime=difftime(T2, T1)
  cat(paste0(round(runtime,3)," ", attr(runtime, "units")))
  if(!is.null(VEM)){
    res=list(Ytest=Y.test,X.test=X.test,Pg=VEM$Pg,beta=VEM$Wg,M=VEM$M,S=VEM$S,
             Omega_hat=VEM$Omega,D=normPLNfit.train$D, theta=normPLNfit.train$theta )
  }else{
    res=list(Ytest=Y.test,X.test=X.test,Pg=NULL,beta=NULL,M=NULL,S=NULL,Omega_hat=NULL,
             D=normPLNfit.train$D, theta=normPLNfit.train$theta )
  }
  return(res)
}

pY_condT<-function(Ytest,X.test,tree,M,S, Omega_hat,D,theta, r=1, plot=FALSE){
  #--calcul de la vraisemblance composite par paires
  p=ncol(Ytest) ;n=nrow(Ytest); O=1:p
  
  #--- calcul critere : pour un T calculer omega, sigma et p(Y)
  if(!is.null(Omega_hat)){
    OmegaT=tree*Omega_hat
    ntrain = nrow(M)
    Rho = cov2cor((t(M)%*%M+diag(colSums(S)))/ntrain)
    quantity<-tree*Rho^2/(1-Rho^2)
    diag(quantity)=NA
    vectorSuml=colSums(quantity, na.rm=TRUE)
    diag(OmegaT) = 1+vectorSuml
    attr(OmegaT,"class")="matrix"
    
    if(r!=0){
      H=(p+1):(p+r)
      Sig=solve(nearPD(OmegaT[O,O]- OmegaT[O,H]%*%solve(OmegaT[H,H])%*%OmegaT[H,O], 
                       eig.tol = 1e-14, posd.tol = 1e-14)$mat)
      SigmaTm=cov2cor(as.matrix(Sig))
    }else{
      Sig=solve(nearPD(OmegaT, eig.tol = 1e-14, posd.tol = 1e-14)$mat)
      SigmaTm=cov2cor(as.matrix(Sig))
    }
    #pour chaque site test i et chaque paire d’especes (j,k),
    #on calcule la densit́e Poissonlog-normale
    mat_dens=matrix(0,p,p)
    if(!is.null(X.test)){
      sapply(1:(p-1), function(i){
        sapply((i+1):p, function(j){
          sig=sqrt(D[i]) ; sig2=sqrt(D[j]) ; rho=SigmaTm[i,j]
          if(rho<0) rho = max(SigmaTm[i,j],-0.9999)
          if(rho>0) rho = min(SigmaTm[i,j],0.9999)
          mu1 = mean(X.test%*%theta[i,]) ; mu2 = mean(X.test%*%theta[j,])
          mat_dens[i,j]<<-sum(pmax(log(dbipoilog(n1=Ytest[,i],n2=Ytest[,j],mu1=mu1,mu2=mu2,
                                                 sig=sig,sig2=sig2,rho=rho)),-709))
          mat_dens[j,i]<<-mat_dens[i,j]
        })
      })
    }else{
      sapply(1:(p-1), function(i){
        sapply((i+1):p, function(j){
          sig=sqrt(D[i]) ; sig2=sqrt(D[j]) ; rho=SigmaTm[i,j]
          if(rho<0) rho = max(SigmaTm[i,j],-0.9999)
          if(rho>0) rho = min(SigmaTm[i,j],0.9999)
          mu1 = theta[i] ; mu2 = theta[j]
          mat_dens[i,j]<<-sum(pmax(log(dbipoilog(n1=Ytest[,i],n2=Ytest[,j],mu1=mu1,mu2=mu2,
                                                 sig=sig,sig2=sig2,rho=rho)),-709))
          mat_dens[j,i]<<-mat_dens[i,j]
        })
      })
    }
      
    
    VCp<-sum(mat_dens[upper.tri(mat_dens, diag=FALSE)])
  }else{
    VCp<-NA
  }
  return(VCp)
}

#cross_val rend le vecteur des vraisemblances pour chaque bloc et arbres échnatillonnés
cross_val<-function(counts,r,X=NULL, nblocs=10,B=30,alpha=0.1,eps=1e-3, cores=3,path){
  T0<-Sys.time()
  blocs.vecPCL=lapply(1:nblocs, function(i){
    cat(paste0("Bloc ",i,":"))
    file=paste0(path,r,"_bloc",i,".rds")
   # if(!file.exists(file)){
      VEMfit=split_fit(Y=counts,X = X,  v=1-(1/nblocs),r=r,alpha=alpha, bloc=i,eps=eps)
      saveRDS(VEMfit, file=file)
   # }  
    cat(paste0("\n   Estimating pairwise composite likelihood..."))
    T1<-Sys.time()
    VEMfit=readRDS(paste0(path,r,"_bloc",i,".rds"))
 
    prob=VEMfit$Pg ; beta=VEMfit$beta ; q=ncol(beta)
    prob[prob>1]=1
    if(!is.null(beta)){
      if(is.finite(SumTree(beta))){
        beta <- beta / SumTree(beta)^(1/(q-1))
      }else{
        beta=beta/100
        beta <- beta / SumTree(beta)^(1/(q-1))
      }
      collec.trees=mclapply(1:B, function(x){
        tree=rSpanTreeV1(beta=beta,prob=prob)$tree
      }, mc.cores=cores)
      collec.trees=collec.trees[which(do.call(rbind, lapply(collec.trees, length))>1)]
      filtre_trees=unique(collec.trees)
      #compute p(Y|T) once for the same T
      logpY=do.call(rbind,lapply(filtre_trees, function(tree){
        pY_condT(Ytest = VEMfit$Ytest,X.test=VEMfit$X.test,tree=tree,M = VEMfit$M,S = VEMfit$S,
                 Omega_hat = VEMfit$Omega_hat,D=VEMfit$D,theta = VEMfit$theta,
                 r = r,plot = FALSE)
      }))

      vec.nb.occur.tree<-do.call(rbind, lapply(filtre_trees, function(tree){
        nb.occur=0
        for(i in 1:length(collec.trees)){#are the two trees the same ?
          nb.occur=nb.occur+1*(sum(tree!=collec.trees[[i]])==0)
        }
        return(nb.occur)
      }))
      vec_logpY = rep(logpY, vec.nb.occur.tree)
    }else{vec_logpY=NA}
    T2<-Sys.time()
    runtime=difftime(T2, T1)
    cat(paste0(round(runtime,3), attr(runtime, "units"),"\n"))
    return(vec_logpY)
  })
  Tn<-Sys.time()
  runtime=difftime(Tn, T0)
  cat(paste0(" in ",round(runtime,3), attr(runtime, "units"),"\n"))
  return(do.call(rbind,blocs.vecPCL))
}
ICL<-function(TrueJ, Pg,Wg ,S, r,d, omega){
  q=ncol(S) ; p=q-r ; n=nrow(S)
  H=(p+1):q
  if(r!=0){
    pen_ZH= 0.5*sum(log(S[,H])) + n*r*0.5*(1+log(2*pi))
  }else{
    pen_ZH= 0
  }
  pen_T=-( sum( Pg * log(Wg+(Wg==0)) ) - logSumTree(Wg)$det) 
  pen_r<-p*(d) + (p*(p+1)/2 +r*p+r)+(q*(q-1)/2 - 1) #d comprends l'intercept
  # norm=n*q*(q-1)/2
  ICL=TrueJ   - (pen_T + pen_ZH + pen_r*log(n)/2)
  return(ICL)
}
VBIC<-function(TrueJ,p,r,d,n){
  q=p+r
  nbparam<-p*(d) + (p*(p+1)/2 +r*p)+(q*(q-1)/2 - 1) #d comprends l'intercept
  vBIC=TrueJ - nbparam*log(n)/2
  return(vBIC)
}
#-- data simulation
# regarder qualité simule des 3 premiers comparé aux autres (et temps)
# n passé à 200 et bootstrap sur 50 au lieu de 200
n=200 ;p=14;type="scale-free" 
sapply(0, function(r){
  start1=ifelse(r==1, 1, 1)
  sapply(start1:3, function(seed){
    path="/Users/raphaellemomal/simulations/CV_example/simu_article/trueR"
    set.seed(seed)
    missing_data<-generate_missing_data(n,p,r,type,plot=FALSE)
    counts=missing_data$Y
    dir=paste0(path,r,"/VEMfit/seed",seed,"_VEMfit_r")
    start2=ifelse(r==0 && seed==10, 2, 0)
    sapply(start2:2, function(testR){
      cat("seed ",seed,", testing r=",testR,":\n")
      cv<-cross_val(counts,r=testR,cores = 3,path=dir) 
      saveRDS(cv, file=paste0(path,r,"/seed",seed,"_cv",testR,".rds"))
    })
  })
})
##Fatala avec covariables
library(ade4)
data(baran95)
counts=as.matrix(baran95$fau)
p=ncol(counts); n=nrow(counts)
set.seed(1)
reorder=sample(1:n, n,replace=FALSE)
counts=counts[reorder,]
X=baran95$plan
X=model.matrix(~X$date+X$site)
X=X[reorder,]
#devtools::load_all()

path="/Users/raphaellemomal/simulations/CV_example/Fatala/"
CV_fatala_0<-cross_val(counts,X = X,r=0,B=100,nblocs=10,cores=3,alpha=0.05,path=path)#1.1 h
saveRDS(CV_fatala_0, file=paste0(path,"CV_fatala_0_covar.rds"))
CV_fatala_1<-cross_val(counts,X = X,r=1,B=100,nblocs=10,alpha=0.05,cores=3,path=path)#5H
saveRDS(CV_fatala_1, file=paste0(path,"CV_fatala_1_covar.rds"))
CV_fatala_2<-cross_val(counts,X = X,r=2,B=100,nblocs=10,alpha=0.05,cores=3,path=path)#8.5h
saveRDS(CV_fatala_2, file=paste0(path,"CV_fatala_2_covar.rds"))


##Barents avec covariables
load("/Users/raphaellemomal/these/Data_SR/BarentsFish.Rdata")
counts=Data$count
p=ncol(counts); n=nrow(counts)
set.seed(1)
reorder=sample(1:n, n,replace=FALSE)
counts=counts[reorder,]
X=data.frame(Data$covariates)
X=model.matrix(~X$Temperature)
X=X[reorder,]
#devtools::load_all()
path="/Users/raphaellemomal/simulations/CV_example/Barents/"
test_covar<-cross_val(counts,X = X,r=0,B=1,
                nblocs=2,cores=3,alpha=0.05,path=path)#50min

test_null<-cross_val(counts,X = NULL,r=0,B=1,
                      nblocs=2,cores=3,alpha=0.05,path=path)#50min
library(PLNmodels)
PLNfit1<-PLN(counts~-1+X , control=list(trace=0))
PLNfit2<-PLN(counts~-1+.,data=data.frame(X) , control=list(trace=0))
PLNfit3<-PLN(counts~1, control=list(trace=0))
PLNfit3$loglik
PLNfit$var_par$
plot(PLNfit1$model_par$Theta,PLNfit3$model_par$Theta)
saveRDS(CV_fatala_0, file=paste0(path,"CV_barents_0_covar.rds"))
CV_fatala_1<-cross_val(counts,X = X,r=1,B=100,nblocs=10,alpha=0.05,cores=3,path=path)#3.3H
saveRDS(CV_fatala_1, file=paste0(path,"CV_barents_1_covar.rds"))
CV_fatala_2<-cross_val(counts,X = X,r=2,B=100,nblocs=10,alpha=0.05,cores=3,path=path)#2.7h
saveRDS(CV_fatala_2, file=paste0(path,"CV_barents_2_covar.rds"))

#---- figures
get_plot<-function(listmodels){
  res=data.frame(VCp=numeric(), r=numeric())
  rMax=length(listmodels)-1
  for(r in 0:rMax){
    data=listmodels[[r+1]]
    VCp=mean(data, na.rm=TRUE)
    res=rbind(res,data.frame(VCp, r=r))
  }
  data = res %>% as_tibble()  %>% group_by(r) %>% 
    summarize(mean.PCL=(VCp), sd.PCL= sd(VCp),
              inf=mean.PCL - 1.96*sd.PCL/sqrt(10),
              sup=mean.PCL + 1.96*sd.PCL/sqrt(10))
  g <-data %>% ggplot( aes(y=mean.PCL,x=r))+ geom_line()+
    geom_point(size=2.5)+theme_light()
  return(list(g , data))
}
path="/Users/raphaellemomal/simulations/CV_example/Barents/"
CV_barents_0=readRDS(paste0(path,"CV_barents_0_covar.rds"))
CV_barents_1=readRDS(paste0(path,"CV_barents_1_covar.rds"))
CV_barents_2=readRDS(paste0(path,"CV_barents_2_covar.rds"))
listmodels_B=list(CV_barents_0,CV_barents_1,CV_barents_2)
path="/Users/raphaellemomal/simulations/CV_example/Fatala/"
CV_fatala_0=readRDS(paste0(path,"CV_fatala_0_covar.rds"))
CV_fatala_1=readRDS(paste0(path,"CV_fatala_1_covar.rds"))
CV_fatala_2=readRDS(paste0(path,"CV_fatala_2_covar.rds"))
listmodels=list(CV_barents_0,CV_barents_1,CV_barents_2)
get_plot(listmodels_B)[[1]]+
  labs(title="Fatala selection model with pairwise composite likelihood")



# ---don't know what's this below
# Results
res=do.call(rbind,lapply(c(0,1), function(r){
  res=do.call(rbind,lapply(1:30, function(seed){
    path="/Users/raphaellemomal/simulations/CV_example/simu_article/trueR"
    res=do.call(rbind,lapply(0:2, function(testR){
      VCp= mean(readRDS(file=paste0(path,r,"/seed",seed,"_cv",testR,".rds")))
      return(data.frame(seed=seed,r=r, testR=testR,VCp=VCp))
    }))
    return(res)
  }))
  return(res)
}))

results=res %>% as_tibble() %>% group_by(seed,r) %>% summarise(resR=which.max(VCp)-1) 
table(results$r, results$resR)
g=results%>% mutate(joliR=ifelse(r==0,"r = 0","r = 1")) %>% 
  ggplot(aes(as.factor(resR)))+geom_bar(fill="#1B9E77") +facet_wrap(~joliR)+mytheme.dark("")+
  labs(x="Selected r", y="", title="Model selection results on 30 scale-free graphs")
ggsave(g, file="/Users/raphaellemomal/these/R/images/selecR.png")

results %>% filter(r==0, resR!=0)

#  calcul des autres critères de sélection
criteria=do.call(rbind,lapply(c(0,1), function(r){
  res= do.call(rbind, lapply(1:30, function(seed){
    cat("\nseed ",seed,": ")
    res=do.call(rbind, lapply(0:2, function(testR){
      cat("testR ",testR," ")
      res=do.call(rbind, lapply(1:10, function(bloc){
        VEMfit <- readRDS(paste0("/Users/raphaellemomal/simulations/CV_example/simu_article/trueR",r,"/VEMfit/seed",
                                 seed,"_VEMfit_r",testR,"_bloc",bloc,".rds"))
        if(!is.null(VEMfit$beta)){
          Pg=VEMfit$Pg ; W=VEMfit$beta ;Omega_hat=VEMfit$Omega_hat
          M=VEMfit$M ;S=VEMfit$S; d=1
          q=ncol(W); p=q-testR; O=1:p; H=(p+1):q; n=nrow(M)
          Pg[Pg>1]=1
          Rho=cov2cor((1/n)*(t(M)%*%M+diag(colSums(S))))
          Wg= computeWg(Rho, Omega_hat, W, testR, n, alpha=1/n,  hist=FALSE)
          # browser()
          # if(is.finite(SumTree(Wg))){
          #   Wg <- Wg / SumTree(Wg)^(1/(q-1))
          # }else{
          #   Wg=Wg/100
          #   Wg <- Wg / SumTree(Wg)^(1/(q-1))
          # }
          logSTW=logSumTree(W)$det ;  logSTWg=logSumTree(Wg)$det
          J=as.numeric(LowerBound(Pg, Omega_hat, M=M, S=S,W=W, Wg=Wg,p, logSTW,logSTWg)["J"])
          icl=ICL(J, Pg, Wg, S,r=testR,d)
          vbic=VBIC(J,p,testR,d,n)
          res=data.frame(seed=seed, r=r,testR=testR,bloc=bloc,icl=icl, J=J,vbic=vbic )
          return(res)
        }}))
      return(res)
    }))
    return(res)
  }))
  return(res)
})) 

resCrit=criteria %>% as_tibble() %>% group_by(r,seed, testR) %>%
  summarise(micl=mean(icl), mJ=mean(J),mvbic=mean(vbic)) %>% ungroup() %>% 
  group_by(r,seed) %>% 
  summarise(resICL=which.max(micl)-1,resJ=which.max(mJ)-1,resVBIC=which.max(mvbic)-1 )
finalRes=left_join(results, resCrit, by=c("r","seed"))

g=finalRes%>% mutate(joliR=ifelse(r==0,"r = 0","r = 1")) %>% rename(PCL=resR, ICL=resICL, `Low.B`=resJ, vBIC=resVBIC) %>% 
  gather(criteria, res, -seed, -r,-joliR) %>% 
  ggplot(aes(as.factor(res), color=fct_relevel(criteria,c("PCL","ICL","Low.B","vBIC")),
             fill=fct_relevel(criteria,c("PCL","ICL","Low.B","vBIC"))))+geom_bar(position=position_dodge()) +
  facet_wrap(~joliR)+mytheme.dark("")+
  labs(x="Selected r", y="", title="Model selection results on 30 scale-free graphs")
ggsave(g, file="/Users/raphaellemomal/these/R/images/selecR_allcrit.png")


finalRes %>% group_by(r) %>% summarise(err.PCL=sum(r!=resR)/n(),
                                       err.ICL=sum(r!=resICL)/n(),
                                       err.vBIC=sum(r!=resVBIC)/n(),
                                       err.J=sum(r!=resJ)/n())



