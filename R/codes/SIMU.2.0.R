source("/Users/raphaellemomal/simulations/codes/Kit_pour_VEM_EM.R")
library(parallel)
library(tidyverse)
##############
# FUNCTIONS
##############
Make_graph<-function(valeur,variable,type,path,nbgraph, d=20, r=10){
  if( variable=="n"){
    graph<-generator_graph(graph=type,d=d,prob=5/d,dens=5/d, r=r)
    param<-generator_param(as.matrix(graph))

    file<-paste0(paste0(path,type,"/n/Sets_param/Graph",nbgraph,".rds"))
    saveRDS(param,file)
  }else{
    graph<-switch(variable, "d"=generator_graph(graph=type,d=valeur,prob=5/valeur,dens=5/valeur, r=r),
                  "prob"=generator_graph(graph=type,prob=valeur, r=r),
                  "r"=generator_graph(graph=type,r=valeur,dens=5/d, r=r),
                  "dens"=generator_graph(graph=type,dens=valeur, r=r))
    param<-generator_param(as.matrix(graph))
    file<-paste0(path,type,"/",variable,"/Sets_param/Graph",nbgraph,"_",valeur,".rds")
    saveRDS(param,file)
  }

}


Make_data<-function(type, variable, covar, path, cores){
  mclapply(parameters[[variable]], function(valeur) {
    lapply(1:100, function(nbgraph) {
      data<-data_from_stored_graphs(type, variable, nbgraph, valeur, covar=covar,path, fixe=FALSE)
      saveRDS(data,paste0(path,type,"/",variable,"/YDataSFn30/Y_",nbgraph,"_",valeur,".rds"))
    })
  }, mc.cores = cores)
}

################################


scores_gcoda_spiec_resid<-function(type, variable, method,cores=1,covar=fixeCovar){#method ou spiecResid

  obj<-lapply(parameters[[variable]], function(valeur) {
    print(paste0("method: ",method, " // type: ", type," // var: ", variable," // valeur :",valeur))
    if(variable=="n") covar<-readRDS(paste0(path, "fixeCovar_n",valeur,".rds"))
    mclapply(1:100, function(nbgraph) {
      Y<-readRDS(paste0(path,type,"/",variable,"/YData/Y_",nbgraph,"_",valeur,".rds"))[[1]]
      #imputation par la mediane, parfois NA, bizarre
      # if(sum(is.na(Y))!=0){
      #   index<- which(is.na(Y),arr.ind=TRUE)
      #   if(sum(is.na(Y))==1){
      #     Y[index]<-median(Y[,index[,2]], na.rm=TRUE)
      #   }else{
      #     Y[index]<-apply(Y[,index[,2]],2,function(x) median(x, na.rm=TRUE))
      #   }
      # }
      if(method=="gcodaResid"){

        T1<-Sys.time()

        out_gcodaResid<-gcoda(Y, counts=T, covar=covar)
        K.score <- Reduce("+",out_gcodaResid$path)
        scores<- K.score / max(K.score)

        T2<-Sys.time()
        time<-difftime(T2,T1)

      }else{
        T1<-Sys.time()
        # browser()
        U<-t(clr.matrix(Y,mar=1))
        m<- model.matrix(~X1+X2+X3,covar)
        model<-lm(U~m)
        inf<- spiec.easi(model$residuals, method="glasso",icov.select = FALSE, nlambda = 50, verbose = FALSE)
        K.score <- Reduce("+",inf$est$path)
        scores<- K.score / max(K.score)

        T2<-Sys.time()
        time<-difftime(T2,T1)
      }
      return(list(scores,time))
    },mc.cores=cores)
  })
  #save scores
  lapply(seq_along(parameters[[variable]]), function(valeur) {
    lapply(1:Bgraph, function(nbgraph) {
      save_path<-paste0(path,type,"/",variable,"/Scores/")
      val<-parameters[[variable]][valeur]
      saveRDS(obj[[valeur]][[nbgraph]][[1]],paste0(save_path,"/",method,"_graph",nbgraph,"_",val,".rds"))
    })
  })
  #save times
  names(obj)<-parameters[[variable]]
  tidyObj<-obj  %>%
    as_tibble() %>%
    gather(valeur_param,listes) %>%
    mutate(listes = map(listes, ~as_tibble(t(.x)))) %>%
    unnest()%>%
    dplyr::select(-V1) %>%
    rename(time=V2) %>%
    unnest()

  saveRDS(tidyObj,paste0(path,type,"/",variable,"/",method,"_features2.rds"))
}

scores_EM<-function(type, variable, covar=fixeCovar, path, step="FALSE", Bgraph, cond.alt,cores=1, maxIter){
  # inferences

  expected_length<-6
  obj<-lapply(parameters[[variable]], function(valeur) {
    if(variable=="n") covar<-readRDS(paste0(path, "fixeCovar_n",valeur,".rds"))
    cat("\ntype: ", type," // var: ", variable," // valeur :",valeur)
    mclapply(1:Bgraph, function(nbgraph) {
      from_stored_graphs(type, variable, nbgraph, valeur,covar=covar,step=step,cond.tol=1e-12,path=path, maxIter=maxIter)
    },mc.cores=cores)
  })
  ## Check if success
  lengths_inlist<-sapply(seq_along(parameters[[variable]]), function(i){
    mean(sapply(obj[[i]],function(x){
      length(x[[1]])
    }))
  })

  success<-mean(lengths_inlist)==expected_length
  if(!success){
    # reperer dans quelle valeur, puis dans quel graph
    vals_fail<-which(lengths_inlist!=expected_length)
    graphs_fail<-lapply(vals_fail,function(valeur){
      which(sapply(obj[[valeur]],function(x){
        length(x[[1]])
      })!=expected_length)
    })
    # redo the failing inferences, without parallelization for easier storage too
    lapply(parameters[[variable]][vals_fail], function(valeur) {
      if(variable=="n") covar<-readRDS(paste0(path, "fixeCovar_n",valeur,".rds"))
      index_val<-which(parameters[[variable]]==valeur)
      seq_graphs<-graphs_fail[[which(parameters[[variable]][vals_fail]==valeur)]]
      lapply(seq_graphs, function(nbgraph) {
        obj[[index_val]][[nbgraph]] <<-from_stored_graphs(type, variable, nbgraph, valeur,
                                                          covar=covar,step=step, cond.tol=cond.alt,path=path)
      })
    })
  }

  ## save scores
  lapply(seq_along(parameters[[variable]]), function(valeur) {
    lapply(1:Bgraph, function(nbgraph) {
      save_path<-paste0(path,type,"/",variable,"/Scores/")
      val<-parameters[[variable]][valeur]
      if(step=="FALSE"){
        #  saveRDS(obj[[valeur]][[nbgraph]][[1]]$P,paste0(save_path,"EMmarg_graph",nbgraph,"_",val,".rds"))
        saveRDS(obj[[valeur]][[nbgraph]][[1]]$ProbaCond,paste0(save_path,"EMCond_graph",nbgraph,"_",val,".rds"))
      }else{
        # saveRDS(obj[[valeur]][[nbgraph]][[1]]$P,paste0(save_path,"OneMarg_graph",nbgraph,"_",val,".rds"))
        saveRDS(obj[[valeur]][[nbgraph]][[1]]$ProbaCond,paste0(save_path,"OneCond_graph",nbgraph,"_",val,".rds"))
      }
    })
  })

  ## check logL increasing
  #  browser()
  if(step=="FALSE"){
    sapply(seq_along(parameters[[variable]]), function(i){
      fails=0
      failedgraphs=c()
      sapply(seq_along(obj[[i]]),function(x){
        logP<-obj[[i]][[x]][[1]]$logpY
        s<- sum(which(logP[-1] - logP[-obj[[i]][[x]][[1]]$maxIter] <0))
        if(s!=0){failedgraphs<<-c(failedgraphs,x)
        fails<<-fails+1}
      })
      if(fails!=0){
        good<<-FALSE
        cat("\nlogLik did not increase for parameter ",parameters[[variable]][i]," for graphs ",failedgraphs)
      }
    })
  }
  #if(!good) browser()


  ##save times alpha and maxIter
  names(obj)<-parameters[[variable]]

  tidyObj<-obj  %>%
    as_tibble() %>%
    gather(valeur_param,listes) %>%
    mutate(listes = map(listes, ~as_tibble(t(.x)))) %>%
    unnest() %>%
    mutate(alpha=map(V1,~.x$alpha),maxIter=map(V1,~.x$maxIter)) %>%
    dplyr::select(-V1) %>%
    unnest()
  if(step=="FALSE"){
    saveRDS(tidyObj,paste0(path,type,"/",variable,"/EM_features.rds"))
  }else{
    saveRDS(tidyObj,paste0(path,type,"/",variable,"/OneStep_features.rds"))
  }
}

#######
# RUN
#######
parameters<-list(d=seq(10,30,2),n=seq(20,120,10), prob=seq(1,5,0.5)/20,
                 r=seq(1,50,5), dens=seq(1,5,0.5)/20)
Bgraph<-100
path<-"/Users/raphaellemomal/simulations/Simu/PLN.2.0/"
fixeCovar<-readRDS(paste0(path, "fixeCovar.rds")) #cbind(rep(c(0,1),each=n/2),rnorm(n,0,0.5),round(runif(n)*5))

n=100
covar<-data.frame(X1=round(runif(n)*5),X2=rnorm(n,0,2),X3=as.character(round(runif(n)*5)))
saveRDS(covar,paste0(path, "covar_SF_n30.rds"))

T1<-Sys.time()
for(type in "scale-free"){
  cparam<-switch(type,"cluster"=c("n","d","dens"),"erdos"=c("n","d"), "scale-free"="d")
  for(variable in cparam){
    # sapply(parameters[[variable]],function(valeur){
    #   cat("variable: ", variable," / valeur: ",valeur,"\n")
    #mclapply(1:Bgraph,function(nbgraph){
    Make_data(type, variable, covar=covar, path, cores=1) # for sF
    # Make_graph(valeur,variable,type,path,nbgraph)
    #   },mc.cores=3)
    # })
  }
}

T2<-Sys.time()
difftime(T2,T1)

T1<-Sys.time() # attention LONG, idem
for(type in c("cluster","erdos","scale-free")){
  cparam<-switch(type,"scale-free"=c("d","n"),"erdos"=c("n","d"), "cluster"=c("d","n","r"))
  for( variable in cparam){
    scores_EM(type, variable, covar=fixeCovar, path, step="TRUE",Bgraph=Bgraph,cond.alt=1e-6, cores=3)
  }
}
T2<-Sys.time()
time_EMtree1<-difftime(T2,T1) # 2.6min
T1<-Sys.time()
for(type in c("cluster","erdos","scale-free")){
  cparam<-switch(type,"scale-free"=c("d","n"),"erdos"=c("n","d"), "cluster"=c("d","n","r"))
  for( variable in cparam){
    scores_EM(type, variable, covar=fixeCovar, path, step="FALSE",Bgraph=Bgraph,cond.alt=1e-6, cores=3, maxIter=5)
  }
}
T2<-Sys.time()
time_EMtree<-difftime(T2,T1)

T1<-Sys.time() # attention LONG, idem
for(type in c("scale-free")){
  cparam<-switch(type,"scale-free"=c("d"),"erdos"=c("n","d"))  #  cparam<-"n"
  for( variable in cparam){

    #  Make_data(type, variable, covar=fixeCovar, path, cores=1)
    scores_EM(type, variable, covar=covar, path, step="TRUE",Bgraph=Bgraph,cond.alt=1e-6, cores=3)

  }
}
T2<-Sys.time()
time_EM1<-difftime(T2,T1)


T1<-Sys.time() #n beaucoup trop long, tester avec des valeurs faciles
for(type in c("cluster")){
  cparam<-switch(type,"scale-free"=c("d"),"erdos"=c("n","d"), "cluster"="d")  #cparam<-"n"
  for( variable in cparam){
    scores_gcoda_spiec_resid(type, variable, method="gcodaResid",cores=3,covar=fixeCovar)
  }
}
T2<-Sys.time()
time_gcoda<-difftime(T2,T1)

T1<-Sys.time() #idem
for(type in c("scale-free")){
  cparam<-switch(type,"scale-free"=c("d"),"erdos"=c("n","d"))  #cparam<-"n"
  for( variable in cparam){
    scores_gcoda_spiec_resid(type, variable, method="spiecResid",cores=3,covar=covar)
  }
}
T2<-Sys.time()
time_spiec<-difftime(T2,T1)





clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

clr.matrix <- function(x.f, mar=2, ...) {
  apply(x.f, mar, clr, ...)
}








