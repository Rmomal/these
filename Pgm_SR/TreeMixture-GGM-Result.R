# Simulation results : comparison Tree / Glasso

rm(list=ls()); par(pch=20, mfrow=c(1, 1)); 

# Parms
p = 20; n = 100; 
Type.list = c('Erdos', 'Cluster', 'Tree')
Connex.list = c(T, F)
Algo.list = c('1step', 'EM', 'GL')
d.list = c(1, 2, 3, 4, 6, 8, 10)/p; # d.nb = length(d.list)
B = 20; 

# Results
Type.vec = Connex.vec = d.vec = Algo.vec = Perf.vec = Col.vec = c()
i = 0
for(Type in Type.list){
   for(Connex in Connex.list){
      for(d in d.list){
         SimDir = paste0('../Simu/', Type, '-connex', Connex, '-p', p, '-n', n, '/')
         SimResName = paste0(SimDir, 'TreeMixture-GGM-', Type, '-p', p, '-d_', round(100*d), '-n', n, '.Rdata')
         if(file.exists(SimResName)){
            print(SimResName); load(SimResName); Index = (i+1):(i+(3*B))
            Type.vec[Index] = rep(Type, B); Connex.vec[Index] = rep(Connex, B); d.vec[Index] = rep(d, B)
            Algo.vec[Index] = rep(Algo.list, each=B)
            Col.vec[Index] = rep(c(1, 2, 4), each=B)
            Perf.vec[Index] = c(Perf.1step, Perf.EM, Perf.GL)
            i = i + (3*B)
         }
      }
   }
}

# Plot
Algo.nb = length(Algo.list); Algo.col = c(1, 2, 4)
for (Type in Type.list){
   for (Connex in Connex.list){
      # Type = 'Erdos'; Connex = T
      print(paste(Type, Connex))
      Select = which((Type.vec==Type) & (Connex.vec==Connex))
      if (length(Select)>0){
         plot(d.vec[Select], rep(0, length(Select)), ylim=c(0, 1), col=0, ylab='AUC', xlab='d', 
              main=paste(Type, Connex))
         for (a in 1:Algo.nb){
            SubSel = which((Type.vec==Type) & (Connex.vec==Connex) & (Algo.vec==Algo.list[a]))
            # points(d.vec[SubSel]+(.0025*(a-(Algo.nb+1)/2)), Perf.vec[SubSel], col=Algo.col[a], pch='o')
            lines(sort(unique(d.vec[SubSel])), as.vector(by(Perf.vec[SubSel], d.vec[SubSel], mean)), col=Algo.col[a], lwd=2)
            lines(sort(unique(d.vec[SubSel])), as.vector(by(Perf.vec[SubSel], d.vec[SubSel], quantile, 0.9)), col=Algo.col[a], lty=2)
            lines(sort(unique(d.vec[SubSel])), as.vector(by(Perf.vec[SubSel], d.vec[SubSel], quantile, 0.1)), col=Algo.col[a], lty=2)
         }
      }
   }
}
