# tests séléection de modèle PLN sur données baran95
library(PLNmodels)
library(ade4)
data("baran95")

sample.fit=sample(1:nrow(Y),round(0.8*nrow(Y)))
sample.test=setdiff(1:nrow(Y), sample.fit)
baran=prepare_data(baran95$fau, baran95$plan)

full<-PLN(Abundance~1, data=baran, control=list(covariance="full"))
diagonal<-PLN(Abundance~1, data=baran,control=list(covariance="diagonal"))
spherical<-PLN(Abundance~1, data=baran,control=list(covariance="spherical"))

BICs<-c(full$BIC,diagonal$BIC,spherical$BIC)
ICLs<-c(full$ICL,diagonal$ICL,spherical$ICL)
logLs<-c(full$loglik,diagonal$loglik,diagonal$loglik)
models<-c("full", "diagonal", "spherical")


data_frame(BICs,ICLs, logLs,models) %>% gather(key, value, -models) %>% 
  ggplot(aes(models, value, color=key))+geom_point()+mytheme

ggimage<-function(data){
  melted_data <- melt(data)
  ggplot(melted_data, aes(x=Var1, y=Var2, fill=value)) + theme_light()+labs(x="",y="")+
    geom_tile() +guides(fill=FALSE)+ theme(plot.title = element_text(size=10, hjust=0.5))+ coord_fixed()
}

ggimage(spherical$model_par$Sigma)
ggimage(full$model_par$Sigma)
library(EMtree)

net.diag=EMtree(PLN.Cor = diagonal,PLN = TRUE)
net.full=EMtree(PLN.Cor = full,PLN = TRUE)
net.spher=EMtree(PLN.Cor = spherical,PLN = TRUE)
draw_network(1*(net.diag$edges_prob>0.05))
draw_network(1*(net.full$edges_prob>0.2))
draw_network(1*(net.spher$edges_prob>0.05))
