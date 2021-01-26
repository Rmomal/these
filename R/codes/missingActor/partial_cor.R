# estimation des corrélations partielles améliorée par la réestimation des covariances avec G
library(ggm)
library(EMtree)
library(nestor)
library(gridExtra)
library(xtable)
n=200 ; p=15 ; S=200
data=data_from_scratch("erdos",p=p, n=n,signed = TRUE,dens = 0.2)
G=1*(data$omega!=0)
diag(G)=0
ggimage(G)
PLN_Y = PLNmodels::PLN(data$data~1)
M=PLN_Y$var_par$M
S2=PLN_Y$var_par$S2
S=t(M)%*%M+diag(colSums(S2))
Prob=EMtree(PLN.Cor=PLN_Y,plot=TRUE)$edges_prob
hist(log(Prob), breaks=30)
Ghat_p=1*(Prob>0.2)
auc(Prob, G)
min(Prob[G==1])
#--------------------
colnames(Ghat_p)=rownames(Ghat_p)=1:p
colnames(G)=rownames(G)=1:p
colnames(S)=rownames(S)=1:p
fitggm=fitConGraph(amat = Ghat_p,S = S,n = n)
fitggm_star=fitConGraph(amat = G,S = S,n = n)
grid.arrange(ggimage(fitggm$Shat),ggimage(S))
plot(S, fitggm$Shat)
corp=parcor(S = fitggm$Shat)
#hist(log(abs(corp)), breaks=40)
#corp[abs(corp)<2e-9]=0
corp_star=parcor(S = fitggm_star$Shat)
#hist(log(abs(corp_star)), breaks=40)
#corp_star[abs(corp_star)<2e-9]=0
corp_naive=parcor(S)
# hist(corp_naive, breaks=40)
# hist(corp_star, breaks=40)
# hist(corp, breaks=40)
# xtable(table(sign(F_Sym2Vec(corp)), -sign(F_Sym2Vec(data$omega))))
# xtable(table(sign(F_Sym2Vec(corp_star)), -sign(F_Sym2Vec(data$omega))))
# xtable(table(sign(F_Sym2Vec(corp_naive)), -sign(F_Sym2Vec(data$omega))))
# Hurray le corp naif qui utilise juste la sortie de PLN sans l'info du graph
# ne peut pas trouver les zéros, et il trie bien les non-nuls.
# Donc seuls les vraies arêtes sont bien estimées, mais on ne peut pas savoir 
# lesquelles c'est
plot(-cov2cor(solve(S)), corp_naive)
#------------------------
# pour faire joli
library(tidyverse)
library(ggridges)
tibble(noG=log(abs(F_Sym2Vec(corp_naive))),
       Ghat=log(abs(F_Sym2Vec(corp))),
       Gstar=log(abs(F_Sym2Vec(corp_star)))) %>%
  gather(key, value) %>%
  mutate(key=fct_recode(key, S="noG", `S+G`="Gstar", `S+Ghat`="Ghat")) %>% 
           ggplot(aes(value, key)) %>%   + geom_density_ridges(
    stat = "binline", bins = 30, scale = 0.95,draw_baseline = FALSE,
    fill="steelblue",color="steelblue",alpha=0.6) +theme_light()+
 # geom_vline(xintercept = -20, col="red", linetype="dashed")+
  labs(x="Log absolute partial correlations", y="")

# upper_tri <- -cov2cor(data$omega)
# upper_tri[lower.tri(upper_tri, diag=TRUE)]=NA
# # Fondre la matrice de corrélation
# library(reshape2)
# melted_cormat <- melt(upper_tri, na.rm = TRUE)
# # Créer un ggheatmap
#  ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "darkcyan", high = "darkorange", mid = "gray90", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab",
#                        name="True partial\nCorrelation") +
#   theme_light()+labs(x="",y="")
# draw_network(G, layout="nicely",btw_rank = 1, nodes_label = 1:15)

# plot(F_Sym2Vec(-cov2cor(data$omega)), F_Sym2Vec(corp_naive))
# plot(F_Sym2Vec(-cov2cor(data$omega)), F_Sym2Vec(corp_star))
# plot(F_Sym2Vec(-cov2cor(data$omega)), F_Sym2Vec(corp))
tibble(truth=F_Sym2Vec(-cov2cor(data$omega)),
       cp_naive=F_Sym2Vec(corp_naive),
       cp_star=F_Sym2Vec(corp_star),
       cp_hat=F_Sym2Vec(corp)) %>% 
  gather(key, value, -truth) %>%  
  mutate(key=fct_recode(key, S="cp_naive", `S+G`="cp_star", `S+Ghat`="cp_hat")) %>% 
  mutate(key=fct_relevel(key, "S", "S+Ghat","S+G")) %>% 
  ggplot(aes(value, truth))+labs(x="Estimated partial correlation",
                                 y="True values")+
  geom_abline(color="grey")+geom_point()+facet_wrap(~key)+theme_light()+
  theme(strip.background=element_rect(fill="gray50",colour ="gray50"))

tibble(truth=F_Sym2Vec(-cov2cor(data$omega)),
       cp_naive=F_Sym2Vec(corp_naive),
       cp_star=F_Sym2Vec(corp_star),
       cp_hat=F_Sym2Vec(corp)) %>% 
  gather(key, value, -truth) %>%  
  mutate(key=fct_recode(key, S="cp_naive", `S+G`="cp_star", `S+Ghat`="cp_hat")) %>% 
  ggplot(aes(value, truth, color=key))+labs(x="Estimated partial correlation",
                                 y="True values")+
  geom_abline(color="grey")+geom_point(alpha=0.6)+theme_light()
