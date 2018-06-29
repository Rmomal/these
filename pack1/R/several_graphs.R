library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)


grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}
auc_cov_no_estim<-readRDS("/home/momal/Git/these/pack1/R/Simu/erdos/prob/auc_1cov.rds")
  auc_covcov<-readRDS("/home/momal/Git/these/pack1/R/Simu/PLN/erdos/prob/auc_covcov.rds")
  auc_nocov<-readRDS("/home/momal/Git/these/pack1/R/Simu/PLN/erdos/prob/auc.rds")
  auc4_nocov_estim<-readRDS("/home/momal/Git/these/pack1/R/Simu/erdos/prob/auc.rds")

  plot_auc <- function(data,title) {
    tab <- data.frame(data)
    lignes <- which(is.na(tab[, 1]))
    if (length(lignes) != 0)
      tab <- tab[-lignes,]
    tab <- gather(tab,
                  key = method,
                  value = value,
                  treeggm,
                  ggm1step,
                  glasso)
    tab <- summarise(
      group_by(tab, var, method),
      mns = median(value),
      inf = quantile(value, 0.25),
      sup = quantile(value, 0.75)
    )
    param <- "prob"
    type <- "erdos"
    variable <- "AUC"

    ggplot(tab, aes(
      y = mns,
      x = as.numeric(as.character(var)),
      group = method,
      color = method
    )) +
      geom_errorbar(aes(ymin = inf, ymax = sup),
                    width = 0,
                    position = position_dodge((max(tab$var) - min(tab$var)) / 100)) +
      #geom_smooth(se=FALSE,size=0.3)+
      geom_point() +
      geom_line(size = 0.2) +
      # geom_linerange(aes(ymin = quantile(value,0.25), ymax = quantile(value,0.75)),group=tab$method)+
      labs(y = variable, x = param) +
      scale_color_manual(
        values = c("#076443", "#56B4E9", "#E69F00"),
        name = "Method:",
        breaks = c("treeggm", "ggm1step", "glasso"),
        labels = c("EM ", "1 step", "SpiecEasi")
      ) +
      scale_y_continuous(limits = c(0.5, 1)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank())+
      labs(title=title,x="Edge probability")
  }

  p1 <- plot_auc(  auc_cov_no_estim , "With unestimated covariable" )
  p2 <-plot_auc(auc_covcov  , "With covariates")
  p3 <- plot_auc(auc_nocov, "Without covariates" )
  p4 <-plot_auc(auc4_nocov_estim, "Estimated covariables unrelated to the data")
  grid_arrange_shared_legend(p3, p2, ncol = 2, nrow = 1)
