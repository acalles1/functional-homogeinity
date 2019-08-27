source("generate_curves.R")

library(ggplot2)
library(cowplot)
library(fda.usc)

theme_set(theme_tufte())
depth.function <- depth.mode
plots <- vector("list", 6) 

for (i in 1:6){
  S <- GenerateCurves()
  S0 <- S[[1]]
  S <- GenerateCurves()
  S1 <- S[[i+1]]
  J <- fdata(S0)
  G <- fdata(S1)
  H <- c(J, G)
  depths.inJ <- depth.function(H, fdataori=J, trim = 0) #with trim = 0.25
  depths.inG <- depth.function(H, fdataori=G, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  deps <- tibble(G_depth = depx,
                 F_Depth = depy)
  plots[[i]] <- local({
    i <- i
    p <- ggplot(deps, aes(G_depth, F_Depth)) + 
      geom_point() +
      xlab("") +
      ylab("") +
      geom_smooth(method = "lm", se=FALSE)
  })
}

plot_grid(plotlist=plots, nrow = 2, ncol = 3)

#theme_minimalist
#theme_tufte