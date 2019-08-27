source("generate_curves.R")

library(tidyfun)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(wesanderson)
library(cowplot)

theme_set(theme_tufte())

cont <- 0
plots <- vector("list", 6) 

for (i in 1:6){
    cont <- cont + 1
    message(cont)
    S <- GenerateCurves()
    S0 <- S[[1]]
    S <- GenerateCurves()
    S1 <- S[[i+1]]
    aux <-list(S0, S1)
    Ss <- ldply(aux)
    Sdata <-data.matrix(Ss)
    cases <- numeric(100)
    for (k in 1:100){
      if (k <= 50){
        cases[k] <- paste("Model", as.character(0))
      }else{
        if (i == 6){
          cases[k] <- paste("Model 0, other sample")
        } else{
        cases[k] <- paste("Model", as.character(i+1))
        }
      }
    }
    sim_df = tibble(
      case = cases,
      f = tfd(Sdata, arg = seq(0,1, l = 30)))
    
    plots[[cont]] <- local({
      cont <- cont
      p <- ggplot(sim_df) + 
        geom_spaghetti(aes(y = f, col = case, alpha = .3), show.legend = FALSE) +
        scale_color_manual(values=wes_palette(name="Darjeeling1", n=2)) +
        ylab("")
      print(p)
    })
}

plot_grid(plotlist=plots, nrow = 2, ncol = 3)

#theme_minimalist
#theme_tufte