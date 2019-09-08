rm(list=ls())

library(fda)
library(tidyfun)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(wesanderson)

data(poblenou)
theme_set(theme_tufte())

cond1 <- as.integer(poblenou$df[,"day.week"])>5
cond2 <- poblenou$df$day.festive==1
cond_wknd_festive <- cond1 | cond2
S0 <- poblenou$nox[cond_wknd_festive]
S1 <- poblenou$nox[!cond_wknd_festive]

cases <- numeric(115)
for (k in 1:115){
  if (cond_wknd_festive[k]){
    cases[k] <- "Festive or Weekend"
  }else{
    cases[k] <- "Week"
  }
}

nox_df = tibble(
  case = cases,
  f = tfd(poblenou$nox$data, arg =poblenou$nox$argvals))

ggplot(nox_df) + 
  geom_spaghetti(aes(y = f, col = case)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", n=2)) +
  xlab("Time (hours)") + ylab(expression(mglm^3)) +
  ggtitle("NOx levels") +
  theme(plot.title = element_text(hjust = 0.5))