library(fda)
library(tidyfun)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(wesanderson)

theme_set(theme_tufte())

male <- t(growth$hgtm)
female <- t(growth$hgtf)
aux <-list(female, male)
Ss <- ldply(aux)
Sdata <-data.matrix(Ss)
cases <- numeric(93)
for (k in 1:93){
  if (k <= 54){
    cases[k] <- "Male"
  }else{
    cases[k] <- "Female"
    }
}

berkley_df = tibble(
  case = cases,
  f = tfd(Sdata, arg =growth$age))

ggplot(berkley_df) + 
    geom_spaghetti(aes(y = f, col = case)) +
    scale_color_manual(values=wes_palette(name="Darjeeling1", n=2)) +
    xlab("Age (years)") + ylab("Height (cm)") +
    ggtitle("Growth curves for boys and girls") +
    theme(plot.title = element_text(hjust = 0.5))

