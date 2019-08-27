library(fda.usc)
library(tidyfun)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(wesanderson)

theme_set(theme_tufte())

MCO_data <- MCO$intact$data
cases <- numeric(89)
for (k in 1:89){
  if (MCO$classintact[k] == 1){
    cases[k] <- "Control Group"
  }else{
    cases[k] <- "Treatment Group"
  }
}

mco_df = tibble(
  case = cases,
  f = tfd(MCO_data, arg = MCO$intact$argvals))

ggplot(mco_df) + 
  geom_spaghetti(aes(y = f, col = case)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", n=2)) +
  xlab("time (s)") + ylab("Ca") +
  ggtitle("Intact cardiac mice cells MCO levels") +
  theme(plot.title = element_text(hjust = 0.5))

