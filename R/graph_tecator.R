library(fda.usc)
library(tidyfun)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(wesanderson)

data(tecator)

theme_set(theme_tufte())

tecator_data <- tecator$absorp.fdata$data
cases <- numeric(215)
for (k in 1:215){
  if (tecator$y$Fat[k] <= 20){
    cases[k] <- "Low Fat"
  }else{
    cases[k] <- "High Fat"
  }
}

tecator_df = tibble(
  case = cases,
  f = tfd(tecator_data, arg = tecator$absorp.fdata$argvals))

ggplot(tecator_df) + 
  geom_spaghetti(aes(y = f, col = case)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", n=2)) +
  xlab("Wavelength (mm)") + ylab("Absorbances") +
  ggtitle("Spectrometric curves for meat samples") +
  theme(plot.title = element_text(hjust = 0.5))

