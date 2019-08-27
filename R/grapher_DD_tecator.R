library(fda)
library(ggplot2)
library(cowplot)
library(fda.usc)

theme_set(theme_tufte())
data(tecator)
hfat <- tecator$absorp.fdata[tecator$y$Fat>=20]
lfat <- tecator$absorp.fdata[tecator$y$Fat<20]
depth.function <- depth.mode

depths.inJ <- depth.function(tecator$absorp.fdata, fdataori=hfat, trim = 0) #with trim = 0.25
depths.inG <- depth.function(tecator$absorp.fdata, fdataori=lfat, trim = 0)
depy <- depths.inJ$dep
depx <- depths.inG$dep
deps <- tibble(G_depth = depx,
               F_Depth = depy)

ggplot(deps, aes(G_depth, F_Depth)) + 
  geom_point() +
  xlab("") +
  ylab("") +
  geom_smooth(method = "lm", se=FALSE)
