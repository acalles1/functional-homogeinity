library(fda)
library(ggplot2)
library(cowplot)
library(fda.usc)

theme_set(theme_tufte())
data(MCO)
control <- MCO$intact[MCO$classintact==1]
treatment <- MCO$intact[MCO$classintact==2]
depth.function <- depth.mode

depths.inJ <- depth.function(MCO$intact, fdataori=control, trim = 0) #with trim = 0.25
depths.inG <- depth.function(MCO$intact, fdataori=treatment, trim = 0)
depy <- depths.inJ$dep
depx <- depths.inG$dep
deps <- tibble(G_depth = depx,
               F_Depth = depy)

ggplot(deps, aes(G_depth, F_Depth)) + 
  geom_point() +
  xlab("") +
  ylab("") +
  geom_smooth(method = "lm", se=FALSE)
