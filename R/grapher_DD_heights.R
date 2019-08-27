library(fda)
library(ggplot2)
library(cowplot)
library(fda.usc)

theme_set(theme_tufte())
depth.function <- depth.mode

male <- t(growth$hgtm)
female <- t(growth$hgtf)
fmale <- fdata(male)
ffemale <- fdata(female)
H <- c(fmale, ffemale)
depths.inJ <- depth.function(H, fdataori=fmale, trim = 0) #with trim = 0.25
depths.inG <- depth.function(H, fdataori=ffemale, trim = 0)
depy <- depths.inJ$dep
depx <- depths.inG$dep
deps <- tibble(G_depth = depx,
               F_Depth = depy)

ggplot(deps, aes(G_depth, F_Depth)) + 
  geom_point() +
  xlab("") +
  ylab("") +
  geom_smooth(method = "lm", se=FALSE)
