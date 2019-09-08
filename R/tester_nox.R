rm(list=ls())

source("DD-plot-test.R")
source("homogeinity_flores2018.R")

data(poblenou)

cond1 <- as.integer(poblenou$df[,"day.week"])>5
cond2 <- poblenou$df$day.festive==1
cond_wknd_festive <- cond1 | cond2
S0 <- poblenou$nox[cond_wknd_festive]
S1 <- poblenou$nox[!cond_wknd_festive]

res1 <- Tester(S0, S1, B=1000, depth.function=depth.FM)
res2 <- Tester(S0, S1, B=1000, depth.function=depth.RP)
res3 <- Tester_Flores(S0, S1, B=1000)