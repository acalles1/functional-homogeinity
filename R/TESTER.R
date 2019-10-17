rm(list = ls())

#
# Author: Alejandro Calle Saldarriaga. 27-06-2019.
#
# This file implements an original homogeinity test for functional data, based on
# the ideas in [Liu et al., 1999].
# The function here are:
# Bootstrapper, which implements a resampling scheme for our statistic.
# Tester, which uses the boostrapped statistics to reject or not reject the test.

library(fda.usc)
library(parallel)
library(car)
library(pracma)
library(ddalpha)


source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_curves.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_hitchcock2007.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/dd-otherdepths.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/DD-plot-test.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/homogeinity_flores2018.R")

K <- 100
#d <- seq(0.05, 0.5, 0.05)
#n_d <- length(d)
results.ddrpd <- matrix(0, 6, 6)
#results.ddfd1 <- numeric(n_d)
#p.ddfm <- numeric(n_d)
#p.ddfd1 <- numeric(n_d)

for (i in 1:6){
    for (j in i:6){
        print('modelo 1')
        print(i)
        print('modelo 2')
        print(j)
        trials <- seq(1:K)
        fx2 <- function(trial){
            S <- GenerateCurves()
            S0 <- S[[i]]
            J <- fdata(S0)
            S <- GenerateCurves()
            S1 <- S[[j]]
            G <- fdata(S1)
            resu <- Tester(J, G, B=1000, depth.function=depth.RPD)
            res <- resu
        }
        fx3 <- function(trial){
            S0 <- Generator()
            S1 <- Generator(delta=d[i])
            resu <- Tester.od(S0, S1, B=10, depth.function=depthf.fd1)
            res <- resu
        }
        tests2 <- mclapply(trials, fx2, mc.cores=8)
        #tests3 <- mclapply(trials, fx3, mc.cores=8)
        suma1 <- 0
        #suma2 <- 0
       # ps1 <- 0
        #ps2 <- 0
        for (k in 1:K){
            if (tests2[[k]]$res){
                suma1 <- suma1 + 1
            }
        #   if (tests3[[k]]$res){
        #       suma2 <- suma2 + 1
        #  }
      #      ps1 <- ps1 + tests2[[k]]$p
        #  ps2 <- ps2 + tests3[[k]]$p
        }
        print('res')
        results.ddrpd[i, j] <- 1 - suma1/K
        #results.ddfd1[i] <- 1 - suma2/K
      #  p.ddfm[i] <- ps1/K
    #  p.ddfd1[i] <- ps2/K
        print(results.ddrpd[i, j])
    # print(results.ddfd1[i])
   }
}

#p.ddrpd
#p.ddfd1
results.ddrpd
#results.ddfd1
