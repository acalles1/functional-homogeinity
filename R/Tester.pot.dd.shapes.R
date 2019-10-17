rm(list = ls())

#
# Author: Alejandro Calle Saldarriaga. 27-06-2019.
#
# This file implements an original homogeinity test for functional data, based on
# the ideas in [Liu et al., 1999].
# The function here are:
# Bootstrapper, which implements a resampling scheme for our statistic.
# Tester, which uses the boostrapped statistics to reject or not reject the test.


source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_curves.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_hitchcock2007.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/dd-otherdepths.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/DD-plot-test.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/homogeinity_flores2018.R")

K <- 100
pot.ddfm <- matrix(0, 4, 4)
pot.ddfd1 <- matrix(0, 4, 4)
for (i in 1:4){
    for (j in i:4){
        print('Modelo 0')
        print(i)
        print('Modelo 1')
        print(j)
        trials <- seq(1:K)
        fx2 <- function(trial){
            S0 <- GenerateCurves.Hitchcock(k=0.1, c=1)
            J <- fdata(S0[[i]])
            S1 <- GenerateCurves.Hitchcock(k=0.1, c=1)
            G <- fdata(S0[[j]])
            resu <- Tester(J, G, B=1000, depth.function=depth.FM)
            res <- resu
        }
        fx3 <- function(trial){
            S0 <- GenerateCurves.Hitchcock(k=0.1, c=1)
            J <- S0[[i]]
            S1 <- GenerateCurves.Hitchcock(k=0.1, c=1)
            G <- S1[[j]]
            resu <- Tester.od(J, G, B=1000, depth.function=depthf.fd1, range=c(0,5))
            res <- resu
        }
        tests2 <- mclapply(trials, fx2, mc.cores=8)
        tests3 <- mclapply(trials, fx3, mc.cores=8)
        suma1 <- 0
        suma2 <- 0
        for (k in 1:K){
            if (tests2[[k]]$res){
                suma1 <- suma1 + 1
            }
            if (tests3[[k]]$res){
            suma2 <- suma2 + 1
            }
        }
        print('resultados potencias')
        pot.ddfm[i,j] <- 1 - suma1/K
        pot.ddfd1[i,j] <- 1 - suma2/K
        print(pot.ddfm[i,j])
        print(pot.ddfd1[i,j])
    }
}

pot.ddfm
pot.ddfd1
