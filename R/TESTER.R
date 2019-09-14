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
source("/home/cauchy/acalles/functional-homogeinity-master/R/dd-otherdepths.R")
K <- 100

#res_mat <- matrix(0, nrow=6, ncol=6)
#for (i in 1:6){
#   for (j in i:6){
#      print('model0')
#      print(i)
#      print('model1')
#      print(j)
#      trials <- seq(1:K)
#      fx1 <- function(trial){
#         S <- GenerateCurves()
#         J <- S[[i]]
#         S <- GenerateCurves()
#         G <- S[[j]]
#         resu <- Tester.od(J, G, B=1000)
#      }
#      tests1 <- mclapply(trials, fx1, mc.cores=8)
#      suma1 <- 0
#      for (k in 1:K){
#        if (tests1[[k]]$res){
#           suma1 <- suma1 + 1
#        }
#      }
#      res_mat[i,j] <- 1 - suma1/K
#      print(res_mat[i,j])
#   }
#}

#res_mat

#for (i in 1:n_deltas){
#      print('delta')
#      print(deltas[i])
#      trials <- seq(1:K)
#      fx1 <- function(trial){
#         S0 <- Generator()
#         J <- S0
#         S1 <- Generator(delta=deltas[i])
#         G <- S1
#         resu <- Tester(J, G, B=1000)
#         res <- resu
#      }
#      tests1 <- mclapply(trials, fx1, mc.cores=8)
#      suma1 <- 0
#     pnow <- 0
#      for (k in 1:K){
#        if (tests1[[k]]$res){
#           suma1 <- suma1 + 1
#        }
#        pnow <- pnow + tests1[[k]]$p
#      }
#      print('res')
#      results[i] <- 1 - suma1/K
#      ps[i] <- pnow/K
#      print(results[i])
#      print('p-val')
#      print(ps[i])
#}

#results
