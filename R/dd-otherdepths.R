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

#source("generate_curves.R")

Bootstrapper.od <- function(J, G, B=1000, depth.function=depthf.fd1){
  # Computes the bootstrap statistics parallely.
  #
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # B: The amount of bootstrap statistics we're going to use.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  # nc: Number of clusters. In my PC I have 4 so this is the default value I
  # will use.
  #
  # Returns:
  # J.fun and G.fun: A list containing this two structures, which are
  # the resampled versions of J and G.
  #
  # Note: Not currently using beta0 for anything, but some versions of this code
  # used it. If you want to experiment with beta1 and beta0 simulataneously
  # its quite easy as the t-values for beta0 are computed here too.
  #
  #cl <- makeCluster(nc)
  Hh <- rbind(J, G)
  Hf <- rawfd2dataf(Hh, c(0,1))
  Jf <- rawfd2dataf(J, c(0,1))
  Gf <- rawfd2dataf(G, c(0,1))
  depths.inJ <- depth.function(Hf, Jf, range=c(0,1), d=30, order = 2, approx=150)
  depths.inG <- depth.function(Hf, Gf, range=c(0,1), d=30, order = 2, approx=150)
  depy <- depths.inJ$Half_FD
  depx <- depths.inG$Half_FD
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  beta0 <- numeric(B)
  beta1 <- numeric(B)
  t0 <- numeric(B)
  t1 <- numeric(B)
  kH <- dim(Hh)[1]
  kN <- dim(J)[1]
  kM <- dim(G)[1]
  for (i in 1:B){
    #print(i)
    #Resample and compute the statistic.
    Hs <- Hh[sample(1:kH,size=kH,replace=TRUE),]
    Hsf <- rawfd2dataf(Hs, c(0,1))
    Js <- Hs[1:kN, ]
    aux <- kN + 1
    #Next kH functions in H are the functions in the resampled version of H.
    Gs <- Hs[aux:kH, ]
    Jfs <- rawfd2dataf(Js, c(0,1))
    Gfs <- rawfd2dataf(Gs, c(0,1))
    #Resampled depth, s
    depths.inJs <- depth.function(Hsf, Jfs, range=c(0,1), d=30, order = 2, approx=150)
    depths.inGs <- depth.function(Hsf, Gfs, range=c(0,1), d=30, order = 2, approx=150)
    depys <- depths.inJs$Half_FD
    depxs <- depths.inGs$Half_FD
    lm.bs <- lm(depys ~ depxs)
    # Betas
    beta1[i] <- lm.bs$coefficients[2]
    beta0[i] <- lm.bs$coefficients[1]
    sum_aux <- sum((depxs - mean(depxs))^2)
    # Formula for the standard deviation of both betas. Standard formula from
    # regression analysis.
    beta1.std <- sqrt((sum(lm.bs$residuals^2))/((kH-2)*sum_aux))
    beta0.std <- sqrt((sum(lm.bs$residuals^2)*sum(depxs^2))/((kH-2)*sum_aux))
    t1[i] <- (beta1[i] - 1)/beta1.std
    t0[i] <- (beta0[i])/beta0.std
  }
  stats <- list()
  stats$b0 <- beta0
  stats$b1 <- beta1
  stats$t0 <- t0
  stats$t1 <- t1
  return(stats)
}

Tester.od <- function(J, G, B=1000, depth.function=depthf.fd1){
  # Tests homogeinity using bootstrap confidence intervals. H0: Homogeinity.
  #
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # B: The amount of bootstrap statistics we're going to use.
  # stat: The statistic we are going to use. You can choose P1-P4.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  # nc: Number of clusters. In my PC I have 4 so this is the default value I
  # will use.
  #
  # Returns:
  # TRUE if we do not reject the H0, meaning that J and G come from the same
  # population. False is we reject the H0, meaning that J and G are not
  # homogeneous and come from different popula 
  Hh <- rbind(J, G)
  kH <- dim(Hh)[1]
  Hf <- rawfd2dataf(Hh, c(0,1))
  Jf <- rawfd2dataf(J, c(0,1))
  Gf <- rawfd2dataf(G, c(0,1))
  depths.inJ <- depth.function(Hf, Jf, range=c(0,1), d=30, order = 2, approx=150)
  depths.inG <- depth.function(Hf, Gf, range=c(0,1), d=30, order = 2, approx=150)
  depy <- depths.inJ$Half_FD
  depx <- depths.inG$Half_FD
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  stats <- Bootstrapper.od(J, G, B=B, depth.function=depth.function)
  # Extract all statistics from the bootstrapper.
  b0 <- stats$b0
  b1 <- stats$b1
  t0 <- stats$t0
  t1 <- stats$t1
  sum_aux <- sum((depx - mean(depx))^2)
  beta1.std <- sqrt((sum(lm.ori$residuals^2))/((kH-2)*sum_aux))
  T.Test1 <- (b1.ori - 1)/beta1.std
  dist.t1 <- ecdf(t1)
  neg.T1 <- min(T.Test1, -T.Test1)
  p1.half.neg <- dist.t1(neg.T1)
  p1.half.pos <- 1.0 - dist.t1(-neg.T1)
  p1 <- p1.half.neg + p1.half.neg
  beta0.std <- sqrt((sum(lm.ori$residuals^2)*sum(depx^2))/((kH-2)*sum_aux))
  T.Test0 <- b0.ori/beta0.std
  dist.t0 <- ecdf(t0)
  neg.T0 <- min(T.Test0, -T.Test0)
  p0.half.neg <- dist.t0(neg.T0)
  p0.half.pos <- 1.0 - dist.t0(-neg.T0)
  p0 <- p0.half.neg + p0.half.pos
  p <- min(p0, p1)
  # 5% cutoff (2.5% upper and lower).
  t0crit <- quantile(t0, probs=c(0.025, 0.975))
  t1crit <- quantile(t1, probs=c(0.025, 0.975))
  beta0.l <- b0.ori - t0crit[2]*beta0.std
  beta0.u <- b0.ori - t0crit[1]*beta0.std
  beta1.l <- b1.ori - t1crit[2]*beta1.std
  beta1.u <- b1.ori - t1crit[1]*beta1.std
  res <- t0crit[1] <= T.Test0 && t0crit[2] >= T.Test0 && t1crit[1] <= T.Test1 && t1crit[2] >= T.Test1
  returns <- list()
  returns$res <- res
  returns$p <- p
  returns$t0crit <- t0crit
  returns$t1crit <- t1crit
  returns$T1 <- T.Test1
  returns$T0 <- T.Test0
  return(returns)
}
