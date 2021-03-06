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

#source("generate_curves.R")

Bootstrapper <- function(J, G, B=1000, depth.function=depth.FM){
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
  H <- c(J, G)
  kN <- length(J)
  kM <- length(G)
  kH <- kN + kM
  #cl <- makeCluster(nc)
  depths.inJ <- depth.function(H, fdataori=J, trim = 0) #with trim = 0.25
  depths.inG <- depth.function(H, fdataori=G, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  beta0 <- numeric(B)
  beta1 <- numeric(B)
  t0 <- numeric(B)
  t1 <- numeric(B)
  for (i in 1:B){
    #Resample and compute the statistic.
    Hs <- H[sample(1:kH,size=kH,replace=TRUE),]
    Js <- Hs[1:kN, ]
    aux <- kN + 1
    #Next kH functions in H are the functions in the resampled version of H.
    Gs <- Hs[aux:kH, ]
    #Resampled depths
    depths.inJs <- depth.function(Hs, fdataori=Js, trim = 0)
    depths.inGs <- depth.function(Hs, fdataori=Gs, trim = 0)
    depys <- depths.inJs$dep
    depxs <- depths.inGs$dep
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

Tester <- function(J, G, B=1000, depth.function=depth.FM){
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
  H <- c(J, G)
  kN <- length(J)
  kM <- length(G)
  kH <- kN + kM
  depths.inJ <- depth.function(H, fdataori=J, trim = 0) #with trim = 0.25
  depths.inG <- depth.function(H, fdataori=G, trim = 0)
  depy <- depths.inJ$dep
  depx <- depths.inG$dep
  lm.ori <- lm(depy ~ depx)
  b1.ori <- lm.ori$coefficients[2]
  b0.ori <- lm.ori$coefficients[1]
  stats <- Bootstrapper(J, G, B=B, depth.function=depth.function)
  # Extract all statistics from the bootstrapper.
  b0 <- stats$b0
  b1 <- stats$b1
  t0 <- stats$t0
  t1 <- stats$t1
  sum_aux <- sum((depx - mean(depx))^2)
  beta1.std <- sqrt((sum(lm.ori$residuals^2))/((kH-2)*sum_aux))
  T.Test1 <- (b1.ori - 1)/beta1.std
  dist.t1 <- ecdf(t1)
  p1.half1 <- dist.t1(T.Test1)
  p1.half2 <- 1.0 - p1.half1
  p1 <- 2*min(p1.half1, p1.half2)
  beta0.std <- sqrt((sum(lm.ori$residuals^2)*sum(depx^2))/((kH-2)*sum_aux))
  T.Test0 <- b0.ori/beta0.std
  dist.t0 <- ecdf(t0)
  p2.half1 <- dist.t0(T.Test0)
  p2.half2 <- 1 - p2.half1
  p2 <- 2*min(p2.half1, p2.half2)
  # ACAAAAA
  p <- c(p1, p2)
  po <- sort(p)
  # Using the holm method
  res <- po[1] >= 0.025 && po[2] >= 0.05
  p.new <- p.adjust(p, "holm")
  p.min <- min(p.new)
  returns <- list()
  returns$res <- res
  returns$p <- p.min
  #HASTA ACAAAAAA
  return(returns)
}
