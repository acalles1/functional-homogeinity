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
    t1[i] <- (beta1[i] - b1.ori)/beta1.std
    t0[i] <- (beta0[i] - b0.ori)/beta0.std
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
  # homogeneous and come from different populations.
  #
  H <- c(J, G)
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
  # 5% cutoff (2.5% upper and lower).
  cvalb0 <- quantile(t0, probs=c(0.025, 0.975))
  cvalb1 <- quantile(t1, probs=c(0.025, 0.975))
  beta0.l <- mean(b0) - cvalb0[2]*sqrt(var(b0))
  beta0.u <- mean(b0) - cvalb0[1]*sqrt(var(b0))
  beta1.l <- mean(b1) - cvalb1[2]*sqrt(var(b1))
  beta1.u <- mean(b1) - cvalb1[1]*sqrt(var(b1))
  return(beta1.l <= 1 && beta1.u >= 1)
}
