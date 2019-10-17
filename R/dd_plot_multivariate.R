rm(list=ls())  
library(depth)
library(ggplot2)
library(MASS)
library(tibble)
library(ggthemes)

theme_set(theme_tufte())
kNc <- 500
kNs <- 2
cov.e <- matrix(c(1, 0, 0, 1), nrow=kNs, ncol=kNs)
cov.h <- matrix(c(1.5, 0.3, 0.3, 1.5))
# Initializ covariance matrix for h(t) with zeroes.
mu1 <- numeric(kNs)
mu2 <- c(0.5, 1.3)
s1 <- mvrnorm(n = kNc, mu=mu1, Sigma=cov.e)
s2 <- mvrnorm(n = kNc, mu=mu2, Sigma=cov.e)
s <- rbind(s1, s2)

depx <- numeric(1000)
depy <- numeric(1000)

for (i in 1:1000){
  depx[i] <- depth(s[i,], s1, method='Liu')
  depy[i] <- depth(s[i,], s2, method='Liu')
}

sam <- tibble(G_depth = depx,
              F_Depth = depy)

ggplot(sam, aes(G_depth, F_Depth)) + 
  geom_point() +
  xlab("") +
  ylab("") +
  geom_smooth(method = "lm", se=FALSE)