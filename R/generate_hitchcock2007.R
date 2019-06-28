library(pracma)
library(MASS)

GenerateCurves.Hitchcock <- function(sigma2=1, beta=1){
  # Generate the curves described in Hitchcock et al. [2007] for homogeinity tests.
  #
  # Args:
  # none
  # 
  # Returns:
  #
  # S, a list with each of the samples, where the first one is the reference 
  # sample (S0), the next 5 ones are the other samples, and the last one
  # is another realization of the reference sample.
  #
  # Number of observations for each curve.
  kNs <- 50
  # Number of curves.
  kNc <- 40
  # 30 equidistant points (ts = timesteps).
  ts <- linspace(0, 5, n = kNs) 
  # Initialize covariance matrix for first e(t) with zeroes.
  cov.e <- matrix(rep(0, len=kNs*kNs), nrow=kNs, ncol=kNs)
  # Initializ covariance matrix for h(t) with zeroes.
  cov.h <- matrix(rep(0, len=kNs*kNs), nrow=kNs, ncol=kNs) 
  mu1 <- numeric(kNs) 
  mu2 <- numeric(kNs)
  mu3 <- numeric(kNs)
  mu4 <- numeric(kNs)
  # Fill up covariances and expected values as described in the paper.
  for (i in 1:kNs) {
    mu1[i] <- -sin(ts[i]-1)*log(ts[i]+0.5)
    mu2[i] <- cos(ts[i])*log(ts[i]+0.5)
    mu3[i] <- -0.25 - 0.1*cos(0.5*(ts[i]-1))*ts[i]^(1.5)*sqrt(5*ts[i]^0.5 + 0.5)
    mu4[i] <- 0.6*cos(ts[i])*log(ts[i]+0.5)*sqrt(ts[i]+0.5)
    for (j in 1:kNs) {
      cov.e[i,j] <- sigma2*((2*beta)^-1)*exp(-beta*abs(ts[i] - ts[j]))
    } 
  }
  err1 <- rnorm(kNs, 0, sigma2)
  err2 <- rnorm(kNs, 0, sigma2)
  err3 <- rnorm(kNs, 0, sigma2)
  err4 <- rnorm(kNs, 0, sigma2)
  
  mu <- numeric(kNs)
  e.S1 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  e.S2 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  e.S3 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  e.S4 <- mvrnorm(n = kNc, mu=mu, Sigma=cov.e)
  
  S1.aux <- sweep(e.S1, 2, mu1, "+")
  S2.aux <- sweep(e.S2, 2, mu2, "+")
  S3.aux <- sweep(e.S3, 2, mu3, "+")
  S4.aux <- sweep(e.S4, 2, mu4, "+")
  
  S1 <- sweep(S1.aux, 2, err1, "+")
  S2 <- sweep(S2.aux, 2, err2, "+")
  S3 <- sweep(S3.aux, 2, err3, "+")
  S4 <- sweep(S4.aux, 2, err4, "+")

  S <- list(S1, S2, S3, S4)
  
  return(S)
  
}

# Debugging
#S <- GenerateCurves.Hitchcock()
#ts <- linspace(0, 5, n = 50) 
#S1 <- S[[1]]
#S2 <- S[[2]]
#S3 <- S[[3]]
#S4 <- S[[4]]
#plot(ts, S1[1,], type='l', col='black')
#for (i in 2:40){
#  lines(ts,S1[i,],col="black")
#}
#for (i in 1:40){
#  lines(ts, S2[i,], col="red")
#}
#for (i in 1:40){
#  lines(ts, S3[i,], col="green")
#}
#for (i in 1:40){
#  lines(ts, S4[i,], col="cyan")
#}

# Hitchcock, D. B., Booth, J. G. & Casella, G. (2007), ‘The effect of pre-smoothing
# functional data on cluster analysis’, Journal of Statistical Computation and 
# Simulation 77(12), 1043–1055