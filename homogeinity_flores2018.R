# 
# Author: Alejandro Calle Saldarriaga. 24-01-2019.
# 
# This file has many functions which implements the main procedures described
# Flores et al. [2018], which implements an homogeinity test for functional
# data and some metrics to measure its performance. The functions here are:
# P1-P4, which are the four statistics described in the paper; ReSample, which
# is the resampling procedure described in the paper; Bootstrapper, which
# computes a confidence interval for P1-P4 using bootstraping procedures; 
# Tester which tests if two samples come from the same population; and an
# script at the end that can be used to see the power and size of the test.
#
# TODO(acalles): depth.RT does not work in P1-P4. Need to raise an error and a
# stop when the user tries to use this depth.

source("generate_curves.R")
library(fda.usc)
library(parallel)

P1 <- function(J, G, depth.function=depth.FM){
  # Computes the statistic P1 defined in Flores et al. [2018],
  # 
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  #
  # Returns:
  # P1: The computed statistic. 
  #
  # Depths of G in G.
  depths.G <- depth.function(G, trim=0)
  # Using the paper's notation. This is the deepest function of G in G. In 
  # other words, it is G's median.
  DgG <- depths.G$median
  # Now we need to compute the depth of this function in J.
  djDgG <- depth.function(DgG, fdataori=J, trim=0)
  return(djDgG$dep)
}

P2 <- function(J, G, depth.function=depth.FM){
  # Computes the statistic P2 defined in Flores et al. [2018],
  # 
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  #
  # Returns:
  # P2: The computed statistic. 
  #
  # Just a normalization of P1.
  kP2 <- abs(P1(J, G, depth.function) - P1(J, J, depth.function))
  return(kP2)
}

P3 <- function(J, G, depth.function=depth.FM){
  # Computes the statistic P3 defined in Flores et al. [2018],
  # 
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001].
  # Note: depth.RT does not work. It might be some error in implementation.
  #
  # Returns:
  # P3: The computed statistic. 
  #
  # Depths of J in G.
  depths.JinG <- depth.function(J, fdataori=G, trim=0)
  # Using the paper's notation. This is the deepest function of J in G.
  DjG <- depths.JinG$median
  # Now we need to compute the depth of this function in J.
  djDjG <- depth.function(DjG, fdataori=J, trim=0)
  return(djDjG$dep)
}

P4 <- function(J, G, depth.function=depth.FM){
  # Computes the statistic P4 defined in Flores et al. [2018],
  # 
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  # depth.function: Pick any functional depth you prefer. The default depth
  # is the Fraiman-Muniz depth, defined in Fraiman and Muniz [2001]. 
  # Note: depth.RT does not work. It might be some error in implementation.
  #
  # Returns:
  # P4: The computed statistic. 
  #
  # A normalization involving P1 and P3. Auxiliares to prevent excesive line
  # length.
  aux1 <- abs(P3(J, G, depth.function) - P1(J, J, depth.function))
  aux2 <- abs(P3(J, G, depth.function) - P1(G, G, depth.function))
  kP4 <- aux1*aux2
  return(kP4)
}

ReSample <- function(J, G){
  # Resamples the functions in the union H = J U G.
  # 
  # Args:
  # J: First sample, normally the reference sample: Sample 0.
  # G: Second sample, noramlly one of the other samples: Sample 1 to 5.
  #
  # Returns:
  # J.fun and G.fun: A list containing this two structures, which are
  # the resampled versions of J and G.
  #
  # Number of functions in J.
  kN <- length(J)
  # Number of points observed for each function in J.
  kP <- length(J[1]$data)
  # Number of functions in G.
  kM <- length(G)
  # Number of functions in the combined sample.
  kH <- kN + kM
  # H = J U G
  H <- c(J, G)
  # n + m random uniforms for resampling.
  r <- (kH)*runif(kH)
  # Initialize the matrix when we will store the resampled functions.
  new.H <-matrix(rep(0, len=kH*kP), nrow=kH, ncol=kP)
  for (i in 1:kH){
    # R to int to choose function in this position
    kr <- ceil(r[i])
    new.H[i, ] <- as.vector(H[kr]$data)
  }
  #First kN functions in H are the functions in the resampled version of J.
  new.J <- new.H[1:kN, ]
  aux <- kN + 1
  #Next kH functions in H are the functions in the resampled version of H. 
  new.G <- new.H[aux:kH, ]
  J.fun <- fdata(new.J)
  G.fun <- fdata(new.G)
  S <- list(J.fun, G.fun)
  return(S)
}

Bootstrapper <- function(J, G, B=1000, stat=P4, 
                         depth.function=depth.FM, nc = 4){
  # Computes the bootstrap statistics parallely.
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
  # J.fun and G.fun: A list containing this two structures, which are
  # the resampled versions of J and G.
  #
  # Number of functions in J.
  cl <- makeCluster(nc)
  bs <- numeric(B)
  #Need to export the global enviroment to the cluster environment
  clusterExport(cl=cl, list("ReSample","J","G","P4","P3",
                            "P2","P1", "depth.function"), envir=environment())
  #Need to export the libraries used to the cluster.
  clusterEvalQ(cl,library(pracma))
  clusterEvalQ(cl,library(fda.usc))
  #parSapply to obtain a vector.
  bs <- parSapply(cl, seq_len(B), function(i){
    #Resample and compute the statistic. 
    Ss <- ReSample(J, G)
    Js <- Ss[[1]]
    Gs <- Ss[[2]]
    x <- stat(Js, Gs, depth.function)
    }
  )
  stopCluster(cl)
  return(bs)
}

Tester <- function(J, G, B=1000, stat=P4, depth.function=depth.FM, nc=4){
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
  # True P4 between J and G.
  trueP4 <- P4(J, G)
  # Compute the bootstrap interval
  bs <- Bootstrapper(J, G, B, stat, depth.function, nc)
  # Sort and trim (upper because of the nature of the computations) to see if
  # the true P4 is in the constructed interval.
  sorted <- sort(bs)
  trimmed <- sorted[1:ceil(B*0.95)]
  lb <- min(trimmed)
  ub <- max(trimmed)
  if (lb <= trueP4 && trueP4 <= ub){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

# Main script. This just test 100 times if two simulated samples
# are homogeneous or not. Assing S0 to S[[1]] always, and change S1 to S[[2]]
# to S[[6]] to test heterogeinity, meaning that the values in tests must be
# FALSE. The power of the test will be (100 - sum(tests))/100
# Use S[[7]] to test homogeinity, meaning that the values in tests must
# be TRUE. Then the size of the test will be sum(tests)/100
tests <- logical(length=100)
ptm <- proc.time()
for (i in 1:100){
  # This routine takes a lot of time (aproximately half an hour in a 4-core PC)
  # so we will see what iteration the algorithm is running to have an
  # approximate measure of how much time the algorithm will run.
  print('iteration')
  print(i)
  # Generate samples.
  S <- GenerateCurves()
  # Choose samples for comparison.
  S0 <- S[[1]]
  S1 <- S[[6]]
  J <- fdata(S0)
  G <- fdata(S1)
  # Test.
  aux <- Tester(J, G)
  print(aux)
  tests[i] <- aux
}
proc.time() - ptm


# References
#
# Flores, Ramón, Lillo, Rosa and Romo, Juan. 2018. Homogeinity test for
# functional data. Journal of Applied Statistics, 45(5), 868-883.
#
# Fraiman, Ricardo and Muniz, Graciela. 2001. Trimmed means for functional
# data. Test, 10(2), 419-440.
