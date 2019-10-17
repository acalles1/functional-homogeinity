source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_curves.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_hitchcock2007.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/dd-otherdepths.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/DD-plot-test.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/homogeinity_flores2018.R")

K <- 100    
d <- seq(0.05, 0.5, 0.05)
n_d <- length(d)
results <- numeric(n_d)

for (i in 1:n_d){
    print('va este d')
    print(d[i])
    res <- numeric(K)
    for (j in 1:K){ 
        S0 <- Generator()
        J <- fdata(S0)
        S1 <- Generator(delta=d[i])
        G <- fdata(S1)
        res[j] <- Tester_Flores(J, G, B=1000, depth.function=depth.FM)
    }
    results[i]  <- 1 - sum(res)/K
    print(results[i])
}
results
