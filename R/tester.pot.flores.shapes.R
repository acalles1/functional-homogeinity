source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_curves.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/generate_hitchcock2007.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/dd-otherdepths.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/DD-plot-test.R")
source("/home/cauchy/acalles/functional-homogeinity-master/R/homogeinity_flores2018.R")

K <- 10
pot.flores <- matrix(0, 4, 4)

for (i in 1:4){
    for (j in i:4){
        res <- logical(K)
        print('Modelo 0')
        print(i)
        print('Modelo 1')
        print(j)
        for (k in 1:K){
            #print('Iteracion')
            #print(k)
            S0 <- GenerateCurves.Hitchcock(k=0.1, c=1)
            J <- fdata(S0[[1]])
            S1 <- GenerateCurves.Hitchcock(k=0.1, c=1)
            G <- fdata(S1[[2]])
            res[k] <- Tester_Flores(J, G, B=1000, depth.function=depth.FM)
        }
        pot.flores[i,j] <- 1 - sum(res)/K
        print(pot.flores[i,j])
    }
}

pot.flores
