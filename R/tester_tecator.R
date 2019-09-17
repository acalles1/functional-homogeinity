rm(list=ls())

source("dd-otherdepths.R")
source("DD-plot-test.R")
source("homogeinity_flores2018.R")
source("generate_curves.R")

data(tecator)

hfat <- tecator$absorp.fdata[tecator$y$Fat>=20]
lfat <- tecator$absorp.fdata[tecator$y$Fat<20]

res_flores <- Tester_Flores(hfat, lfat, B=1000, stat=P4, depth.function = depth.FM)
res_dd <- Tester(hfat, lfat, B=1000, depth.function = depth.FM)
res_dd_rp <- Tester(hfat, lfat, B=1000, depth.function = depth.RP)
res_dd_fd1 <- Tester.od(hfat$data, lfat$data, B=1000)

print('The result with the Flores Test is')
print(res_flores)
print('The result with the DD-FM test is')
print(res_dd)
print('The result with the DD-RP test is')
print(res_dd_rp)
print('The result with the DD-FD1 test is')
print(res_dd_fd1)