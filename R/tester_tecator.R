rm(list=ls())

source("DD-plot-test.R")
source("homogeinity_flores2018.R")
source("generate_curves.R")

data(tecator)
hfat <- tecator$absorp.fdata[tecator$y$Fat>=20]
lfat <- tecator$absorp.fdata[tecator$y$Fat<20]

res_flores <- Tester_Flores(hfat, lfat, B=1000, stat=P4, depth.function = depth.FM, nc=3)
res_dd <- Tester(hfat, lfat, B=1000, depth.function = depth.FM, nc = 3)
res_dd_hmodal <- Tester(hfat, lfat, B=1000, depth.function = depth.mode, nc = 3)

print('The result with Flores Test is')
print(res_flores)
print('The result with DD-FM test is')
print(res_dd)
print('The result with DD-hm test is')
print(res_dd_hmodal)