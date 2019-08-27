rm(list=ls())

library(fda)

source("DD-plot-test.R")
source("homogeinity_flores2018.R")
source("generate_curves.R")

male <- t(growth$hgtm)
female <- t(growth$hgtf)
fmale <- fdata(male)
ffemale <- fdata(female)

res_flores <- Tester_Flores(fmale, ffemale, B=1000, stat=P4, depth.function = depth.FM, nc=3)
res_dd <- Tester(fmale, ffemale, B=1000, depth.function = depth.FM, nc = 3)
res_dd_hmodal <- Tester(fmale, ffemale, B=1000, depth.function = depth.mode, nc = 3)

print('The result with Flores Test is')
print(res_flores)
print('The result with DD-FM test is')
print(res_dd)
print('The result with DD-hm test is')
print(res_dd_hmodal)