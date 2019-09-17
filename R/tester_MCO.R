rm(list=ls())

source("DD-plot-test.R")
source("homogeinity_flores2018.R")
source("generate_curves.R")
source("dd-otherdepths.R")

data(MCO)
control <- MCO$intact[MCO$classintact==1]
treatment <- MCO$intact[MCO$classintact==2]

res_flores <- Tester_Flores(control, treatment, B=1000, stat=P4, depth.function = depth.FM)
res_dd <- Tester(control, treatment, B=1000, depth.function = depth.FM)
res_dd_rp <- Tester(control, treatment, B=1000, depth.function = depth.RP)
res_dd_fd1 <- Tester.od(control$data, treatment$data, B=1000)

print('The result with Flores Test is')
print(res_flores)
print('The result with the DD-FM test is')
print(res_dd)
print('The result with the DD-RP test is')
print(res_dd_rp)
print('The result with the DD-fd1 test is')
print(res_dd_fd1)