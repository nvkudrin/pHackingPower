#####################################################################
# Monte Carlo simulation
# Paper: (When) can we detect p-hacking?
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
#####################################################################

rm(list = ls())

library("fdrtool")
library("pracma")
library(gdata)
library(spatstat)
#library(haven)
library("rddensity")
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
library(foreach)
library(doParallel)
tic()
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

H = c("0", "1", "2", "3", "4")
K = c("3", "5", "7")

a_names_sel = c()
b_names_sel = c()
a_names_iv = c()
b_names_iv = c()
a_names_var = c()
b_names_var = c()

b_names_sel_min = c()

b_names_iv_min = c()

b_names_var_min = c()
for (j in 1:4){
 
  a_names_var = c(a_names_var, paste0("P0var", H[j]))
  b_names_var = c(b_names_var, paste0("P1var", H[j]))
  b_names_var_min = c(b_names_var_min, paste0("P1var", H[j], "min"))
  for (k in 1:3){
    a_names_sel = c(a_names_sel, paste0("P0sel", H[j], K[k]))
    b_names_sel = c(b_names_sel, paste0("P1sel", H[j], K[k]))
    b_names_sel_min = c(b_names_sel_min, paste0("P1sel", H[j], K[k], "min"))
    if (k<3){
      a_names_iv = c(a_names_iv, paste0("P0iv", H[j], K[k]))
      b_names_iv = c(b_names_iv, paste0("P1iv", H[j], K[k]))
      b_names_iv_min = c(b_names_iv_min, paste0("P1iv", H[j], K[k], "min"))
    }
  }
}
a_names = c(a_names_sel, a_names_iv, a_names_var)
b_names = c(b_names_sel, b_names_iv, b_names_var)
b_names_min = c(b_names_sel_min, b_names_iv_min, b_names_var_min)

a_names = c(a_names, a_names)
b_names = c(b_names, b_names_min)

finalMatrix <- foreach(i=1:length(a_names), .combine=cbind) %dopar% {
  source("MC_Tests.R")
  source("MC_power.R")
  a <- read.csv(file=paste0("DGPs/", a_names[i], ".csv"), header=FALSE, sep=",")
  a = a[,1]
  b <- read.csv(file=paste0("DGPs/",b_names[i], ".csv"), header=FALSE, sep=",")
  b = b[,1]
  tempMatrix = MC(a, b, 5000, 0.15, 15, 1, i)
  
  tempMatrix
}
#stop cluster
stopCluster(cl)
toc()
write.csv(finalMatrix, file = "PowerCurves/RejectionRates_April6new.csv")


