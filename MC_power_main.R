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

H = c("0", "1", "2", "3")
K = c("3", "5", "7")
a_names_temp = c()
b_names_temp = c()
for (s in c("1", "2")){
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
 
  a_names_var = c(a_names_var, paste0("P0var", H[j], s, "sided"))
  b_names_var = c(b_names_var, paste0("P1var", H[j], s, "sided"))
  b_names_var_min = c(b_names_var_min, paste0("P1var", H[j], "min", s, "sided"))
  for (k in 1:3){
    a_names_sel = c(a_names_sel, paste0("P0sel", H[j], K[k], s, "sided"))
    b_names_sel = c(b_names_sel, paste0("P1sel", H[j], K[k], s, "sided"))
    b_names_sel_min = c(b_names_sel_min, paste0("P1sel", H[j], K[k], "min", s, "sided"))
    if (k<3){
      a_names_iv = c(a_names_iv, paste0("P0iv", H[j], K[k], s, "sided"))
      b_names_iv = c(b_names_iv, paste0("P1iv", H[j], K[k], s, "sided"))
      b_names_iv_min = c(b_names_iv_min, paste0("P1iv", H[j], K[k], "min", s, "sided"))
    }
  }
}
a_names = c(a_names_sel, a_names_iv, a_names_var)
b_names = c(b_names_sel, b_names_iv, b_names_var)
b_names_min = c(b_names_sel_min, b_names_iv_min, b_names_var_min)
a_names_temp = c(a_names_temp, a_names, a_names)
b_names_temp = c(b_names_temp, b_names, b_names_min)
}
a_names = a_names_temp
b_names = b_names_temp

finalMatrix <- foreach(i=1:length(a_names), .combine=cbind) %dopar% {
  source("MC_Tests.R")
  source("MC_power.R")
  sided = strtoi(substr(a_names[i], nchar(a_names[i]) - 5, nchar(a_names[i]) - 5))
  #if (sided == 2){
  #a <- read.csv(file=paste0("DGPs/", a_names[i], ".csv"), header=FALSE, sep=",")
  a <- read.csv(file=paste0(paste0("DGPs", a_names[i]), ".csv"), header=FALSE, sep=",")
  a = a[,1]
  #b <- read.csv(file=paste0("DGPs/",b_names[i], ".csv"), header=FALSE, sep=",")
  b <- read.csv(file=paste0(paste0("DGPs",b_names[i]), ".csv"), header=FALSE, sep=",")
  b = b[,1]
  tempMatrix = MC(a, b, 5000, 0.15, 15, 5000,sided, i)
  tempMatrix
  #}
}
#stop cluster
stopCluster(cl)
toc()
write.csv(finalMatrix, file = "PowerCurves/RejectionRates_Apr5.csv")

# for (i in 1:length(a_names)){
#   print(i)
#   sided = strtoi(substr(a_names[i], nchar(a_names[i]) - 5, nchar(a_names[i]) - 5))
#   #a <- read.csv(file=paste0("DGPs/", a_names[i], ".csv"), header=FALSE, sep=",")
#   a <- read.csv(file=paste0(paste0("DGPs", a_names[i]), ".csv"), header=FALSE, sep=",")
#   a = a[,1]
#   #b <- read.csv(file=paste0("DGPs/",b_names[i], ".csv"), header=FALSE, sep=",")
#   b <- read.csv(file=paste0(paste0("DGPs",b_names[i]), ".csv"), header=FALSE, sep=",")
#   b = b[,1]
#   tempMatrix = MC(a, b, 5000, 0.15, 15, 1,sided, i)
# }
