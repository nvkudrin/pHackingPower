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

#H = c("0", "1", "2", "3", "4")
H = c("0", "1", "2", "3") #Pi_RCT oly
K = c("3", "5", "7")
a_names_temp = c()
b_names_temp = c()
for (s in c("1", "2")){
  #for (s in c("2")){
a_names_sel = c()
b_names_sel = c()
a_names_iv = c()
b_names_iv = c()
a_names_iv_F = c()
b_names_iv_F = c()
a_names_var = c()
b_names_var = c()
a_names_clust = c()
b_names_clust = c()

b_names_sel_min = c()
b_names_iv_min = c()
b_names_iv_F_min = c()
b_names_var_min = c()
b_names_clust_min = c()

names_for_scattters_sel = c()
names_for_scattters_sel_min = c()
names_for_scattters_iv = c()
names_for_scattters_iv_min = c()

for (j in 1:4){
  for (g2s in 0:1){
    if ((s == "2")&(g2s==1)){
  a_names_var = c(a_names_var, paste0("P0var", H[j], s, "sided",g2s))
  b_names_var = c(b_names_var, paste0("P1var", H[j], s, "sided",g2s))
  b_names_var_min = c(b_names_var_min, paste0("P1var", H[j], "min", s, "sided",g2s))
  
  
    a_names_clust = c(a_names_clust, paste0("P0clust", H[j], s, "sided",g2s))
    b_names_clust = c(b_names_clust, paste0("P1clust", H[j], s, "sided",g2s))
    b_names_clust_min = c(b_names_clust_min, paste0("P1clust", H[j], "min", s, "sided",g2s))
    }
  for (k in 1:3){
    
    if (((g2s==1)&(s == "2"))|((g2s==0)&(k==1)&(s=="2"))|((g2s==1)&(k==1)&(s=="1"))){
    a_names_sel = c(a_names_sel, paste0("P0sel", H[j], K[k], s, "sided",g2s))
    b_names_sel = c(b_names_sel, paste0("P1sel", H[j], K[k], s, "sided",g2s))
    b_names_sel_min = c(b_names_sel_min, paste0("P1sel", H[j], K[k], "min", s, "sided",g2s))
    }
    if ((g2s==1)&(s == "2")){
    names_for_scattters_sel = c(names_for_scattters_sel, paste0("P1sel", H[j], K[k], s, "sided",g2s))
    names_for_scattters_sel_min = c(names_for_scattters_sel_min, paste0("P1sel", H[j], K[k], "min", s, "sided",g2s))
    }
    if ((s == "2")&(g2s==1)){
    if (k<3){
      a_names_iv = c(a_names_iv, paste0("P0iv", H[j], K[k], s, "sided",g2s))
      b_names_iv = c(b_names_iv, paste0("P1iv", H[j], K[k], s, "sided",g2s))
      b_names_iv_min = c(b_names_iv_min, paste0("P1iv", H[j], K[k], "min", s, "sided",g2s))
      
      a_names_iv_F = c(a_names_iv_F, paste0("P0iv_F", H[j], K[k], s, "sided",g2s))
      b_names_iv_F = c(b_names_iv_F, paste0("P1iv_F", H[j], K[k], s, "sided",g2s))
      b_names_iv_F_min = c(b_names_iv_F_min, paste0("P1iv_F", H[j], K[k], "min", s, "sided",g2s))
      
      names_for_scattters_iv = c(names_for_scattters_iv, paste0("P1iv", H[j], K[k], s, "sided",g2s))
      names_for_scattters_iv_min = c(names_for_scattters_iv_min,paste0("P1iv", H[j], K[k], "min", s, "sided",g2s))
    }
    }
  }
  }
}
a_names = c(a_names_sel, a_names_iv,a_names_iv_F, a_names_var, a_names_clust)
b_names = c(b_names_sel, b_names_iv,b_names_iv_F, b_names_var,b_names_clust)
b_names_min = c(b_names_sel_min, b_names_iv_min,b_names_iv_F_min, b_names_var_min,b_names_clust_min)
a_names_temp = c(a_names_temp, a_names, a_names)
b_names_temp = c(b_names_temp, b_names, b_names_min)
}
a_names = a_names_temp
b_names = b_names_temp

finalMatrix <- foreach(i=1:length(a_names), .combine=cbind) %dopar% {
  source("MC_Tests.R")
  source("MC_power.R")
  sided = strtoi(substr(a_names[i], nchar(a_names[i]) - 5-1, nchar(a_names[i]) - 5-1))
  #if (sided == 2){
  #a <- read.csv(file=paste0("DGPs/", a_names[i], ".csv"), header=FALSE, sep=",")
  a <- read.csv(file=paste0(paste0("DGPs_csv/selected/DGPs", a_names[i]), ".csv"), header=FALSE, sep=",")
  a = a[,1]
  #b <- read.csv(file=paste0("DGPs/",b_names[i], ".csv"), header=FALSE, sep=",")
  b <- read.csv(file=paste0(paste0("DGPs_csv/selected/DGPs",b_names[i]), ".csv"), header=FALSE, sep=",")
  b = b[,1]
  tempMatrix = MC(a, b, 5000, 0.15, 15, 5000,sided, i)
  tempMatrix
  #}
}
#stop cluster
stopCluster(cl)
toc()
write.csv(finalMatrix, file = "PowerCurves/RejectionRates_March18.csv")

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

m0 = c("P1sel032sided1", "P1sel03min2sided1", "P1iv032sided1", "P1iv03min2sided1", "P1var02sided1", "P1var0min2sided1", "P1clust02sided1", "P1clust0min2sided1")
m1 = c("P1sel132sided1", "P1sel13min2sided1", "P1iv132sided1", "P1iv13min2sided1", "P1var12sided1", "P1var1min2sided1", "P1clust12sided1", "P1clust1min2sided1")
m2 = c("P1sel232sided1", "P1sel23min2sided1", "P1iv232sided1", "P1iv23min2sided1", "P1var22sided1", "P1var2min2sided1", "P1clust22sided1", "P1clust2min2sided1")
m3 = c("P1sel332sided1", "P1sel33min2sided1", "P1iv332sided1", "P1iv33min2sided1", "P1var32sided1", "P1var3min2sided1", "P1clust32sided1", "P1clust3min2sided1")
M0=match(m0, b_names)
M1=match(m1, b_names)
M2=match(m2, b_names)
M3=match(m3, b_names)

S1 = match(names_for_scattters_sel, b_names)
S2 = match(names_for_scattters_sel_min, b_names)
S3 = match(names_for_scattters_iv, b_names)
S4 = match(names_for_scattters_iv_min, b_names)


