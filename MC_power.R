#####################################################################
# This function computes rejection rates in MC simulations
#a - non-p-hacked distribition, b - p-hacked, N - sample size, p_max_1/2 - subinterval upperbound
#bins_1/2 the number of bis to use in CS tests, M - number of MC replications, iter - dgp index
# Paper: Detecting p-hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
#####################################################################
MC <- function(a, b, N, p_max_2, bins_2, M,  iter){
  set.seed(123)
  alpha_avg = 0.05
  num_sz=1
  K <- 21
  
  fraction_a = mean(a<=0.15)
  fraction_b = mean(b<=0.15)
  fraction_diff = fraction_b - fraction_a
    
    B0_2 = Bound0_1s(p_max_2, bins_2)
    B1_2 = Bound1_1s(p_max_2, bins_2) 
    B2_2 = Bound2_1s(p_max_2, bins_2)
    
  LocBin = matrix(0, nrow=M*num_sz, ncol = K)
  rdd = matrix(0, nrow=M*num_sz, ncol = K)
  LCM8 = matrix(0, nrow=M*num_sz, ncol = K)
  FM = matrix(0, nrow=M*num_sz, ncol = K)
  CS_1 = matrix(0, nrow=M*num_sz, ncol = K)
  CS_UB = matrix(0, nrow=M*num_sz, ncol = K)
  CS_2B = matrix(0, nrow=M*num_sz, ncol = K)
  fail = matrix(0, nrow=M*num_sz, ncol = K)
  
  for (m in 1:M){
    p0 = sample(a, N, replace = TRUE)
    p1 = sample(b, N, replace = TRUE)
    d2 = runif(N, min = 0, max = 1)
    
    
    for (k in 1:K){
      tau = (k-1)/(K-1)
      d3 = 1*(d2<tau)
      P = (d3*p1 + (1-d3)*p0)
      
      LocBin[m,k] =LocBin[m,k]+ 1*(Binomial(P,0.04,0.05,"c")<=alpha_avg)
      rdd[m,k] =rdd[m,k]+ 1*(Discontinuity_test(P, 0.05)<=alpha_avg)
      LCM8[m,k] =LCM8[m,k]+ 1*(LCM(P,0,p_max_2)<=alpha_avg)
      FM[m,k] =FM[m,k]+ 1*(Fisher(P,0,p_max_2)<=alpha_avg)
      cs1 = CoxShi(P,0, 0,p_max_2, bins_2, 1, 0,0,0,0,1)
      CS_1[m,k] =CS_1[m,k]+ 1*(cs1<=alpha_avg)
      csub = CoxShi(P,0, 0,p_max_2, bins_2, 2, 1,B0_2,B1_2,B2_2,0) #modified
      CS_UB[m,k] =CS_UB[m,k]+ 1*(csub<=alpha_avg)
      cs2b = CoxShi(P,0, 0,p_max_2, bins_2, 2, 1,B0_2,B1_2,B2_2,1)
      CS_2B[m,k] =CS_2B[m,k]+ 1*(cs2b<=alpha_avg)
      
      fail[m,k] = 1*(cs1==999)+1*(cs1==888)
      
    }
  }
  

  LocBin = matrix(colMeans(LocBin), K, 1)
  LCM8 = matrix(colMeans(LCM8),K,1)
  FM = matrix(colMeans(FM),K,1)
  CS_1 = matrix(colMeans(CS_1),K,1)
  CS_UB = matrix(colMeans(CS_UB),K,1)
  CS_2B = matrix(colMeans(CS_2B),K,1)
  rdd = matrix(colMeans(rdd),K,1)
  fail = matrix(colMeans(fail),K,1)
  
  it = matrix(iter, size(rdd))
  
  #concat <- cbind(it,LocBin, FM, rdd, CS_1, CS_2B, LCM8, fail, LocBin_min, FM_min, rdd_min, CS_1_min, CS_2B_min, LCM8_min, fail_min)
  concat <- cbind(it,LocBin, FM, rdd, CS_1, CS_UB, CS_2B, LCM8, fail)
  #return(matrix(concat, K, 15))
  return(matrix(concat, K, 9))
}