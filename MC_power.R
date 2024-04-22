#####################################################################
# This function computes rejection rates in MC simulations
#a - non-p-hacked distribition, b - p-hacked, N - sample size, p_max_1/2 - subinterval upperbound
#bins_1/2 the number of bis to use in CS tests, M - number of MC replications, iter - dgp index
# Paper: Detecting p-hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
#####################################################################
MC <- function(a, b, N, p_max_2, bins_2, M, sided, iter){
  set.seed(123)
  alpha_avg = 0.05
  num_sz=1
  K <- 21
  
  fraction_a = mean(a<=0.15)
  fraction_b = mean(b<=0.15)
  fraction_diff = fraction_b - fraction_a
    
  if (sided == 1){
    B0_2 = Bound0_1s(p_max_2, bins_2)
    B1_2 = Bound1_1s(p_max_2, bins_2) 
    B2_2 = Bound2_1s(p_max_2, bins_2)
  }
  if (sided == 2){
    B0_2 = Bound0(0,p_max_2, bins_2)
    B1_2 = Bound1(0,p_max_2, bins_2) 
    B2_2 = Bound2(0,p_max_2, bins_2)
  }
  K = 3*K
  LocBin = matrix(0, nrow=M*num_sz, ncol = K)
  rdd = matrix(0, nrow=M*num_sz, ncol = K)
  LCM8 = matrix(0, nrow=M*num_sz, ncol = K)
  FM = matrix(0, nrow=M*num_sz, ncol = K)
  CS_1 = matrix(0, nrow=M*num_sz, ncol = K)
  CS_UB = matrix(0, nrow=M*num_sz, ncol = K)
  CS_2B = matrix(0, nrow=M*num_sz, ncol = K)
  fail_1 = matrix(0, nrow=M*num_sz, ncol = K)
  fail_ub = matrix(0, nrow=M*num_sz, ncol = K)
  fail_2b = matrix(0, nrow=M*num_sz, ncol = K)
  K = K/3
  for (m in 1:M){
    p0 = sample(a, N, replace = TRUE)
    p1 = sample(b, N, replace = TRUE)
    d2 = runif(N, min = 0, max = 1)
    
    pub_bias_draw = runif(N, min = 0, max = 1)
    pub_bias_sharp = 1*(pub_bias_draw<0.1)
    
    
    for (k in 1:K){
      tau = (k-1)/(K-1)
      d3 = 1*(d2<tau)
      P = (d3*p1 + (1-d3)*p0)
      
      publish_sharp = 1*(P<=0.05) + (P>0.05)*pub_bias_sharp
      
      P_published_sharp = P[publish_sharp == 1]
      
      P_published_smooth = P[pub_bias_draw < exp(-8.45*P)]
      
      for (pb_mode in 1:3){
        
        if ((pb_mode == 1)|(tau == 0.5)){
        if (pb_mode == 2){
          P = P_published_sharp
        }
        if (pb_mode == 3){
          P = P_published_smooth
        }
      LocBin[m,k+(pb_mode-1)*K] =LocBin[m,k+(pb_mode-1)*K]+ 1*(Binomial(P,0.04,0.05,"c")<=alpha_avg)
      rdd[m,k+(pb_mode-1)*K] =rdd[m,k+(pb_mode-1)*K]+ 1*(Discontinuity_test(P, 0.05)<=alpha_avg)
      LCM8[m,k+(pb_mode-1)*K] =LCM8[m,k+(pb_mode-1)*K]+ 1*(LCM(P,0,p_max_2)<=alpha_avg)
      FM[m,k+(pb_mode-1)*K] =FM[m,k+(pb_mode-1)*K]+ 1*(Fisher(P,0,p_max_2)<=alpha_avg)
      cs1 = CoxShi(P,0, 0,p_max_2, bins_2, 1, 0,0,0,0,1)
      CS_1[m,k+(pb_mode-1)*K] =CS_1[m,k+(pb_mode-1)*K]+ 1*(cs1<=alpha_avg)
      csub = CoxShi(P,0, 0,p_max_2, bins_2, 2, 1,B0_2,B1_2,B2_2,0) #modified
      CS_UB[m,k+(pb_mode-1)*K] =CS_UB[m,k+(pb_mode-1)*K]+ 1*(csub<=alpha_avg)
      cs2b = CoxShi(P,0, 0,p_max_2, bins_2, 2, 1,B0_2,B1_2,B2_2,1)
      CS_2B[m,k+(pb_mode-1)*K] =CS_2B[m,k+(pb_mode-1)*K]+ 1*(cs2b<=alpha_avg)
      fail_1[m,k+(pb_mode-1)*K] = 1*(cs1==999)+1*(cs1==888)
      fail_ub[m,k+(pb_mode-1)*K] = 1*(csub==999)+1*(csub==888)
      fail_2b[m,k+(pb_mode-1)*K] = 1*(cs2b==999)+1*(cs2b==888)
      }
      }
    }
  }
  
  K = K*3

  LocBin = matrix(colMeans(LocBin), K, 1)
  LCM8 = matrix(colMeans(LCM8),K,1)
  FM = matrix(colMeans(FM),K,1)
  CS_1 = matrix(colMeans(CS_1),K,1)
  CS_UB = matrix(colMeans(CS_UB),K,1)
  CS_2B = matrix(colMeans(CS_2B),K,1)
  rdd = matrix(colMeans(rdd),K,1)
  fail_1 = matrix(colMeans(fail_1),K,1)
  fail_ub = matrix(colMeans(fail_ub),K,1)
  fail_2b = matrix(colMeans(fail_2b),K,1)
  it = matrix(iter, size(rdd))
  
  #concat <- cbind(it,LocBin, FM, rdd, CS_1, CS_2B, LCM8, fail, LocBin_min, FM_min, rdd_min, CS_1_min, CS_2B_min, LCM8_min, fail_min)
  concat <- cbind(it,LocBin, FM, rdd, CS_1, CS_UB, CS_2B, LCM8, fail_1, fail_ub, fail_2b)
  #return(matrix(concat, K, 15))
  return(matrix(concat, K, 11))
}