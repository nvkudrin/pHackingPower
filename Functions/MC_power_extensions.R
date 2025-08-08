#####################################################################
# These functions compute rejection rates in MC simulations for Appendices F-I
# Paper: The power of tests for detecting p-hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
#####################################################################

#####################################################################
# Appendix F
#####################################################################
# rho - vector of correlation coefficients
# h_type - value of h (0,1,2,3)
# N - sample size
# p_max - right end of testing interval
# bins - number of bins for CS tests
# M - number of MC repetitions
# sided - 1 or 2 depending of whether researchers use 1- or 2-sided tests, respectively
# tau - value of tau (fraction of p-hackers)
#####################################################################
MC_examine_rho <- function(rho, h_type, N, p_max, bins, M, sided, tau) {

  # Set random seed for reproducibility
  set.seed(123)
  
  # Significance level
  alpha <- 0.05
  
  # Number of correlation values
  K <- length(rho)
  
  # Set boundary functions based on test sidedness
  if (sided == 1) {
    B0_2 <- Bound0_1s(p_max, bins)
    B1_2 <- Bound1_1s(p_max, bins)
    B2_2 <- Bound2_1s(p_max, bins)
  }
  
  if (sided == 2) {
    B0_2 <- Bound0(0, p_max, bins)
    B1_2 <- Bound1(0, p_max, bins)
    B2_2 <- Bound2(0, p_max, bins)
  }
  
  # Initialize result matrices for thresholding p-values
  LocBin <- matrix(0, nrow = M, ncol = K)
  rdd    <- matrix(0, nrow = M, ncol = K)
  LCM8   <- matrix(0, nrow = M, ncol = K)
  FM     <- matrix(0, nrow = M, ncol = K)
  CS_1   <- matrix(0, nrow = M, ncol = K)
  CS_UB  <- matrix(0, nrow = M, ncol = K)
  CS_2B  <- matrix(0, nrow = M, ncol = K)
  fail_1 <- matrix(0, nrow = M, ncol = K)
  fail_ub <- matrix(0, nrow = M, ncol = K)
  fail_2b <- matrix(0, nrow = M, ncol = K)
  
  # Initialize result matrices for minimum p-values
  LocBin_min <- matrix(0, nrow = M, ncol = K)
  rdd_min    <- matrix(0, nrow = M, ncol = K)
  LCM8_min   <- matrix(0, nrow = M, ncol = K)
  FM_min     <- matrix(0, nrow = M, ncol = K)
  CS_1_min   <- matrix(0, nrow = M, ncol = K)
  CS_UB_min  <- matrix(0, nrow = M, ncol = K)
  CS_2B_min  <- matrix(0, nrow = M, ncol = K)
  fail_1_min <- matrix(0, nrow = M, ncol = K)
  fail_ub_min <- matrix(0, nrow = M, ncol = K)
  fail_2b_min <- matrix(0, nrow = M, ncol = K)
  
  # Main Monte Carlo loop
  for (m in 1:M) {
    
    # Generate base random errors
    e0 <- rnorm(N, 0, 1)
    u0 <- rnorm(N, 0, 1)
    
    # Loop over correlation values
    for (k in 1:K) {
      
      # Generate correlated errors
      e1 <- rho[k] * e0 + sqrt(1 - rho[k]^2) * u0
      
      # Generate heterogeneity parameter
      if (h_type < 3) {
        h_sample <- h_type
      }
      
      if (h_type == 3) {
        h_sample <- rgamma(N, shape = 0.8547432, scale = 1.8690772)
      }
      
      # Generate test statistics and p-values
      T0 <- h_sample + e0
      T1 <- h_sample + e1
      P0 <- 1 - pnorm(T0)
      P1 <- 1 - pnorm(T1)
      
      # Threshold p-value selection
      p0 <- P0
      p1 <- P0 * (P0 <= 0.05) + P1 * (P0 > 0.05)
      
      # Minimum p-value selection
      p1min <- 0.5 * (P0 + P1 - abs(P0 - P1))
      
      # Generate selection indicator
      d2 <- runif(N, min = 0, max = 1)
      d3 <- 1 * (d2 < tau)
      
      # Final p-values
      P    <- (d3 * p1 + (1 - d3) * p0)
      Pmin <- (d3 * p1min + (1 - d3) * p0)
      
      # Test statistics for thresholding p-values
      LocBin[m, k] <- LocBin[m, k] + 1 * (Binomial(P, 0.04, 0.05, "c") <= alpha)
      rdd[m, k]    <- rdd[m, k] + 1 * (Discontinuity_test(P, 0.05) <= alpha)
      LCM8[m, k]   <- LCM8[m, k] + 1 * (LCM(P, 0, p_max) <= alpha)
      FM[m, k]     <- FM[m, k] + 1 * (Fisher(P, 0, p_max) <= alpha)
      
      # Cox-Shi tests for thresholding p-values
      cs1 <- CoxShi(P, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
      CS_1[m, k] <- CS_1[m, k] + 1 * (cs1 <= alpha)
      
      csub <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0) # modified
      CS_UB[m, k] <- CS_UB[m, k] + 1 * (csub <= alpha)
      
      cs2b <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
      CS_2B[m, k] <- CS_2B[m, k] + 1 * (cs2b <= alpha)
      
      # Track failures for thresholding p-values
      fail_1[m, k] <- 1 * (cs1 == 999) + 1 * (cs1 == 888)
      fail_ub[m, k] <- 1 * (csub == 999) + 1 * (csub == 888)
      fail_2b[m, k] <- 1 * (cs2b == 999) + 1 * (cs2b == 888)
      
      # Test statistics for minimum p-values
      LocBin_min[m, k] <- LocBin_min[m, k] + 1 * (Binomial(Pmin, 0.04, 0.05, "c") <= alpha)
      rdd_min[m, k]    <- rdd_min[m, k] + 1 * (Discontinuity_test(Pmin, 0.05) <= alpha)
      LCM8_min[m, k]   <- LCM8_min[m, k] + 1 * (LCM(Pmin, 0, p_max) <= alpha)
      FM_min[m, k]     <- FM_min[m, k] + 1 * (Fisher(Pmin, 0, p_max) <= alpha)
      
      # Cox-Shi tests for minimum p-values
      cs1_min <- CoxShi(Pmin, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
      CS_1_min[m, k] <- CS_1_min[m, k] + 1 * (cs1_min <= alpha)
      
      csub_min <- CoxShi(Pmin, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0) # modified
      CS_UB_min[m, k] <- CS_UB_min[m, k] + 1 * (csub_min <= alpha)
      
      cs2b_min <- CoxShi(Pmin, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
      CS_2B_min[m, k] <- CS_2B_min[m, k] + 1 * (cs2b_min <= alpha)
      
      # Track failures for minimum p-values
      fail_1_min[m, k] <- 1 * (cs1_min == 999) + 1 * (cs1_min == 888)
      fail_ub_min[m, k] <- 1 * (csub_min == 999) + 1 * (csub_min == 888)
      fail_2b_min[m, k] <- 1 * (cs2b_min == 999) + 1 * (cs2b_min == 888)
    }
  }
  
  # Calculate averages across Monte Carlo iterations for thresholding p-values
  LocBin  <- matrix(colMeans(LocBin), K, 1)
  LCM8    <- matrix(colMeans(LCM8), K, 1)
  FM      <- matrix(colMeans(FM), K, 1)
  CS_1    <- matrix(colMeans(CS_1), K, 1)
  CS_UB   <- matrix(colMeans(CS_UB), K, 1)
  CS_2B   <- matrix(colMeans(CS_2B), K, 1)
  rdd     <- matrix(colMeans(rdd), K, 1)
  fail_1  <- matrix(colMeans(fail_1), K, 1)
  fail_ub <- matrix(colMeans(fail_ub), K, 1)
  fail_2b <- matrix(colMeans(fail_2b), K, 1)
  
  # Calculate averages across Monte Carlo iterations for minimum p-values
  LocBin_min  <- matrix(colMeans(LocBin_min), K, 1)
  LCM8_min    <- matrix(colMeans(LCM8_min), K, 1)
  FM_min      <- matrix(colMeans(FM_min), K, 1)
  CS_1_min    <- matrix(colMeans(CS_1_min), K, 1)
  CS_UB_min   <- matrix(colMeans(CS_UB_min), K, 1)
  CS_2B_min   <- matrix(colMeans(CS_2B_min), K, 1)
  rdd_min     <- matrix(colMeans(rdd_min), K, 1)
  fail_1_min  <- matrix(colMeans(fail_1_min), K, 1)
  fail_ub_min <- matrix(colMeans(fail_ub_min), K, 1)
  fail_2b_min <- matrix(colMeans(fail_2b_min), K, 1)
  
  # Combine all results into output matrix
  concat <- cbind(rho, LocBin, FM, rdd, CS_1, CS_UB, CS_2B, LCM8, 
                  fail_1, fail_ub, fail_2b, LocBin_min, FM_min, rdd_min, 
                  CS_1_min, CS_UB_min, CS_2B_min, LCM8_min, 
                  fail_1_min, fail_ub_min, fail_2b_min)
  
  return(matrix(concat, K, 21))
}


#####################################################################
# Appendix G
#####################################################################
# pmis - vector of misclassification probabilities
# h_type - value of h (0,1,2,3)
# N - sample size
# p_max - right end of testing interval
# bins - number of bins for CS tests
# M - number of MC repetitions
# sided - 1 or 2 depending of whether researchers use 1- or 2-sided tests, respectively
# rho - correlation between t-stats
# tau - value of tau (fraction of p-hackers)
#####################################################################
MC_examine_pmis <- function(pmis, h_type, N, p_max, bins, M, sided, rho, tau) {
  
  # Set random seed for reproducibility
  set.seed(123)
  
  # Significance level
  alpha <- 0.05
  
  # Number of misclassification probability values
  K <- length(pmis)
  
  # Set boundary functions based on test sidedness
  if (sided == 1) {
    B0_2 <- Bound0_1s(p_max, bins)
    B1_2 <- Bound1_1s(p_max, bins)
    B2_2 <- Bound2_1s(p_max, bins)
  }
  
  if (sided == 2) {
    B0_2 <- Bound0(0, p_max, bins)
    B1_2 <- Bound1(0, p_max, bins)
    B2_2 <- Bound2(0, p_max, bins)
  }
  
  # Initialize result matrices for thresholding p-values
  LocBin <- matrix(0, nrow = M, ncol = K)
  rdd    <- matrix(0, nrow = M, ncol = K)
  LCM8   <- matrix(0, nrow = M, ncol = K)
  FM     <- matrix(0, nrow = M, ncol = K)
  CS_1   <- matrix(0, nrow = M, ncol = K)
  CS_UB  <- matrix(0, nrow = M, ncol = K)
  CS_2B  <- matrix(0, nrow = M, ncol = K)
  fail_1 <- matrix(0, nrow = M, ncol = K)
  fail_ub <- matrix(0, nrow = M, ncol = K)
  fail_2b <- matrix(0, nrow = M, ncol = K)
  
  # Initialize result matrices for minimum p-values
  LocBin_min <- matrix(0, nrow = M, ncol = K)
  rdd_min    <- matrix(0, nrow = M, ncol = K)
  LCM8_min   <- matrix(0, nrow = M, ncol = K)
  FM_min     <- matrix(0, nrow = M, ncol = K)
  CS_1_min   <- matrix(0, nrow = M, ncol = K)
  CS_UB_min  <- matrix(0, nrow = M, ncol = K)
  CS_2B_min  <- matrix(0, nrow = M, ncol = K)
  fail_1_min <- matrix(0, nrow = M, ncol = K)
  fail_ub_min <- matrix(0, nrow = M, ncol = K)
  fail_2b_min <- matrix(0, nrow = M, ncol = K)
  
  # Set up correlation structure for instrumental variable approach
  gamma_z <- sqrt(1 - rho)
  R <- matrix(c(1, gamma_z, gamma_z, 
                gamma_z, 1, gamma_z^2, 
                gamma_z, gamma_z^2, 1), 3, 3)
  
  # Main Monte Carlo loop
  for (m in 1:M) {
    
    # Generate correlated random variables using Cholesky decomposition
    WXU  <- rnorm(N, 0, 1)
    WZ1U <- rnorm(N, 0, 1)
    WZ2U <- rnorm(N, 0, 1)
    Ws   <- matrix(c(WXU, WZ1U, WZ2U), N, 3) %*% chol(R)
    
    # Generate heterogeneity parameter
    if (h_type < 3) {
      h_sample <- h_type
    }
    
    if (h_type == 3) {
      h_sample <- rgamma(N, shape = 0.8547432, scale = 1.8690772)
    }
    
    # Generate test statistics for main outcomes
    T0 <- h_sample + (Ws[, 1] - gamma_z * Ws[, 2]) / sqrt(rho)
    T1 <- h_sample + (Ws[, 1] - gamma_z * Ws[, 3]) / sqrt(rho)
    
    # Generate test statistics for control variables
    Tz0 <- (Ws[, 2] - gamma_z * Ws[, 1]) / sqrt(rho)
    Tz1 <- (Ws[, 3] - gamma_z * Ws[, 1]) / sqrt(rho)
    
    # Convert to two-sided p-values
    P0 <- 2 * (1 - pnorm(abs(T0)))
    P1 <- 2 * (1 - pnorm(abs(T1)))
    p0 <- P0
    
    # Threshold p-value selection for main outcomes
    p1 <- P0 * (P0 <= 0.05) + P1 * (P0 > 0.05)
    
    # Generate p-values
    p_z0 <- 2 * (1 - pnorm(abs(Tz0)))
    T_z  <- Tz0 * (p0 == p1) + Tz1 * (p0 != p1)
    p_z1 <- 2 * (1 - pnorm(abs(T_z)))
    
    # Minimum p-value selection for main outcomes
    p0     <- P0
    p1min  <- 0.5 * (P0 + P1 - abs(P0 - P1))
    p_z0_min <- 2 * (1 - pnorm(abs(Tz0)))
    T_z    <- Tz0 * (p0 == p1min) + Tz1 * (p0 != p1min)
    p_z1_min <- 2 * (1 - pnorm(abs(T_z)))
    
    # Generate selection indicator
    d2 <- runif(N, min = 0, max = 1)
    d3 <- 1 * (d2 < tau)
    
    # Generate misclassification indicators
    d_mis <- runif(N, min = 0, max = 1)
    
    # Loop over misclassification probability values
    for (k in 1:K) {
      misclass <- 1*(d_mis<=pmis[k])
      
    
    # Apply misclassification to thresholding p-values
    p0_t <- p_z0 * (misclass) + p0 * (1 - misclass)
    p1_t <- p_z1 * (misclass) + p1 * (1 - misclass)
    
    
    # Apply misclassification to minimum p-values
    p0_m <- p_z0_min * (misclass) + p0 * (1 - misclass)
    p1_m <- p_z1_min * (misclass) + p1min * (1 - misclass)
    
      
      # Final p-values with selection
      P    <- (d3 * p1_t + (1 - d3) * p0_t)
      Pmin <- (d3 * p1_m + (1 - d3) * p0_m)
      
      # Test statistics for thresholding p-values
      LocBin[m, k] <- LocBin[m, k] + 1 * (Binomial(P, 0.04, 0.05, "c") <= alpha)
      rdd[m, k]    <- rdd[m, k] + 1 * (Discontinuity_test(P, 0.05) <= alpha)
      LCM8[m, k]   <- LCM8[m, k] + 1 * (LCM(P, 0, p_max) <= alpha)
      FM[m, k]     <- FM[m, k] + 1 * (Fisher(P, 0, p_max) <= alpha)
      
      # Cox-Shi tests for thresholding p-values
      cs1 <- CoxShi(P, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
      CS_1[m, k] <- CS_1[m, k] + 1 * (cs1 <= alpha)
      
      csub <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0) # modified
      CS_UB[m, k] <- CS_UB[m, k] + 1 * (csub <= alpha)
      
      cs2b <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
      CS_2B[m, k] <- CS_2B[m, k] + 1 * (cs2b <= alpha)
      
      # Track failures for thresholding p-values
      fail_1[m, k]  <- 1 * (cs1 == 999) + 1 * (cs1 == 888)
      fail_ub[m, k] <- 1 * (csub == 999) + 1 * (csub == 888)
      fail_2b[m, k] <- 1 * (cs2b == 999) + 1 * (cs2b == 888)
      
      # Test statistics for minimum p-values
      LocBin_min[m, k] <- LocBin_min[m, k] + 1 * (Binomial(Pmin, 0.04, 0.05, "c") <= alpha)
      rdd_min[m, k]    <- rdd_min[m, k] + 1 * (Discontinuity_test(Pmin, 0.05) <= alpha)
      LCM8_min[m, k]   <- LCM8_min[m, k] + 1 * (LCM(Pmin, 0, p_max) <= alpha)
      FM_min[m, k]     <- FM_min[m, k] + 1 * (Fisher(Pmin, 0, p_max) <= alpha)
      
      # Cox-Shi tests for minimum p-values
      cs1_min <- CoxShi(Pmin, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
      CS_1_min[m, k] <- CS_1_min[m, k] + 1 * (cs1_min <= alpha)
      
      csub_min <- CoxShi(Pmin, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0) # modified
      CS_UB_min[m, k] <- CS_UB_min[m, k] + 1 * (csub_min <= alpha)
      
      cs2b_min <- CoxShi(Pmin, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
      CS_2B_min[m, k] <- CS_2B_min[m, k] + 1 * (cs2b_min <= alpha)
      
      # Track failures for minimum p-values
      fail_1_min[m, k]  <- 1 * (cs1_min == 999) + 1 * (cs1_min == 888)
      fail_ub_min[m, k] <- 1 * (csub_min == 999) + 1 * (csub_min == 888)
      fail_2b_min[m, k] <- 1 * (cs2b_min == 999) + 1 * (cs2b_min == 888)
    }
  }
  
  # Calculate averages across Monte Carlo iterations for threshold p-values
  LocBin  <- matrix(colMeans(LocBin), K, 1)
  LCM8    <- matrix(colMeans(LCM8), K, 1)
  FM      <- matrix(colMeans(FM), K, 1)
  CS_1    <- matrix(colMeans(CS_1), K, 1)
  CS_UB   <- matrix(colMeans(CS_UB), K, 1)
  CS_2B   <- matrix(colMeans(CS_2B), K, 1)
  rdd     <- matrix(colMeans(rdd), K, 1)
  fail_1  <- matrix(colMeans(fail_1), K, 1)
  fail_ub <- matrix(colMeans(fail_ub), K, 1)
  fail_2b <- matrix(colMeans(fail_2b), K, 1)
  
  # Calculate averages across Monte Carlo iterations for minimum p-values
  LocBin_min  <- matrix(colMeans(LocBin_min), K, 1)
  LCM8_min    <- matrix(colMeans(LCM8_min), K, 1)
  FM_min      <- matrix(colMeans(FM_min), K, 1)
  CS_1_min    <- matrix(colMeans(CS_1_min), K, 1)
  CS_UB_min   <- matrix(colMeans(CS_UB_min), K, 1)
  CS_2B_min   <- matrix(colMeans(CS_2B_min), K, 1)
  rdd_min     <- matrix(colMeans(rdd_min), K, 1)
  fail_1_min  <- matrix(colMeans(fail_1_min), K, 1)
  fail_ub_min <- matrix(colMeans(fail_ub_min), K, 1)
  fail_2b_min <- matrix(colMeans(fail_2b_min), K, 1)
  
  # Combine all results into output matrix
  concat <- cbind(pmis, LocBin, FM, rdd, CS_1, CS_UB, CS_2B, LCM8, 
                  fail_1, fail_ub, fail_2b, LocBin_min, FM_min, rdd_min, 
                  CS_1_min, CS_UB_min, CS_2B_min, LCM8_min, 
                  fail_1_min, fail_ub_min, fail_2b_min)
  
  return(matrix(concat, K, 21))
}

#####################################################################
# Appendix H
#####################################################################
# P0 - null distributions
# P1 - p-hacked thresholding distribution
# P1min - p-hacked minimum distribution
# N - sample size
# p_max - right end of testing interval
# bins_small_n - number of bins for CS tests with small sample
# bins_large_n - number of bins for CS tests with large sample
# M - number of MC repetitions
# sided - 1 or 2 depending of whether researchers use 1- or 2-sided tests, respectively
# tau - value of tau (fraction of p-hackers)
#####################################################################
MC_examine_n <- function(P0, P1, P1min, N, p_max, bins_small_n, bins_large_n, M, sided, tau) {
  
  # Set random seed for reproducibility
  set.seed(123)
  
  # Significance level
  alpha <- 0.05
  
  # Number of sample size values to examine
  K <- length(N)
  maxN <- max(N)
  
  # Set boundary functions based on test sidedness
  if (sided == 1) {
    B0_2_small <- Bound0_1s(p_max, bins_small_n)
    B1_2_small <- Bound1_1s(p_max, bins_small_n)
    B2_2_small <- Bound2_1s(p_max, bins_small_n)
  }
  
  if (sided == 2) {
    B0_2_small <- Bound0(0, p_max, bins_small_n)
    B1_2_small <- Bound1(0, p_max, bins_small_n)
    B2_2_small <- Bound2(0, p_max, bins_small_n)
  }
  
  if (sided == 1) {
    B0_2_large <- Bound0_1s(p_max, bins_small_n)
    B1_2_large <- Bound1_1s(p_max, bins_small_n)
    B2_2_large <- Bound2_1s(p_max, bins_small_n)
  }
  
  if (sided == 2) {
    B0_2_large <- Bound0(0, p_max, bins_large_n)
    B1_2_large <- Bound1(0, p_max, bins_large_n)
    B2_2_large <- Bound2(0, p_max, bins_large_n)
  }
  
  # Initialize result matrices for thresholding p-values
  LocBin  <- matrix(0, nrow = M, ncol = K)
  Disc    <- matrix(0, nrow = M, ncol = K)
  LCM8    <- matrix(0, nrow = M, ncol = K)
  FM      <- matrix(0, nrow = M, ncol = K)
  CS_1    <- matrix(0, nrow = M, ncol = K)
  CS_UB   <- matrix(0, nrow = M, ncol = K)
  CS_2B   <- matrix(0, nrow = M, ncol = K)
  fail_1  <- matrix(0, nrow = M, ncol = K)
  fail_ub <- matrix(0, nrow = M, ncol = K)
  fail_2b <- matrix(0, nrow = M, ncol = K)
  
  # Initialize result matrices for minimum p-values
  LocBin_min  <- matrix(0, nrow = M, ncol = K)
  Disc_min    <- matrix(0, nrow = M, ncol = K)
  LCM8_min    <- matrix(0, nrow = M, ncol = K)
  FM_min      <- matrix(0, nrow = M, ncol = K)
  CS_1_min    <- matrix(0, nrow = M, ncol = K)
  CS_UB_min   <- matrix(0, nrow = M, ncol = K)
  CS_2B_min   <- matrix(0, nrow = M, ncol = K)
  fail_1_min  <- matrix(0, nrow = M, ncol = K)
  fail_ub_min <- matrix(0, nrow = M, ncol = K)
  fail_2b_min <- matrix(0, nrow = M, ncol = K)
  
  # Main Monte Carlo loop
  for (m in 1:M) {
    
    # Sample indices to draw p-values from input distributions
    draw  <- randi(length(P0), maxN, 1)
    p0    <- P0[draw]
    p1    <- P1[draw]
    p1min <- P1min[draw]
    
    # Generate selection indicators
    d2     <- runif(maxN, min = 0, max = 1)
    d3     <- 1 * (d2 < tau)
    
    # Create final p-value vectors with selection
    P     <- d3 * p1 + (1 - d3) * p0
    P_min <- d3 * p1min + (1 - d3) * p0
    
    # Loop over different sample sizes
    for (k in 1:K) {
      if (N[k]>1000){
        B0_2 = B0_2_large
        B1_2 = B1_2_large 
        B2_2 = B2_2_large
        bins = bins_large_n
      }
      if (N[k]<=1000){
        B0_2 = B0_2_small
        B1_2 = B1_2_small 
        B2_2 = B2_2_small
        bins = bins_small_n
      }
      # Test statistics for thresholding p-values (using first N[k] observations)
      LocBin[m, k] <- LocBin[m, k] + 1 * (Binomial(P[1:N[k]], 0.04, 0.05, "c") <= alpha)
      Disc[m, k]   <- Disc[m, k] + 1 * (Discontinuity_test(P[1:N[k]], 0.05) <= alpha)
      LCM8[m, k]   <- LCM8[m, k] + 1 * (LCM(P[1:N[k]], 0, p_max) <= alpha)
      FM[m, k]     <- FM[m, k] + 1 * (Fisher(P[1:N[k]], 0, p_max) <= alpha)
      
      # Cox-Shi tests for thresholding p-values
      cs1 <- CoxShi(P[1:N[k]], 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
      CS_1[m, k] <- CS_1[m, k] + 1 * (cs1 <= alpha)
      
      csub <- CoxShi(P[1:N[k]], 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0)
      CS_UB[m, k] <- CS_UB[m, k] + 1 * (csub <= alpha)
      
      cs2b <- CoxShi(P[1:N[k]], 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
      CS_2B[m, k] <- CS_2B[m, k] + 1 * (cs2b <= alpha)
      
      # Track failures for thresholding p-values
      fail_1[m, k]  <- 1 * (cs1 == 999) + 1 * (cs1 == 888)
      fail_ub[m, k] <- 1 * (csub == 999) + 1 * (csub == 888)
      fail_2b[m, k] <- 1 * (cs2b == 999) + 1 * (cs2b == 888)
      
      # Test statistics for minimum p-values (using first N[k] observations)
      LocBin_min[m, k] <- LocBin_min[m, k] + 1 * (Binomial(P_min[1:N[k]], 0.04, 0.05, "c") <= alpha)
      Disc_min[m, k]   <- Disc_min[m, k] + 1 * (Discontinuity_test(P_min[1:N[k]], 0.05) <= alpha)
      LCM8_min[m, k]   <- LCM8_min[m, k] + 1 * (LCM(P_min[1:N[k]], 0, p_max) <= alpha)
      FM_min[m, k]     <- FM_min[m, k] + 1 * (Fisher(P_min[1:N[k]], 0, p_max) <= alpha)
      
      # Cox-Shi tests for minimum p-values
      cs1_min <- CoxShi(P_min[1:N[k]], 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
      CS_1_min[m, k] <- CS_1_min[m, k] + 1 * (cs1_min <= alpha)
      
      csub_min <- CoxShi(P_min[1:N[k]], 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0)
      CS_UB_min[m, k] <- CS_UB_min[m, k] + 1 * (csub_min <= alpha)
      
      cs2b_min <- CoxShi(P_min[1:N[k]], 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
      CS_2B_min[m, k] <- CS_2B_min[m, k] + 1 * (cs2b_min <= alpha)
      
      # Track failures for minimum p-values
      fail_1_min[m, k]  <- 1 * (cs1_min == 999) + 1 * (cs1_min == 888)
      fail_ub_min[m, k] <- 1 * (csub_min == 999) + 1 * (csub_min == 888)
      fail_2b_min[m, k] <- 1 * (cs2b_min == 999) + 1 * (cs2b_min == 888)
    }
  }
  
  # Calculate averages across Monte Carlo iterations for thresholding p-values
  LocBin  <- matrix(colMeans(LocBin), K, 1)
  LCM8    <- matrix(colMeans(LCM8), K, 1)
  FM      <- matrix(colMeans(FM), K, 1)
  CS_1    <- matrix(colMeans(CS_1), K, 1)
  CS_UB   <- matrix(colMeans(CS_UB), K, 1)
  CS_2B   <- matrix(colMeans(CS_2B), K, 1)
  Disc    <- matrix(colMeans(Disc), K, 1)
  fail_1  <- matrix(colMeans(fail_1), K, 1)
  fail_ub <- matrix(colMeans(fail_ub), K, 1)
  fail_2b <- matrix(colMeans(fail_2b), K, 1)
  
  # Calculate averages across Monte Carlo iterations for minimum p-values
  LocBin_min  <- matrix(colMeans(LocBin_min), K, 1)
  LCM8_min    <- matrix(colMeans(LCM8_min), K, 1)
  FM_min      <- matrix(colMeans(FM_min), K, 1)
  CS_1_min    <- matrix(colMeans(CS_1_min), K, 1)
  CS_UB_min   <- matrix(colMeans(CS_UB_min), K, 1)
  CS_2B_min   <- matrix(colMeans(CS_2B_min), K, 1)
  Disc_min    <- matrix(colMeans(Disc_min), K, 1)
  fail_1_min  <- matrix(colMeans(fail_1_min), K, 1)
  fail_ub_min <- matrix(colMeans(fail_ub_min), K, 1)
  fail_2b_min <- matrix(colMeans(fail_2b_min), K, 1)
  
  # Combine all results into output matrix
  concat <- cbind(N, LocBin, FM, Disc, CS_1, CS_UB, CS_2B, LCM8, 
                  fail_1, fail_ub, fail_2b, LocBin_min, FM_min, Disc_min, 
                  CS_1_min, CS_UB_min, CS_2B_min, LCM8_min, 
                  fail_1_min, fail_ub_min, fail_2b_min)
  
  return(matrix(concat, K, 21))
}



#####################################################################
# Appendix I
# P0 - null distributions
# P1 - p-hacked thresholding distribution
# P1min - p-hacked minimum distribution
# N - sample size
# p_max - right end of testing interval
# bins - number of bins for CS tests
# M - number of MC repetitions
# sided - 1 or 2 depending of whether researchers use 1- or 2-sided tests, respectively
# tau_t - value of tau (fraction of p-hackers) for thresholding example
# tau_min - value of tau (fraction of p-hackers) for minimum example
#####################################################################
MC_for_power_bias <- function(P0, P1, P1min, N, p_max, bins, M, sided, tau_t, tau_min) {
  set.seed(123)
  alpha <- 0.05  # significance level
  
  K <- 1
  
  # Calculate bounds on proportions for 1-/2-sided tests
  if (sided == 1) {
    B0_2 <- Bound0_1s(p_max, bins)
    B1_2 <- Bound1_1s(p_max, bins)
    B2_2 <- Bound2_1s(p_max, bins)
  }
  if (sided == 2) {
    B0_2 <- Bound0(0, p_max, bins)
    B1_2 <- Bound1(0, p_max, bins)
    B2_2 <- Bound2(0, p_max, bins)
  }
  
  # Initialize output matrices
  LocBin      <- matrix(0, nrow = M, ncol = 1)
  Disc        <- matrix(0, nrow = M, ncol = 1)
  LCM8        <- matrix(0, nrow = M, ncol = 1)
  FM          <- matrix(0, nrow = M, ncol = 1)
  CS_1        <- matrix(0, nrow = M, ncol = 1)
  CS_UB       <- matrix(0, nrow = M, ncol = 1)
  CS_2B       <- matrix(0, nrow = M, ncol = 1)
  fail_1      <- matrix(0, nrow = M, ncol = 1)
  fail_ub     <- matrix(0, nrow = M, ncol = 1)
  fail_2b     <- matrix(0, nrow = M, ncol = 1)
  
  LocBin_min  <- matrix(0, nrow = M, ncol = 1)
  Disc_min    <- matrix(0, nrow = M, ncol = 1)
  LCM8_min    <- matrix(0, nrow = M, ncol = 1)
  FM_min      <- matrix(0, nrow = M, ncol = 1)
  CS_1_min    <- matrix(0, nrow = M, ncol = 1)
  CS_UB_min   <- matrix(0, nrow = M, ncol = 1)
  CS_2B_min   <- matrix(0, nrow = M, ncol = 1)
  fail_1_min  <- matrix(0, nrow = M, ncol = 1)
  fail_ub_min <- matrix(0, nrow = M, ncol = 1)
  fail_2b_min <- matrix(0, nrow = M, ncol = 1)
  
  for (m in 1:M) {
    draw    <- randi(length(P0), N, 1)
    p0      <- P0[draw]
    p1      <- P1[draw]
    p1min   <- P1min[draw]
    
    
      # Generate selection indicators
      d2     <- runif(N, min = 0, max = 1)
      d3      <- 1 * (d2 < tau_t)
      d3_min      <- 1 * (d2 < tau_min)
      P       <- d3 * p1 + (1 - d3) * p0
      P_min   <- d3_min * p1min + (1 - d3_min) * p0

      
          LocBin[m, 1]      <- LocBin[m, 1]      + 1 * (Binomial(P, 0.04, 0.05, "c") <= alpha)
          Disc[m, 1]        <- Disc[m, 1]        + 1 * (Discontinuity_test(P, 0.05) <= alpha)
          LCM8[m, 1]        <- LCM8[m, 1]        + 1 * (LCM(P, 0, p_max) <= alpha)
          FM[m, 1]          <- FM[m, 1]          + 1 * (Fisher(P, 0, p_max) <= alpha)
          cs1                 <- CoxShi(P, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
          CS_1[m, 1]        <- CS_1[m, 1]        + 1 * (cs1 <= alpha)
          csub                <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0)
          CS_UB[m, 1]       <- CS_UB[m, 1]       + 1 * (csub <= alpha)
          cs2b                <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
          CS_2B[m, 1]       <- CS_2B[m, 1]       + 1 * (cs2b <= alpha)
          fail_1[m, 1]      <- 1 * (cs1 == 999)  + 1 * (cs1 == 888)
          fail_ub[m, 1]     <- 1 * (csub == 999) + 1 * (csub == 888)
          fail_2b[m, 1]     <- 1 * (cs2b == 999) + 1 * (cs2b == 888)
          
          LocBin_min[m, 1]  <- LocBin_min[m, 1]  + 1 * (Binomial(P_min, 0.04, 0.05, "c") <= alpha)
          Disc_min[m, 1]    <- Disc_min[m, 1]    + 1 * (Discontinuity_test(P_min, 0.05) <= alpha)
          LCM8_min[m, 1]    <- LCM8_min[m, 1]    + 1 * (LCM(P_min, 0, p_max) <= alpha)
          FM_min[m, 1]      <- FM_min[m, 1]      + 1 * (Fisher(P_min, 0, p_max) <= alpha)
          cs1_min             <- CoxShi(P_min, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
          CS_1_min[m, 1]    <- CS_1_min[m, 1]    + 1 * (cs1_min <= alpha)
          csub_min            <- CoxShi(P_min, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0)
          CS_UB_min[m, 1]   <- CS_UB_min[m, 1]   + 1 * (csub_min <= alpha)
          cs2b_min            <- CoxShi(P_min, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
          CS_2B_min[m, 1]   <- CS_2B_min[m, 1]   + 1 * (cs2b_min <= alpha)
          fail_1_min[m, 1]  <- 1 * (cs1_min == 999)  + 1 * (cs1_min == 888)
          fail_ub_min[m, 1] <- 1 * (csub_min == 999) + 1 * (csub_min == 888)
          fail_2b_min[m, 1] <- 1 * (cs2b_min == 999) + 1 * (cs2b_min == 888)
        
  }

  
  # Aggregate results (means over Monte Carlo reps)
  LocBin      <- matrix(colMeans(LocBin),      K, 1)
  LCM8        <- matrix(colMeans(LCM8),        K, 1)
  FM          <- matrix(colMeans(FM),          K, 1)
  CS_1        <- matrix(colMeans(CS_1),        K, 1)
  CS_UB       <- matrix(colMeans(CS_UB),       K, 1)
  CS_2B       <- matrix(colMeans(CS_2B),       K, 1)
  Disc        <- matrix(colMeans(Disc),        K, 1)
  fail_1      <- matrix(colMeans(fail_1),      K, 1)
  fail_ub     <- matrix(colMeans(fail_ub),     K, 1)
  fail_2b     <- matrix(colMeans(fail_2b),     K, 1)
  
  LocBin_min  <- matrix(colMeans(LocBin_min),  K, 1)
  LCM8_min    <- matrix(colMeans(LCM8_min),    K, 1)
  FM_min      <- matrix(colMeans(FM_min),      K, 1)
  CS_1_min    <- matrix(colMeans(CS_1_min),    K, 1)
  CS_UB_min   <- matrix(colMeans(CS_UB_min),   K, 1)
  CS_2B_min   <- matrix(colMeans(CS_2B_min),   K, 1)
  Disc_min    <- matrix(colMeans(Disc_min),    K, 1)
  fail_1_min  <- matrix(colMeans(fail_1_min),  K, 1)
  fail_ub_min <- matrix(colMeans(fail_ub_min), K, 1)
  fail_2b_min <- matrix(colMeans(fail_2b_min), K, 1)
  
  concat1 <- cbind(LocBin, FM, Disc, CS_1, CS_UB, CS_2B, LCM8, fail_1, fail_ub, fail_2b)
  concat2 <- cbind(LocBin_min, FM_min, Disc_min, CS_1_min, CS_UB_min, CS_2B_min, LCM8_min, fail_1_min, fail_ub_min, fail_2b_min)
  
  concat <- rbind(concat1, concat2)
  return(matrix(concat, 2, 10))
}
