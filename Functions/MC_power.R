#####################################################################
# This function computes rejection rates in MC simulations
# Paper: The power of tests for detecting p-hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
#####################################################################
# P0 - non-p-hacked distribition
# P1 - p-hacked thresholding 
# P1min - p-hacked minimum
# N - sample size
# p_max - subinterval upperbound
# bins - the number of bins to use in CS tests
# M - number of MC replications
# Tau - vector of fractions of p-hackers
#####################################################################
MC <- function(P0, P1, P1min, N, p_max, bins, M, sided, Tau) {
  set.seed(123)
  alpha <- 0.05  # significance level
  
  K <- length(Tau)
  
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
  LocBin      <- matrix(0, nrow = M, ncol = 3 * K)
  Disc        <- matrix(0, nrow = M, ncol = 3 * K)
  LCM8        <- matrix(0, nrow = M, ncol = 3 * K)
  FM          <- matrix(0, nrow = M, ncol = 3 * K)
  CS_1        <- matrix(0, nrow = M, ncol = 3 * K)
  CS_UB       <- matrix(0, nrow = M, ncol = 3 * K)
  CS_2B       <- matrix(0, nrow = M, ncol = 3 * K)
  fail_1      <- matrix(0, nrow = M, ncol = 3 * K)
  fail_ub     <- matrix(0, nrow = M, ncol = 3 * K)
  fail_2b     <- matrix(0, nrow = M, ncol = 3 * K)
  
  LocBin_min  <- matrix(0, nrow = M, ncol = 3 * K)
  Disc_min    <- matrix(0, nrow = M, ncol = 3 * K)
  LCM8_min    <- matrix(0, nrow = M, ncol = 3 * K)
  FM_min      <- matrix(0, nrow = M, ncol = 3 * K)
  CS_1_min    <- matrix(0, nrow = M, ncol = 3 * K)
  CS_UB_min   <- matrix(0, nrow = M, ncol = 3 * K)
  CS_2B_min   <- matrix(0, nrow = M, ncol = 3 * K)
  fail_1_min  <- matrix(0, nrow = M, ncol = 3 * K)
  fail_ub_min <- matrix(0, nrow = M, ncol = 3 * K)
  fail_2b_min <- matrix(0, nrow = M, ncol = 3 * K)
  
  for (m in 1:M) {
    draw    <- randi(length(P0), N, 1)
    p0      <- P0[draw]
    p1      <- P1[draw]
    p1min   <- P1min[draw]
    
    d2              <- runif(N, min = 0, max = 1)
    pub_bias_draw   <- runif(N, min = 0, max = 1)
    pub_bias_sharp  <- 1 * (pub_bias_draw < 0.1)
    
    for (k in 1:K) {
      tau     <- Tau[k]
      d3      <- 1 * (d2 < tau)
      P       <- d3 * p1 + (1 - d3) * p0
      P_min   <- d3 * p1min + (1 - d3) * p0
      
      # Publication bias corrections
      publish_sharp        <- 1 * (P <= 0.05) + (P > 0.05) * pub_bias_sharp
      P_published_sharp    <- P[publish_sharp == 1]
      P_published_smooth   <- P[pub_bias_draw < exp(-8.45 * P)]
      
      publish_sharp_min      <- 1 * (P_min <= 0.05) + (P_min > 0.05) * pub_bias_sharp
      P_published_sharp_min  <- P_min[publish_sharp_min == 1]
      P_published_smooth_min <- P_min[pub_bias_draw < exp(-8.45 * P_min)]
      
      for (pb_mode in 1:3) {
        if ((pb_mode == 1) || (tau == 1) || (tau == 0.5) || (tau == 0)) {
          if (pb_mode == 2) {
            P     <- P_published_sharp
            P_min <- P_published_sharp_min
          }
          if (pb_mode == 3) {
            P     <- P_published_smooth
            P_min <- P_published_smooth_min
          }
          
          col <- k + (pb_mode - 1) * K
          
          LocBin[m, col]      <- LocBin[m, col]      + 1 * (Binomial(P, 0.04, 0.05, "c") <= alpha)
          Disc[m, col]        <- Disc[m, col]        + 1 * (Discontinuity_test(P, 0.05) <= alpha)
          LCM8[m, col]        <- LCM8[m, col]        + 1 * (LCM(P, 0, p_max) <= alpha)
          FM[m, col]          <- FM[m, col]          + 1 * (Fisher(P, 0, p_max) <= alpha)
          cs1                 <- CoxShi(P, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
          CS_1[m, col]        <- CS_1[m, col]        + 1 * (cs1 <= alpha)
          csub                <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0)
          CS_UB[m, col]       <- CS_UB[m, col]       + 1 * (csub <= alpha)
          cs2b                <- CoxShi(P, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
          CS_2B[m, col]       <- CS_2B[m, col]       + 1 * (cs2b <= alpha)
          fail_1[m, col]      <- 1 * (cs1 == 999)  + 1 * (cs1 == 888)
          fail_ub[m, col]     <- 1 * (csub == 999) + 1 * (csub == 888)
          fail_2b[m, col]     <- 1 * (cs2b == 999) + 1 * (cs2b == 888)
          
          LocBin_min[m, col]  <- LocBin_min[m, col]  + 1 * (Binomial(P_min, 0.04, 0.05, "c") <= alpha)
          Disc_min[m, col]    <- Disc_min[m, col]    + 1 * (Discontinuity_test(P_min, 0.05) <= alpha)
          LCM8_min[m, col]    <- LCM8_min[m, col]    + 1 * (LCM(P_min, 0, p_max) <= alpha)
          FM_min[m, col]      <- FM_min[m, col]      + 1 * (Fisher(P_min, 0, p_max) <= alpha)
          cs1_min             <- CoxShi(P_min, 0, 0, p_max, bins, 1, 0, 0, 0, 0, 1)
          CS_1_min[m, col]    <- CS_1_min[m, col]    + 1 * (cs1_min <= alpha)
          csub_min            <- CoxShi(P_min, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0)
          CS_UB_min[m, col]   <- CS_UB_min[m, col]   + 1 * (csub_min <= alpha)
          cs2b_min            <- CoxShi(P_min, 0, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1)
          CS_2B_min[m, col]   <- CS_2B_min[m, col]   + 1 * (cs2b_min <= alpha)
          fail_1_min[m, col]  <- 1 * (cs1_min == 999)  + 1 * (cs1_min == 888)
          fail_ub_min[m, col] <- 1 * (csub_min == 999) + 1 * (csub_min == 888)
          fail_2b_min[m, col] <- 1 * (cs2b_min == 999) + 1 * (cs2b_min == 888)
        }
      }
    }
  }
  
  K <- K * 3
  
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
  
  FirstCol    <- as.vector(outer(Tau, c("No Pub Bias", "Sharp Pub Bias", "Smooth Pub Bias"), function(tau, case) {
    paste0("tau = ", tau, " / ", case)}))
  concat <- cbind(FirstCol, LocBin, FM, Disc, CS_1, CS_UB, CS_2B, LCM8, fail_1, fail_ub, fail_2b,
                  LocBin_min, FM_min, Disc_min, CS_1_min, CS_UB_min, CS_2B_min, LCM8_min,
                  fail_1_min, fail_ub_min, fail_2b_min)
  
  #concat <- rbind(concat[1:21, ], concat[21 + c(1, 11, 21), ], concat[42 + c(1, 11, 21), ])
  patterns <- c("/ No Pub Bias", "tau = 0 / Sharp", "tau = 0.5 / Sharp", "tau = 1 / Sharp", "tau = 0 / Smooth", "tau = 0.5 / Smooth", "tau = 1 / Smooth")  
  combined_pattern <- paste(patterns, collapse = "|") # "Figure_1|Figure_2|Figure_10"
  filtered_concat <- concat[grepl(combined_pattern, concat[,1]), ]
  return(matrix(filtered_concat, size(filtered_concat)[1], 21))
}
