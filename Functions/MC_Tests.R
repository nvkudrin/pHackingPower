#####################################################################
# Tests for Monte Carlo Simulation
# Paper: Power of Tests for Detecting p-Hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
#####################################################################

# Required packages
# install.packages(c("NlcOptim", "fdrtool", "pracma", "gdata", "spatstat", "rddensity"))


library("fdrtool")
library("pracma")
library("gdata")
library("spatstat")
library("rddensity")
library("ggplot2")
library("matrixcalc")

# Helper: round and format numeric output to k decimal places
specify_decimal <- function(x, k) {
  trimws(format(round(x, k), nsmall = k))
}

set.seed(123)

#####################################################################
# Binomial Test
# type = "c": test on closed interval [p_min, p_max]
# type = "o": test on open interval   (p_min, p_max)
#####################################################################

Binomial <- function(P, p_min, p_max, type) {
  if (type == "c") {
    P <- P[P <= p_max & P >= p_min]
  }
  if (type == "o") {
    P <- P[P < p_max & P > p_min]
  }
  nn <- length(P)
  kk <- sum(P > (p_max + p_min) / 2)
  return(1 - pbinom(kk - 1, nn, 0.5))
}

#####################################################################
# LCM Test
#####################################################################

# Simulate Brownian Bridge (BB) paths and computes ||LCM(BB) - BB||
SimBB <- function(M) {
  N <- 10000
  x <- 1:N / N
  BBsup <- matrix(0, nrow = M, ncol = 1)
  
  for (m in 1:M) {
    eps <- rnorm(n = N)
    eps <- eps / sqrt(N)
    W <- cumsum(eps)
    B <- W - x * W[N]  # Brownian bridge
    C <- c(0, x)
    B <- c(0, B)
    
    lcmaj <- gcmlcm(C, B, type = "lcm")
    f_lcm <- approxfun(lcmaj$x.knots, lcmaj$y.knots)
    y <- f_lcm(C)
    
    BBsup[m, 1] <- max(abs(y - B))
  }
  
  return(BBsup)
}

# Empirical CDF of sup-norm deviations under null
F_LCMsup <- ecdf(SimBB(10000))

# LCM test on [p_min, p_max]; P is a vector of p-values
LCM <- function(P, p_min, p_max) {
  P <- P[P <= p_max & P >= p_min]
  nn <- length(P)
  
  f <- ecdf(P)
  x <- seq(0, 1, length = 1000)
  y <- f(x * (p_max - p_min) + p_min)
  
  lcmaj <- gcmlcm(x, y, type = "lcm")
  ff_lcm <- approxfun(lcmaj$x.knots, lcmaj$y.knots)
  z <- as.numeric(ff_lcm(x))
  
  Test_stat <- sqrt(nn) * max(abs(y - z))
  return(1 - F_LCMsup(Test_stat))
}

#####################################################################
# Fisher's Method for combining p-values
# Tests for excess mass near lower bound
#####################################################################

Fisher <- function(P, p_min, p_max) {
  P <- P[P < p_max & P >= p_min]
  nn <- length(P)
  statFM <- -2 * sum(log(1 - (P - p_min) / (p_max - p_min)))
  return(1 - pchisq(statFM, df = 2 * nn))
}

#####################################################################
# Discontinuity Test using rddensity package
#####################################################################

Discontinuity_test <- function(P, c) {
  res <- rddensity(P, c = c)
  return(res$test$p_jk)
}

#####################################################################
# Cox & Shi Test Helper Functions
# Used to define bounds under null (monotonicity or k-monotonicity)
#####################################################################

# Two-sided critical value function for a given p
cv2 <- function(p) {
  qnorm(1 - p / 2)
}

# Lambda function for computing theoretical bounds
lambda2 <- function(x1, x2, h) {
  pnorm(cv2(x1) - h) - pnorm(cv2(x2) - h) +
    pnorm(cv2(x1) + h) - pnorm(cv2(x2) + h)
}

#####################################################################
# Bounds for Two-Sided t-tests
# These are used when applying monotonicity tests on [0, p_max]
#####################################################################

# Bound on raw bin probabilities
Bound0 <- function(p_min, p_max, J) {
  h <- seq(0, 100, by = 0.001)
  X <- linspace(p_min, p_max, J + 1)
  B <- matrix(0, J, 1)
  
  for (j in 1:J) {
    Obj1 <- lambda2(X[j], X[j + 1], h)
    B[j] <- max(Obj1)
  }
  
  if (p_min == 0) {
    B[1] <- 1
  }
  
  return(B)
}

# Bound on first differences of bin probabilities
Bound1 <- function(p_min, p_max, J) {
  h <- seq(0, 100, by = 0.001)
  X <- linspace(p_min, p_max, J + 1)
  B <- matrix(0, J - 1, 1)
  
  for (j in 1:(J - 1)) {
    Obj1 <- lambda2(X[j], X[j + 1], h)
    Obj2 <- lambda2(X[j + 1], X[j + 2], h)
    A <- Obj2 - Obj1
    B[j] <- max(abs(A))
  }
  
  if (p_min == 0) {
    B[1] <- 1
  }
  
  return(B)
}

# Bound on second differences of bin probabilities
Bound2 <- function(p_min, p_max, J) {
  h <- seq(0, 100, by = 0.001)
  X <- linspace(p_min, p_max, J + 1)
  B <- matrix(0, J - 2, 1)
  
  for (j in 1:(J - 2)) {
    Obj1 <- lambda2(X[j], X[j + 1], h)
    Obj2 <- lambda2(X[j + 1], X[j + 2], h)
    Obj3 <- lambda2(X[j + 2], X[j + 3], h)
    A <- Obj3 - 2 * Obj2 + Obj1
    B[j] <- max(abs(A))
  }
  
  if (p_min == 0) {
    B[1] <- 1
  }
  
  return(B)
}

#####################################################################
# Bounds for One-Sided t-tests
# These are used when applying monotonicity tests on [0, p_max]
#####################################################################

# Bound on bin probabilities (0th difference)
Bound0_1s <- function(pmax, J) {
  X <- linspace(0, pmax, J + 1)
  B <- matrix(0, J, 1)
  for (j in 1:J) {
    B[j] <- 2 * pnorm((qnorm(1 - X[j]) - qnorm(1 - X[j + 1])) / 2) - 1
  }
  B[1] <- 1
  return(B)
}

# Bound on first differences
Bound1_1s <- function(pmax, J) {
  h <- linspace(0, 100, n = 100000)
  X <- linspace(0, pmax, J + 1)
  B <- matrix(0, J - 1, 1)
  
  for (j in 1:(J - 1)) {
    Obj1 <- pnorm(qnorm(1 - X[j]) - h) - pnorm(qnorm(1 - X[j + 1]) - h)
    Obj2 <- pnorm(qnorm(1 - X[j + 1]) - h) - pnorm(qnorm(1 - X[j + 2]) - h)
    A <- Obj2 - Obj1
    B[j] <- max(abs(A))
  }
  
  B[1] <- 1
  return(B)
}

# Bound on second differences
Bound2_1s <- function(pmax, J) {
  h <- linspace(0, 100, n = 100000)
  X <- linspace(0, pmax, J + 1)
  B <- matrix(0, J - 2, 1)
  
  for (j in 1:(J - 2)) {
    Obj1 <- pnorm(qnorm(1 - X[j]) - h) - pnorm(qnorm(1 - X[j + 1]) - h)
    Obj2 <- pnorm(qnorm(1 - X[j + 1]) - h) - pnorm(qnorm(1 - X[j + 2]) - h)
    A1 <- Obj2 - Obj1
    
    Obj3 <- pnorm(qnorm(1 - X[j + 1]) - h) - pnorm(qnorm(1 - X[j + 2]) - h)
    Obj4 <- pnorm(qnorm(1 - X[j + 2]) - h) - pnorm(qnorm(1 - X[j + 3]) - h)
    A2 <- Obj4 - Obj3
    
    B[j] <- max(abs(A2 - A1))
  }
  
  B[1] <- 1
  return(B)
}

#####################################################################
# Cox & Shi Test for Monotonicity and Higher-Order Shape Constraints
# Inputs:
# - Q: vector of p-values
# - ind: vector of paper IDs
# - J: number of histogram bins
# - K: order of shape constraint (e.g., 1 = monotone, 2 = convex)
# - B: use bounds (1 = yes, 0 = no)
# - Bnd0, Bnd1, Bnd2: theoretical bounds
# - Include_LB: whether to include lower bounds (0 or 1)
#####################################################################

CoxShi <- function(Q, ind, p_min, p_max, J, K, B, Bnd0, Bnd1, Bnd2, Include_LB) {
  B0 <- Bnd0
  B1 <- Bnd1
  B2 <- Bnd2
  
  P <- Q[Q <= p_max & Q >= p_min]
  
  if (length(ind) > 1) {
    ind <- ind[Q <= p_max & Q >= p_min]
    indu <- unique(ind)
  }
  
  Bnd_adj <- length(P) / length(Q)
  N <- length(P)
  bin <- seq(p_min, p_max, length = J + 1)
  Phat <- matrix(0, nrow = J - 1, ncol = 1)
  
  # Compute histogram bin frequencies
  for (s in 1:(J - 1)) {
    Phat[s] <- sum((P > bin[s]) * (P <= bin[s + 1])) / N
  }
  Phat[1] <- Phat[1] + sum(P == bin[1]) / N
  
  # Adjust bounds
  if (B == 0) {
    B0 <- -matrix(1, nrow = J, ncol = 1)
  }
  
  if (B == 1) {
    B0 <- -B0 / Bnd_adj
    B1 <- -B1 / Bnd_adj
    B2 <- -B2 / Bnd_adj
    if (p_min == 0) {
      B0[1] <- -1
      B1[1] <- -1
      B2[1] <- -1
    }
  }
  
  # Compute covariance matrix Omega
  if (length(ind) > 1) {
    Omega <- matrix(0, J - 1, J - 1)
    for (i in indu) {
      X <- P[ind == i]
      mq <- repmat(matrix(X, 1, length(X)), J - 1, 1)
      a1 <- (mq <= bin[2:J])
      a2 <- (mq > bin[1:(J - 1)])
      mq0 <- 1 * (mq == 0)
      mq0[2:(J - 1), ] <- 0
      mq <- 1 * (a1 * a2) + mq0
      mq <- mq - repmat(Phat, 1, length(X))
      Omega <- Omega + mq %*% matrix(1, length(X), length(X)) %*% t(mq)
    }
    Omega <- Omega / length(P)
  }
  
  if (length(ind) == 1) {
    if (min(Phat) == 0) {
      Qhat <- Phat * N / (N + 1) + 1 / (J * (N + 1))
      #Qhat <- Phat * (N - J) / (N) + 1 / (N )
      Omega <- diag(c(Qhat)) - Qhat %*% t(Qhat)
    }
    if (min(Phat) > 0) {
      Omega <- diag(c(Phat)) - Phat %*% t(Phat)
    }
  }
  
  # Construct first-order difference matrix D
  D <- matrix(0, J - 1, J)
  for (i in 1:(J - 1)) {
    D[i, i]     <- -1
    D[i, i + 1] <- 1
  }
  
  Dk <- -D
  
  # Construct higher-order difference matrix (up to order K)
  if (K > 1) {
    d <- D
    for (k in 2:K) {
      d <- D[1:(J - k), 1:(J - k + 1)] %*% d
      Dk <- rbind(Dk, (-1)^k * d)
    }
  }
  
  # Constraint matrix depending on bounds and lower bound inclusion
  if (B == 0) {
    Dk <- rbind(-diag(J), diag(J), Dk)
  }
  if (B == 1 && Include_LB == 1) {
    Dk <- rbind(-diag(J), -Dk, diag(J), Dk)
  }
  if (B == 1 && Include_LB == 0) {
    Dk <- rbind(-diag(J), -Dk, diag(J))
  }
  if (K == 0) {
    Dk <- rbind(-diag(J), diag(J))
  }
  
  # Construct target vectors
  eJ <- matrix(0, J, 1)
  eJ[J] <- 1
  
  F1 <- rbind(-diag(J - 1), matrix(1, 1, J - 1))
  
  # Right-hand side vector c
  c <- matrix(0, (K + 1) * (J - K / 2), 1)
  
  if (Include_LB == 0) {
    c <- matrix(0, J, 1)
  }
  
  if (B == 0) {
    c <- rbind(matrix(-1, J, 1), c)
  }
  
  if (B == 1) {
    if (K == 0) {
      c <- rbind(B0, c)
    }
    if (K == 1) {
      c <- rbind(B0, B1, c)
    }
    if (K == 2) {
      c <- rbind(B0, B1, B2, c)
    }
  }
  
  # Final constraint system: A %*% t <= b
  A <- Dk %*% F1
  b <- Dk %*% eJ - c
  
  ###################################################################
  # Optimization
  # QLR
  # under A %*% t <= b
  ###################################################################
  
  # Invert Omega safely
  myQ <- tryCatch(solve(Omega), error = function(e) "error")
  
  if (myQ[1] != "error") {
    fn <- function(t) {
      diff <- Phat - t
      return(N * t(diff) %*% myQ %*% diff)
    }
    
    t0 <- matrix(1, J - 1, 1) / J
    res <- tryCatch(fmincon(t0, fn, A = A, b = b), error = function(e) "error")
    iter <- 0
    
    while ((res[1] == "error") && (iter < 10)) {
      iter <- iter + 1
      ru <- runif(J)
      t0 <- matrix(ru[1:(J - 1)] / sum(ru), J - 1, 1)
      res <- tryCatch(fmincon(t0, fn, A = A, b = b), error = function(e) "error")
    }
    
    if (iter == 10) {
      return(999)  # Optimization failed after 10 attempts
    } else {
      t_opt <- t(t(res$par))
      T <- fn(t_opt)
      
      # Degrees of freedom = rank of active constraints
      Ba <- A[which(res$info$lambda$ineqlin > 0), ]
      JX <- qr(Ba)$rank
      
      if (res$convergence == 0) {
        return(1 - pchisq(T, df = JX) * (JX > 0))
      } else {
        return(999)
      }
    }
  } else {
    return(888)  # Omega not invertible
  }
}
