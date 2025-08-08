##############################################################################
# Work related to Brodeur et al. (2022) "Methods Matter" dataset
# This code calibrates Pi to Gamma(a,b) distribution based on the RCT subset of the data
# and replicates Empirical Application of the paper
# Paper: The Power of Tests for Detecting p-Hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
##############################################################################
# Required packages
# install.packages(c("NlcOptim", "fdrtool", "pracma", "gdata", "spatstat", "rddensity", "stringr", "haven", "stats"))
##############################################################################
# Set Working Directory to Script Location
#
# Ensures all relative paths in the script work as expected.
# WARNING:
# - Works when run with source("scriptname.R") or in RStudio's "Source" button.
# - May NOT work in interactive sessions or R Markdown.
# - If this fails, manually set your working directory.
##############################################################################
setwd(dirname(normalizePath(sys.frame(1)$ofile)))

##############################################################################
# Basic Environment Setup
##############################################################################
rm(list = ls())
# Load libraries
library("fdrtool")
library("pracma")
library("gdata")
library("spatstat")
library("rddensity")
library("foreach")
library("doParallel")
library("stringr")
library(stats)
library(haven)
library(stats)

# Helper: Get decimal places in string/number
decimal_places <- function(x) {
  ifelse(grepl("\\.", x), nchar(sub(".*\\.", "", x)), 0)
}
##############################################################################
# Load Methods Matter Data
##############################################################################
df <- read_dta("MethodsMatterData/MM Data.dta")

# Log-likelihood for convolution of N(0,1) with Gamma(a,b)
GammaLogLikelihood <- function(theta, X) {
  a <- theta[1]
  b <- theta[2]
  
  # Return large negative value if parameters are invalid
  if (a <= 0 || b <= 0) {
    return(1e100)
  }
  
  l <- sapply(X, function(xj) {
    integrand <- function(h) {
      (dnorm(xj - h) + dnorm(xj + h)) * dgamma(h, shape = a, scale = b)
    }
    val <- tryCatch(
      integrate(integrand, lower = 0, upper = Inf)$value,
      error = function(e) 0
    )
    return(val)
  })
  
  if (any(l == 0)) {
    return(1e100)
  } else {
    return(-sum(log(l)))
  }
}


# Optimize
res <- optim(
  par = c(1, 1),
  fn = GammaLogLikelihood,
  X = df$t[df$method == "RCT" & df$t<10], # t-stats for RCT studies and under 10
  method = "Nelder-Mead",
  lower = c(0.2, 0.2) 
)

# Display results
print(paste0("Pi_hat = Gamma(", res[1]$par[1],", ", res[1]$par[2], ")"))

##############################################################################
question <- readline("Would you like to run further to replicate Empirical Application? (Yes/No): ")

if (tolower(question) %in% c("yes", "y")) {
  cat("Running the Empirical Application replication...\n")
  # Your code here
} else {
  stop("Empirical Application replication skipped by user.")
}

##############################################################################
# Empirical Illustration
##############################################################################
report <- df$report

# Extract decimals for each reported variable
decimal_places_t   <- decimal_places(df$t)
decimal_places_pv  <- decimal_places(df$pv)
decimal_places_mu  <- decimal_places(df$mu_orig)
decimal_places_se  <- decimal_places(df$sd_orig)

# Numeric versions (NAs introduced if not numeric)
t_num  <- as.numeric(df$t)
pv_num <- as.numeric(df$pv)
mu_num <- as.numeric(df$mu_orig)
se_num <- as.numeric(df$sd_orig)

# Compute p-values from t-statistics
p_orig <- 2 * (1 - pnorm(abs(t_num)))

# Impute missing (NA) values with 999
pv_num[is.na(pv_num)]   <- 999
mu_num[is.na(mu_num)]   <- 999
se_num[is.na(se_num)]   <- 999

##############################################################################
# Clean/Extract F-statistic
##############################################################################
# Remove non-numeric fstat, replace with ""
df$fstat[!grepl("^\\d+(\\.\\d+)?$", df$fstat)] <- ""
Fstat <- as.numeric(df$fstat)
Fstat <- as.numeric(gsub("[^0-9\\.]", "", df$fstat)) # Ensures only numbers and dots remain

##############################################################################
# More Identifiers & Meta Data
##############################################################################
ID <- as.numeric(factor(paste0(df$journal, df$article)))
N <- length(ID)
method <- as.numeric(factor(df$method, levels = c("DID", "IV", "RCT", "RDD"), labels = c(1,2,3,4)))
journal <- df$journal
year <- as.numeric(df$year)
top5 <- journal %in% c("American Economic Review", "Journal of Political Economy",
                       "Quarterly Journal of Economics", "Econometrica", "Review of Economic Studies")
AJQ <- journal %in% c("American Economic Review", "Journal of Political Economy", "Quarterly Journal of Economics")

##############################################################################
# Source External Functions and Set Parameters
##############################################################################
source("Functions/MC_Tests.R")
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall = k))
p_max <- 0.15
bins <- 15
# Bounds for 2-sided tests
B0_2 <- Bound0(0, p_max, bins)
B1_2 <- Bound1(0, p_max, bins)
B2_2 <- Bound2(0, p_max, bins)

##############################################################################
# Load "Star Wars" Data and Compute IDs/p-values
##############################################################################
df_sw <- read_dta("MethodsMatterData/Star Wars Data.dta")
report_sw <- matrix(NA, nrow = nrow(df_sw), ncol = 1)
p_orig_sw <- matrix(NA, nrow = nrow(df_sw), ncol = 1)
ID_sw <- as.numeric(factor(paste0(df_sw$journal, df_sw$article_page, df_sw$year, df_sw$issue)))
N_sw <- length(ID_sw)

for (j in 1:nrow(df_sw)) {
  if (df_sw$t_stat[j] != "") {
    report_sw[j,1] <- "t"
    p_orig_sw[j,1] <- 2 * (1 - pnorm(abs(as.numeric(df_sw$t_stat[j]))))
  } else if (df_sw$p_value[j] != "") {
    report_sw[j,1] <- "p"
    p_orig_sw[j,1] <- as.numeric(df_sw$p_value[j])
  } else if (df_sw$coefficient[j] != "" & df_sw$standard_deviation[j] != "") {
    report_sw[j,1] <- "s"
    p_orig_sw[j,1] <- 2 * (1 - pnorm(abs(as.numeric(df_sw$coefficient[j]) / as.numeric(df_sw$standard_deviation[j]))))
  }
}

# Decimal places for SW
decimal_places_t_sw   <- decimal_places(df_sw$t)
decimal_places_pv_sw  <- decimal_places(df_sw$p_value)
decimal_places_mu_sw  <- decimal_places(df_sw$coefficient)
decimal_places_se_sw  <- decimal_places(df_sw$standard_deviation)

t_num_sw  <- as.numeric(df_sw$t)
pv_num_sw <- as.numeric(df_sw$p_value)
mu_num_sw <- as.numeric(df_sw$coefficient)
se_num_sw <- as.numeric(df_sw$standard_deviation)

t_num_sw[is.na(t_num_sw)]   <- 999
pv_num_sw[is.na(pv_num_sw)] <- 999
mu_num_sw[is.na(mu_num_sw)] <- 999
se_num_sw[is.na(se_num_sw)] <- 999

##############################################################################
# MAIN MONTE CARLO LOOP: DE-ROUNDING AND TESTING
##############################################################################
Results_avg <- matrix(0, nrow = 7, ncol = 14)
Results_orig <- matrix(0, nrow = 7, ncol = 14)
M <- 1000 # Number of MC repetitions

for (m in 1:M) {
  cat("Iteration:", m, "\n")
  
  # De-rounding: Add uniform noise to reported values at the implied decimal precision
  t_dr   <- abs(t_num + (runif(N) - 0.5) * 10^(-decimal_places_t))
  pv_dr  <- abs(pv_num + (runif(N) - 0.5) * 10^(-decimal_places_pv))
  mu_dr  <- mu_num + (runif(N) - 0.5) * 10^(-decimal_places_mu)
  se_dr  <- abs(se_num + (runif(N) - 0.5) * 10^(-decimal_places_se))
  
  # Reconstruct p-values for each type of report
  p_t  <- 2 * pnorm(-t_dr)
  p_p  <- pv_dr
  p_s  <- 2 * pnorm(-abs(mu_dr / se_dr))
  p_ci <- 2 * pnorm(-abs(mu_dr / se_num))
  
  p <- p_t * (report == "t") + p_p * (report == "p") + p_s * (report == "s") + p_ci * (report == "ci")
  p[p > 1] <- 1
  
  # De-rounding for Star Wars
  t_dr_sw   <- abs(t_num_sw + (runif(N_sw) - 0.5) * 10^(-decimal_places_t_sw))
  pv_dr_sw  <- abs(pv_num_sw + (runif(N_sw) - 0.5) * 10^(-decimal_places_pv_sw))
  mu_dr_sw  <- mu_num_sw + (runif(N_sw) - 0.5) * 10^(-decimal_places_mu_sw)
  se_dr_sw  <- abs(se_num_sw + (runif(N_sw) - 0.5) * 10^(-decimal_places_se_sw))
  
  p_t_sw  <- 2 * pnorm(-t_dr_sw)
  p_p_sw  <- pv_dr_sw
  p_s_sw  <- 2 * pnorm(-abs(mu_dr_sw / se_dr_sw))
  p_ci_sw <- 2 * pnorm(-abs(mu_dr_sw / se_num_sw))
  
  p_sw <- p_t_sw * (report_sw == "t") + p_p_sw * (report_sw == "p") +
    p_s_sw * (report_sw == "s") + p_ci_sw * (report_sw == "ci")
  p_sw[p_sw > 1] <- 1
  
  # -- Main loop over "groups"/subsamples --
  for (j in 1:13) {
    # Select observations by method/sample
    if (j == 1)      { idx <- method == 1 }                    # DID
    else if (j == 2) { idx <- method == 3 }                    # RCT
    else if (j == 3) { idx <- method == 4 }                    # RDD
    else if (j == 4) { idx <- method == 2 }                    # IV
    else if (j == 5) { idx <- (Fstat < 30) & (method == 2) & !is.na(Fstat) }   # IV F<30
    else if (j == 6) { idx <- (Fstat >= 30) & (method == 2) & !is.na(Fstat) }  # IV F>=30
    else if (j == 7) { idx <- top5 }                           # Top 5 journals
    else if (j == 8) { idx <- !top5 }                          # Non-top 5
    else if (j == 9) { idx <- year == 2015 }
    else if (j == 10){ idx <- year == 2018 }
    else if (j == 11){ idx <- AJQ }
    else if (j == 12){ idx <- !is.na(report_sw) & !is.na(p_orig_sw)
    id  <- ID_sw[idx]; pval <- p_sw[idx]
    } # Star Wars, special
    else if (j == 13){ idx <- rep(TRUE, N) }                   # Full sample
    
    if (j != 12) {
      id <- ID[idx]
      pval <- p[idx]
    }
    # For original results, use last MC draw
    if (m == M) {
      if (j == 12) {
        id_orig <- ID_sw[idx]
        pval_orig <- p_orig_sw[idx]
      } else {
        id_orig <- ID[idx]
        pval_orig <- p_orig[idx]
      }
    }
    
    # Run the battery of tests
    Results_avg[1, j] <- Results_avg[1, j] + (Binomial(pval, 0.04, 0.05, "c") < 0.05) / M
    Results_avg[2, j] <- Results_avg[2, j] + (Discontinuity_test(pval, 0.05) < 0.05) / M
    Results_avg[3, j] <- Results_avg[3, j] + (CoxShi(pval, id, 0, p_max, bins, 1, 0, 0, 0, 0, 1) < 0.05) / M
    Results_avg[4, j] <- Results_avg[4, j] + (CoxShi(pval, id, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0) < 0.05) / M
    Results_avg[5, j] <- Results_avg[5, j] + (CoxShi(pval, id, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1) < 0.05) / M
    Results_avg[6, j] <- Results_avg[6, j] + (LCM(pval, 0, p_max) < 0.05) / M
    Results_avg[7, j] <- length(pval)
    
    # For original (un-derounded) p-values, save just once
    if (m == M) {
      Results_orig[1, j] <- specify_decimal(Binomial(pval_orig, 0.04, 0.05, "c"), 3)
      Results_orig[2, j] <- specify_decimal(Discontinuity_test(pval_orig, 0.05), 3)
      Results_orig[3, j] <- specify_decimal(CoxShi(pval_orig, id_orig, 0, p_max, bins, 1, 0, 0, 0, 0, 1), 3)
      Results_orig[4, j] <- specify_decimal(CoxShi(pval_orig, id_orig, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 0), 3)
      Results_orig[5, j] <- specify_decimal(CoxShi(pval_orig, id_orig, 0, p_max, bins, 2, 1, B0_2, B1_2, B2_2, 1), 3)
      Results_orig[6, j] <- specify_decimal(LCM(pval_orig, 0, p_max), 3)
      Results_orig[7, j] <- as.character(length(pval_orig))
    }
  }
}

##############################################################################
# Overall Rejection Rates
##############################################################################

Results_avg[1:7, 14] = rowMeans(Results_avg[1:7, 1:13])
Results_orig[1:7, 14] = rowMeans(Results_orig[1:7, 1:13]<=0.05)

Results_avg[7, 14] <- " "
Results_orig[7, 14] <- " "
##############################################################################
# Fill Results Table Labels
##############################################################################
test_names <- c("LocBin", "Discont.", "CS1", "CSUB", "CS2B", "LCM", "Obs")
subsample_names <- c("DID", "RCT", "RDD", "IV", "F<30", "F>=30", 
               "Top 5", "!Top 5", "2015", "2018", "AJQ", "Star Wars", "All", "Overall RR")

colnames(Results_orig) <- subsample_names
rownames(Results_orig) <- test_names

colnames(Results_avg) <- subsample_names
rownames(Results_avg) <- test_names
##############################################################################
# Write Output Files
##############################################################################
if (!dir.exists("csvFiles")) dir.create("csvFiles")
if (!dir.exists("csvFiles/Empirical_Application")) dir.create("csvFiles/Empirical_Application")
write.csv(Results_orig, file = "csvFiles/Empirical_Application/Table3a.csv",  row.names = TRUE, col.names = TRUE)
write.csv(Results_avg,      file = "csvFiles/Empirical_Application/Table3b.csv", row.names = TRUE, col.names = TRUE)
