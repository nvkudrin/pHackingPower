#####################################################################
# Monte Carlo simulation: main
# Paper: The Power of Tests for Detecting p-Hacking
# Authors: G. Elliott, N. Kudrin, K. Wuthrich
#####################################################################
# Required packages
# install.packages(c("NlcOptim", "fdrtool", "pracma", "gdata", "spatstat", "rddensity", "stringr", "haven", "stats"))
#####################################################################
# Automatically set the working directory to the location of this script.
# This ensures that all relative file paths in the script work as expected, 
# regardless of where the script is sourced from.
# 
# WARNING:
# - This method works when the script is run using source("scriptname.R") or via RStudio's "Source" button.
# - It may NOT work as expected in interactive sessions or R Markdown files.
# - Alternatively, set the working directory to the source file location
setwd(dirname(normalizePath(sys.frame(1)$ofile)))

# Clear environment
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

# Helper function to round to fixed decimal places
specify_decimal <- function(x, k) {
  trimws(format(round(x, k), nsmall = k))
}

# Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores[1] - 2)  # leave two cores free
registerDoParallel(cl)

# Load file paths for different designs
CovariateSelection    <- list.files("csvFiles/Distributions/CovariateSelection", pattern = "P0_", full.names = TRUE)
IVSelection           <- list.files("csvFiles/Distributions/IVSelection", pattern = "P0_", full.names = TRUE)
LagLengthSelection    <- list.files("csvFiles/Distributions/LagLengthSelection", pattern = "P0_", full.names = TRUE)
ClusterSelection      <- list.files("csvFiles/Distributions/ClusterSelection", pattern = "P0_", full.names = TRUE)

P <- c(CovariateSelection, IVSelection, LagLengthSelection, ClusterSelection)

# Define header row
header <- matrix(c(
  "Tau/Pub Bias", "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b",
  "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b"
), 1, 21)

# Start simulation
tic()

Tau = seq(0, 1, by = 0.05) # fraction of p-hackers
finalMatrix <- foreach(i = 1:length(P), .combine = cbind) %dopar% {
  source("Functions/MC_Tests.R")
  source("Functions/MC_power.R")
  library("stringr")
  
  # Parse file path
  P_parts <- strsplit(P[i], "_")[[1]]
  sided <- strtoi(substr(P_parts[4], 1, 1))
  
  # Load null distribution
  P0 <- read.csv(file = P[i], header = FALSE)[, 1]
  
  # Load alternative distribution (Thresholding)
  str_sub(P_parts[1], -1, -1) <- "1"
  P1 <- read.csv(file = paste(P_parts, collapse = "_"), header = FALSE)[, 1]
  
  # Load alternative distribution (Minimum)
  str_sub(P_parts[1], -1, -1) <- "1min"
  P1min <- read.csv(file = paste(P_parts, collapse = "_"), header = FALSE)[, 1]
  
  # Run Monte Carlo power simulations
  tempMatrix <- MC(P0, P1, P1min, 5000, 0.15, 15, 5000, sided, Tau)
  
  # Create top headers
  P_parts[5] <- substr(P_parts[5], 1, 1)
  
  P_parts[6] <- "Threshold"
  header_top1 <- matrix(paste(P_parts[2:6], collapse = "_"), 1, 10)
  
  P_parts[6] <- "Minimum"
  header_top2 <- matrix(paste(P_parts[2:6], collapse = "_"), 1, 10)
  
  header_top <- cbind("", header_top1, header_top2)
  
  # Combine everything
  tempMatrix <- rbind(header_top, header, tempMatrix)
  
  tempMatrix
}

toc()

# Shut down parallel cluster
stopCluster(cl)
# Extract column names from first row (excluding first element)
column_names <- finalMatrix[1, ]  # Skip first element, take rest of first row
# Extract the actual data (excluding first row and first column)
data_matrix <- finalMatrix[-1, ]
# Set the names
colnames(data_matrix) <- column_names
if (!dir.exists("csvFiles/Power_Calculations")) dir.create("csvFiles/Power_Calculations")
# Write to CSV
write.csv(data_matrix, file = "csvFiles/Power_Calculations/RejectionRates_main.csv", row.names = FALSE, col.names = TRUE)
