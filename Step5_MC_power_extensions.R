#####################################################################
# Monte Carlo simulation: Extensions (Appendices F-I)
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

####################### Examine the effect of changing rho (Covariate Selection Analytical example) on power #################################
print("Appendix F")
tic()
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

# Define header row
header <- matrix(c(
  "rho", "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b",
  "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b"
), 1, 21)

rho = seq(0.1, 1, by=0.1)
finalMatrix <- foreach(i=0:3, .combine=cbind) %dopar% {
  source("Functions/MC_Tests.R")
  source("Functions/MC_power_extensions.R")

  tempMatrix = MC_examine_rho(rho, i, 5000, 0.15, 15, 1000, 1, 0.5)

  header_top1 <- matrix(paste(i, "Threshold"), 1, 10)
  header_top2 <- matrix(paste(i, "Minimum"), 1, 10)
  header_top <- cbind("rho", header_top1, header_top2)
  
  # Combine everything
  tempMatrix <- rbind(header_top, header, tempMatrix)
  tempMatrix
  
  
}

stopCluster(cl)
toc()
# Extract column names from first row (excluding first element)
column_names <- finalMatrix[1, ]  # Skip first element, take rest of first row
# Extract the actual data (excluding first row and first column)
data_matrix <- finalMatrix[-1, ]
# Set the names
colnames(data_matrix) <- column_names
if (!dir.exists("csvFiles/Power_Calculations")) dir.create("csvFiles/Power_Calculations")
# Write to CSV
write.csv(data_matrix, file = "csvFiles/Power_Calculations/RejectionRates_different_rho.csv", row.names = FALSE, col.names = TRUE)

####################### Examine the effect of collecting non-focal t-statistics (Covariate Selection Analytical example) on power #################################
print("Appendix G")
tic()
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

# Define header row
header <- matrix(c(
  "pmis", "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b",
  "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b"
), 1, 21)

pmis = seq(0, 1, by=0.1)
finalMatrix <- foreach(i=0:3, .combine=cbind) %dopar% {
  source("Functions/MC_Tests.R")
  source("Functions/MC_power_extensions.R")

  tempMatrix = MC_examine_pmis(pmis, i, 5000, 0.15, 15, 1000, 2, 0.5, 0.5)

  header_top1 <- matrix(paste(i, "Threshold"), 1, 10)
  header_top2 <- matrix(paste(i, "Minimum"), 1, 10)
  header_top <- cbind("pmis", header_top1, header_top2)
  
  # Combine everything
  tempMatrix <- rbind(header_top, header, tempMatrix)
  tempMatrix
}
stopCluster(cl)
toc()
# Extract column names from first row (excluding first element)
column_names <- finalMatrix[1, ]  # Skip first element, take rest of first row
# Extract the actual data (excluding first row and first column)
data_matrix <- finalMatrix[-1, ]
# Set the names
colnames(data_matrix) <- column_names

if (!dir.exists("csvFiles/Power_Calculations")) dir.create("csvFiles/Power_Calculations")
# Write to CSV
write.csv(data_matrix, file = "csvFiles/Power_Calculations/RejectionRates_different_pmis.csv", row.names = FALSE, col.names = TRUE)


####################### Examine the effect of sampel size (Covariate Selection Analytical example (2-sided tests)) on power #################################
print("Appendix H")
tic()
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


folder <- "csvFiles/Distributions/CovariateSelection/"
P <- c("Covariate_03_2sided_1.csv","Covariate_13_2sided_1.csv","Covariate_23_2sided_1.csv","Covariate_33_2sided_1.csv")

# Define header row
header <- matrix(c(
  "N", "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b",
  "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b"
), 1, 21)

sample_sizes <- seq(400, 5000, by=400)

finalMatrix <- foreach(i=0:3, .combine=cbind) %dopar% {
    source("Functions/MC_Tests.R")
    source("Functions/MC_power_extensions.R")
    library("stringr")
    
    # Load null distribution
    P0 <- read.csv(file = paste0(folder, "P0_", P[i+1]), header = FALSE)[, 1]
    # Load alternative distribution (Thresholding)
    P1 <- read.csv(file = paste0(folder, "P1_", P[i+1]), header = FALSE)[, 1]
    # Load alternative distribution (Minimum)
    P1min <- read.csv(file = paste0(folder, "P1min_", P[i+1]), header = FALSE)[, 1]
   
    
    # Run Monte Carlo power simulations
    tempMatrix <- MC_examine_n(P0, P1, P1min, sample_sizes, 0.15, 6, 15, 1000, 2, 0.5)
    
    header_top1 <- matrix(paste(i, "Threshold"), 1, 10)
    header_top2 <- matrix(paste(i, "Minimum"), 1, 10)
    header_top <- cbind("N", header_top1, header_top2)
    
    # Combine everything
    tempMatrix <- rbind(header_top, header, tempMatrix)
    tempMatrix
}
stopCluster(cl)
toc()
# Extract column names from first row (excluding first element)
column_names <- finalMatrix[1, ]  # Skip first element, take rest of first row
# Extract the actual data (excluding first row and first column)
data_matrix <- finalMatrix[-1, ]
# Set the names
colnames(data_matrix) <- column_names
if (!dir.exists("csvFiles/Power_Calculations")) dir.create("csvFiles/Power_Calculations")
# Write to CSV
write.csv(data_matrix, file = "csvFiles/Power_Calculations/RejectionRates_different_N.csv", row.names = FALSE, col.names = TRUE)

####################### Examine the relationship between bias and power #################################
print("Appendix I")

# Set up parallel backend
cores <- detectCores()
cl <- makeCluster(cores[1] - 2)  # leave two cores free
registerDoParallel(cl)

# Load file paths for different designs
CovariateSelection    <- list.files("csvFiles/ForAppendixI/Distributions/CovariateSelection", pattern = "P0_", full.names = TRUE)

P <- c(CovariateSelection)

# Define header row
header <- matrix(c("DGP/Test", "LocBin", "FM", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM",
  "fail_cs1", "fail_csub", "fail_cs2b"), 1, 11)

# Start simulation
tic()

finalMatrix <- foreach(i = 1:length(P), .combine = rbind) %dopar% {
  source("Functions/MC_Tests.R")
  source("Functions/MC_power_extensions.R")
  library("stringr")
  
  # Parse file path
  P_parts <- strsplit(P[i], "_")[[1]]
  sided <- strtoi(substr(P_parts[5], 1, 1))
  
  # Load null distribution
  P0 <- read.csv(file = P[i], header = FALSE)[, 1]
  
  # Load alternative distribution (Thresholding)
  str_sub(P_parts[1], -1, -1) <- "1"
  P1 <- read.csv(file = paste(P_parts, collapse = "_"), header = FALSE)[, 1]
  
  # Load alternative distribution (Minimum)
  str_sub(P_parts[1], -1, -1) <- "1min"
  P1min <- read.csv(file = paste(P_parts, collapse = "_"), header = FALSE)[, 1]
  
  # Run Monte Carlo power simulations
  tempMatrix <- MC_for_power_bias(P0, P1, P1min, 5000, 0.15, 15, 5000, sided, 0.25, 0.9)
  
  # Create top headers
  P_parts[6] <- substr(P_parts[6], 1, 1)
  
  P_parts[7] <- "Threshold"
  row1 <- matrix(paste(P_parts[2:7], collapse = "_"), 1, 1)
  
  P_parts[7] <- "Minimum"
  row2 <- matrix(paste(P_parts[2:7], collapse = "_"), 1, 1)
  
  rows <- rbind(row1, row2)
  
  # Combine everything
  tempMatrix <- cbind(rows, tempMatrix)
  
  tempMatrix
}
finalMatrix <- rbind(header, finalMatrix)

toc()
# Shut down parallel cluster
stopCluster(cl)

# Extract column names from first row (excluding first element)
column_names <- finalMatrix[1, -1]  # Skip first element, take rest of first row
# Extract row names from first column (excluding first element)  
row_names <- finalMatrix[-1, 1]     # Skip first element, take rest of first column
# Extract the actual data (excluding first row and first column)
data_matrix <- finalMatrix[-1, -1]
# Set the names
colnames(data_matrix) <- column_names
rownames(data_matrix) <- row_names
if (!dir.exists("csvFiles/Power_Calculations")) dir.create("csvFiles/Power_Calculations")
# Write to CSV
write.csv(data_matrix, file = "csvFiles/Power_Calculations/RejectionRates_many_h.csv", row.names = TRUE, col.names = TRUE)
