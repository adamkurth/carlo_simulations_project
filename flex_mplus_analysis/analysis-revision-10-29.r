# analysis.r
rm(list = ls()) 
library(dplyr) 
library(readr)
library(tidyr)

# -------------------------------------------------------------------
# 0) Working Directory and File Paths
# -------------------------------------------------------------------
cwd <- setwd("/Users/adamkurth/Documents/vscode/code/IRT_MIRT_project")
data.path <- paste0(cwd, "/flex_mplus_analysis/data/", sep = "")

path.flex <- paste0(data.path, "results_fllexmirt_intercepts_slopes.csv")
path.mplus <- paste0(data.path, "results_mplus_intercepts_slopes.csv")
path.mirt <- paste0(data.path, "mirt_processed_MR_copy.csv")

# -------------------------------------------------------------------
# 1) Read Data
# -------------------------------------------------------------------
data.flex  <- read.csv(path.flex,  header = TRUE, sep = ",")
data.mplus <- read.csv(path.mplus, header = TRUE, sep = ",")
data.mirt  <- read.csv(path.mirt,  header = TRUE, sep = ",")

# -------------------------------------------------------------------
# 2) Proper MIRT Data Processing
#    We'll build a 40-row DF: intercept_1..20, slope_1..20
#    then reorder them to interleaved: (intercept_1, slope_1, intercept_2, slope_2, ...)
# -------------------------------------------------------------------
data.mirt.processed <- data.frame(
  MEAN = c(data.mirt$intsemeans, data.mirt$slopesemeans),
  RMSD = c(data.mirt$intRMSSD,  data.mirt$slopeRMSSD),
  VAR  = c(data.mirt$intsevar,  data.mirt$slopesevar)
)

# Temporarily label as intercept_1..20, slope_1..20
temp_rownames <- c(
  paste0("intercept_", 1:20),
  paste0("slope_", 1:20)
)
rownames(data.mirt.processed) <- temp_rownames

# Interleave them: intercept_1, slope_1, intercept_2, slope_2, ...
desired_order <- as.vector(rbind(
  paste0("intercept_", 1:20),
  paste0("slope_", 1:20)
))
data.mirt.processed <- data.mirt.processed[desired_order, ]
rownames(data.mirt.processed) <- desired_order

cat("\n----- MIRT Data (head) to Verify Intercept_i / Slope_i Pairing -----\n")
print(head(data.mirt.processed, 4))

# -------------------------------------------------------------------
# 3) Utility Functions
# -------------------------------------------------------------------
calc.var <- function(x){
  R <- length(x)
  mean.x <- mean(x)
  sum((x - mean.x)^2) / (R - 1)  # sample variance
}
calc.mean <- function(x) mean(x)
calc.rmsd <- function(x) {
  R <- length(x)
  mean.x <- mean(x)
  sqrt(sum((x - mean.x)^2) / R) # RMSD
}

# -------------------------------------------------------------------
# 4) Process data for FlexMIRT and Mplus
#    (both in wide format: columns intercept_1..20, slope_1..20)
# -------------------------------------------------------------------
process.data <- function(data) {
  results <- list()
  for (i in 1:20) {
    # intercept_i
    int_col <- paste0("intercept_", i)
    int_data <- data[[int_col]]
    results[[int_col]] <- c(
      VAR  = calc.var(int_data),
      RMSD = calc.rmsd(int_data),
      MEAN = calc.mean(int_data)
    )
    # slope_i
    slope_col <- paste0("slope_", i)
    slope_data <- data[[slope_col]]
    results[[slope_col]] <- c(
      VAR  = calc.var(slope_data),
      RMSD = calc.rmsd(slope_data),
      MEAN = calc.mean(slope_data)
    )
  }
  out <- as.data.frame(do.call(rbind, results))
  colnames(out) <- c("VAR","RMSD","MEAN")
  
  cat("\nFirst two pairs from processed data:\n")
  print(head(out, 4))
  return(out)
}

data.flex.processed  <- process.data(data.flex)
data.mplus.processed <- process.data(data.mplus)

# -------------------------------------------------------------------
# 5) Write out processed data
# -------------------------------------------------------------------
write.csv(data.flex.processed,  file = paste0(data.path, "flex_processed.csv"),  row.names = TRUE)
write.csv(data.mplus.processed, file = paste0(data.path, "mplus_processed.csv"), row.names = TRUE)
write.csv(data.mirt.processed,  file = paste0(data.path, "mirt_processed.csv"),  row.names = TRUE)

# -------------------------------------------------------------------
# 6) Calculate Ratios
#    We'll build a 40-row data frame in the same row order:
#    row 1 = intercept_1, row 2 = slope_1, row 3 = intercept_2, row 4 = slope_2, ...
# -------------------------------------------------------------------
safe_ratio <- function(numerator, denominator) {
  if (length(numerator) != length(denominator)) {
    warning("Lengths of numerator and denominator do not match")
    return(rep(NA, max(length(numerator), length(denominator))))
  }
  ifelse(denominator == 0, NA, numerator / denominator)
}

# We'll create "Parameter" and "Item" in *interleaved* order:
#   row1 => Intercept, item=1
#   row2 => Slope,     item=1
#   row3 => Intercept, item=2
#   row4 => Slope,     item=2, etc.
param.order <- rep(c("Intercept", "Slope"), 20)
item.order  <- as.vector(rbind(1:20, 1:20))

ratios <- tryCatch({
  data.frame(
    Parameter = param.order,
    Item      = item.order,
    Ratio_MEAN_Flex_MIRT  = safe_ratio(data.flex.processed$MEAN,  data.mirt.processed$MEAN),
    Ratio_MEAN_Flex_Mplus = safe_ratio(data.flex.processed$MEAN,  data.mplus.processed$MEAN),
    Ratio_MEAN_MIRT_Mplus = safe_ratio(data.mirt.processed$MEAN,  data.mplus.processed$MEAN),
    
    Ratio_RMSD_Flex_MIRT  = safe_ratio(data.flex.processed$RMSD,  data.mirt.processed$RMSD),
    Ratio_RMSD_Flex_Mplus = safe_ratio(data.flex.processed$RMSD,  data.mplus.processed$RMSD),
    Ratio_RMSD_MIRT_Mplus = safe_ratio(data.mirt.processed$RMSD,  data.mplus.processed$RMSD)
  )
}, error = function(e) {
  cat("Error in creating ratios data frame:", e$message, "\n")
  return(NULL)
})

# -------------------------------------------------------------------
# 7) Insert True Parameter Values (If Known)
# -------------------------------------------------------------------
true.parameters <- data.frame(
  Item      = 1:20,
  intercept = rep(NA, 20), # or actual intercept values if known
  slope     = rep(NA, 20)  # or actual slope values if known
)

# For each row i in the final table:
#   if row i is an intercept, we use true.parameters$intercept
#   if row i is a slope, we use true.parameters$slope
true_vals <- numeric(40)
for(i in seq_len(40)) {
  if(param.order[i] == "Intercept") {
    true_vals[i] <- true.parameters$intercept[item.order[i]]
  } else {
    true_vals[i] <- true.parameters$slope[item.order[i]]
  }
}

# -------------------------------------------------------------------
# 8) Build final.table (intercept_1, slope_1, intercept_2, slope_2, ... => 40 rows)
# -------------------------------------------------------------------
final.table <- data.frame(
  Parameter = param.order,
  Item      = item.order,
  True_Value = true_vals,
  
  FlexMIRT_Mean = data.flex.processed$MEAN,
  MIRT_Mean     = data.mirt.processed$MEAN,
  Mplus_Mean    = data.mplus.processed$MEAN,
  
  Ratio_MEAN_Flex_MIRT  = ratios$Ratio_MEAN_Flex_MIRT,
  Ratio_MEAN_Flex_Mplus = ratios$Ratio_MEAN_Flex_Mplus,
  Ratio_MEAN_MIRT_Mplus = ratios$Ratio_MEAN_MIRT_Mplus,
  
  FlexMIRT_RMSD = data.flex.processed$RMSD,
  MIRT_RMSD     = data.mirt.processed$RMSD,
  Mplus_RMSD    = data.mplus.processed$RMSD,
  
  FlexMIRT_VAR  = data.flex.processed$VAR,
  MIRT_VAR      = data.mirt.processed$VAR,
  Mplus_VAR     = data.mplus.processed$VAR
)

# NOTE: We do NOT reorder or group_by here. We keep the row order
# that matches data.flex.processed / data.mplus.processed / data.mirt.processed.
# final.table is now 40 rows, interleaved.

# -------------------------------------------------------------------
# 9) Create Separate Tables for Means, RMSD, Variance
# -------------------------------------------------------------------
means.table <- data.frame(
  Parameter  = final.table$Parameter,
  Item       = final.table$Item,
  True_Value = final.table$True_Value,
  FlexMIRT   = final.table$FlexMIRT_Mean,
  MIRT       = final.table$MIRT_Mean,
  Mplus      = final.table$Mplus_Mean,
  Ratio_Flex_MIRT   = final.table$Ratio_MEAN_Flex_MIRT,
  Ratio_Flex_Mplus  = final.table$Ratio_MEAN_Flex_Mplus,
  Ratio_MIRT_Mplus  = final.table$Ratio_MEAN_MIRT_Mplus
)

rmsd.table <- data.frame(
  Parameter  = final.table$Parameter,
  Item       = final.table$Item,
  True_Value = final.table$True_Value,
  FlexMIRT   = final.table$FlexMIRT_RMSD,
  MIRT       = final.table$MIRT_RMSD,
  Mplus      = final.table$Mplus_RMSD,
  Ratio_Flex_MIRT   = final.table$FlexMIRT_RMSD / final.table$MIRT_RMSD,
  Ratio_Flex_Mplus  = final.table$FlexMIRT_RMSD / final.table$Mplus_RMSD,
  Ratio_MIRT_Mplus  = final.table$MIRT_RMSD     / final.table$Mplus_RMSD
)

var.table <- data.frame(
  Parameter  = final.table$Parameter,
  Item       = final.table$Item,
  True_Value = final.table$True_Value,
  FlexMIRT   = final.table$FlexMIRT_VAR,
  MIRT       = final.table$MIRT_VAR,
  Mplus      = final.table$Mplus_VAR,
  Ratio_Flex_MIRT   = final.table$FlexMIRT_VAR / final.table$MIRT_VAR,
  Ratio_Flex_Mplus  = final.table$FlexMIRT_VAR / final.table$Mplus_VAR,
  Ratio_MIRT_Mplus  = final.table$MIRT_VAR     / final.table$Mplus_VAR
)

# -------------------------------------------------------------------
# 10) Write CSV Outputs
# -------------------------------------------------------------------
write.csv(means.table, file = paste0(data.path, "means_table.csv"), row.names = FALSE)
write.csv(rmsd.table,  file = paste0(data.path, "rmsd_table.csv"),  row.names = FALSE)
write.csv(var.table,   file = paste0(data.path, "var_table.csv"),   row.names = FALSE)
write.csv(final.table, file = paste0(data.path, "final_table.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 11) Calculate & Write Overall Ratios
# -------------------------------------------------------------------
calculate_overall_ratios <- function(table) {
  summary.df <- table %>%
    group_by(Parameter) %>%
    summarize(
      Overall_Ratio_Flex_MIRT  = mean(Ratio_Flex_MIRT,  na.rm = TRUE),
      Overall_Ratio_Flex_Mplus = mean(Ratio_Flex_Mplus, na.rm = TRUE),
      Overall_Ratio_MIRT_Mplus = mean(Ratio_MIRT_Mplus, na.rm = TRUE)
    )
  as.data.frame(summary.df)
}

means.overall <- calculate_overall_ratios(means.table)
rmsd.overall  <- calculate_overall_ratios(rmsd.table)
var.overall   <- calculate_overall_ratios(var.table)

write.csv(means.overall, file = paste0(data.path, "means_overall_ratios.csv"),    row.names = FALSE)
write.csv(rmsd.overall,  file = paste0(data.path, "rmsd_overall_ratios.csv"),     row.names = FALSE)
write.csv(var.overall,   file = paste0(data.path, "variance_overall_ratios.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 12) Verification
# -------------------------------------------------------------------
verify_pairing <- function(data, name="Dataset") {
  intercepts <- data[data$Parameter == "Intercept", ]
  slopes     <- data[data$Parameter == "Slope",     ]
  
  # Check same number
  if(nrow(intercepts) != nrow(slopes)) {
    warning(paste(name, ": Mismatch in number of intercepts and slopes"))
    return(FALSE)
  }
  # Check item alignment
  if(!all(intercepts$Item == slopes$Item)) {
    warning(paste(name, ": Mismatch in item numbers between intercepts and slopes"))
    return(FALSE)
  }
  
  cat("\nFirst two pairs from", name, ":\n")
  for(i in 1:2) {
    cat("Row = ", rownames(intercepts)[i], " | Intercept row:\n")
    print(intercepts[i, ])
    cat("Row = ", rownames(slopes)[i], " | Slope row:\n")
    print(slopes[i, ])
  }
  TRUE
}

verify_pairing(means.table,  "Means Table")
verify_pairing(rmsd.table,   "RMSD Table")
verify_pairing(var.table,    "Variance Table")
verify_pairing(final.table,  "Final Table")

# -------------------------------------------------------------------
# 13) Print Sample Output
# -------------------------------------------------------------------
cat("\n----- Means Table (first few rows) -----\n")
print(head(means.table, 10))

cat("\n----- RMSD Table (first few rows) -----\n")
print(head(rmsd.table, 10))

cat("\n----- Variance Table (first few rows) -----\n")
print(head(var.table, 10))

cat("\n----- Overall Ratios: Means -----\n")
print(means.overall)

cat("\n----- Overall Ratios: RMSD -----\n")
print(rmsd.overall)

cat("\n----- Overall Ratios: Variance -----\n")
print(var.overall)

cat("\nScript completed.\n")