# =============================================================================
# File: 05_Mantel_analysis_fixed.R
# Description: Fixed Mantel correlation analysis for AMF and Rhizobia with soil properties
# =============================================================================

# Clean workspace
rm(list = ls())

# Load required libraries
library(linkET)
library(dplyr)
library(reshape2)
library(tidyverse)

# Set paths
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

cat("Loading and preparing data...\n")

# Load soil properties
trait <- read.csv(file.path(data_dir, "Mantel_env.csv"), header = TRUE)

# Load microbial data
AMF <- read.csv(file.path(data_dir, "AMF_ASV_rarefaction.csv"), header = TRUE, row.names = 1)
Rhizobia <- read.csv(file.path(data_dir, "Rhizobia_ASV_rarefaction.csv"), header = TRUE, row.names = 1)

# Select samples with N treatments
AMF <- AMF[, grepl("N0|N2|N4|N6", colnames(AMF))]
Rhizobia <- Rhizobia[, grepl("N0|N2|N4|N6", colnames(Rhizobia))]

# Prepare soil properties data
trait <- trait[, c(1, 5:25)]
rownames(trait) <- trait$sample
trait$sample <- NULL

# Transpose microbial data
TAMF <- t(AMF)
TRhizobia <- t(Rhizobia)

# Combine AMF and Rhizobia data
TASV <- cbind(TAMF, TRhizobia)

# Split data by niche with corrected patterns
ORI <- TASV[grepl("\\.O\\.", rownames(TASV)), ]      # Original soil
BS <- TASV[grepl("\\.Z\\.", rownames(TASV)), ]       # Bulk soil (Rootzone)
RH <- TASV[grepl("\\.RH\\.", rownames(TASV)), ]      # Rhizosphere
RT <- TASV[grepl("\\.RT\\.", rownames(TASV)), ]      # Root
NOD <- TASV[grepl("\\.N\\.", rownames(TASV)), ]      # Nodule

cat("\nSample counts by niche:\n")
cat("ORI (Original soil):", nrow(ORI), "samples\n")
cat("BS (Bulk soil):", nrow(BS), "samples\n")
cat("RH (Rhizosphere):", nrow(RH), "samples\n")
cat("RT (Root):", nrow(RT), "samples\n")
cat("NOD (Nodule):", nrow(NOD), "samples\n")

# =============================================================================
# 2. FIXED MANTEL ANALYSIS FUNCTION
# =============================================================================

perform_mantel_analysis <- function(niche_name, micro_data, trait_data) {
  cat(paste("\n=== Mantel analysis for", niche_name, "===\n"))
  
  # Define pattern mapping for each niche
  pattern_map <- list(
    "ORI" = "\\.O\\.",
    "BS" = "\\.Z\\.",
    "RH" = "\\.RH\\.",
    "RT" = "\\.RT\\.",
    "NOD" = "\\.N\\."
  )
  
  # Get the correct pattern for this niche
  pattern <- pattern_map[[niche_name]]
  
  if (is.null(pattern)) {
    cat("Error: No pattern defined for niche", niche_name, "\n")
    return(NULL)
  }
  
  # Filter trait data for the niche using the correct pattern
  trait_niche <- trait_data[grepl(pattern, rownames(trait_data)), ]
  
  cat("Trait samples found for", niche_name, ":", nrow(trait_niche), "\n")
  cat("Microbial samples for", niche_name, ":", nrow(micro_data), "\n")
  
  # Ensure samples match
  common_samples <- intersect(rownames(micro_data), rownames(trait_niche))
  
  cat("Common samples found:", length(common_samples), "\n")
  
  if (length(common_samples) == 0) {
    cat("No common samples found. Checking sample names...\n")
    cat("Microbial samples (first 5):", head(rownames(micro_data), 5), "\n")
    cat("Trait samples (first 5):", head(rownames(trait_niche), 5), "\n")
    return(NULL)
  }
  
  # Subset data to common samples
  micro_data <- micro_data[common_samples, ]
  trait_niche <- trait_niche[common_samples, ]
  
  if (nrow(trait_niche) < 2) {
    cat("Not enough trait samples for correlation analysis.\n")
    return(NULL)
  }
  
  # Calculate Mantel correlation
  cat("Calculating Mantel correlations...\n")
  mantel_result <- mantel_test(
    spec = micro_data, 
    env = trait_niche, 
    spec_select = list(
      AMF = 1:ncol(TAMF),
      Rhizobia = (ncol(TAMF) + 1):(ncol(TAMF) + ncol(TRhizobia))
    ), 
    mantel_fun = 'mantel'
  )
  
  # Add significance categories
  mantel_result <- mantel_result %>%
    mutate(
      rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
               labels = c('< 0.2', '0.2 - 0.4', '>= 0.4')),
      pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), 
               labels = c('< 0.01', '0.01 - 0.05', '>= 0.05'))
    )
  
  # Save results
  results_file <- file.path(results_dir, paste0("Mantel_results_", niche_name, ".csv"))
  write.csv(as.data.frame(mantel_result), results_file, row.names = FALSE)
  
  cat(paste("Mantel results saved as:", results_file, "\n"))
  
  # Print summary
  cat("\nSignificant correlations (p < 0.05):\n")
  sig_results <- mantel_result[mantel_result$p < 0.05, ]
  if (nrow(sig_results) > 0) {
    for (i in 1:nrow(sig_results)) {
      cat(paste0(
        sig_results$spec[i], " - ", sig_results$env[i], 
        ": r = ", round(sig_results$r[i], 3), 
        ", p = ", ifelse(sig_results$p[i] < 0.001, "<0.001", round(sig_results$p[i], 3)), "\n"
      ))
    }
  } else {
    cat("No significant correlations found.\n")
  }
  
  return(mantel_result)
}

# =============================================================================
# 3. PERFORM MANTEL ANALYSIS FOR EACH NICHE
# =============================================================================

cat("\n" , strrep("=", 60) , "\n")
cat("PERFORMING MANTEL ANALYSES\n")
cat(strrep("=", 60) , "\n")

niches <- list(
  "ORI" = ORI,
  "BS" = BS,
  "RH" = RH,
  "RT" = RT,
  "NOD" = NOD
)

all_results <- list()

for (niche_name in names(niches)) {
  if (nrow(niches[[niche_name]]) > 0) {
    results <- perform_mantel_analysis(niche_name, niches[[niche_name]], trait)
    if (!is.null(results)) {
      all_results[[niche_name]] <- results
    }
  } else {
    cat(paste("\nNo microbial data for", niche_name, "niche\n"))
  }
}

# =============================================================================
# 4. COMBINE AND PROCESS ALL RESULTS
# =============================================================================

cat("\n" , strrep("=", 60) , "\n")
cat("COMBINING ALL MANTEL RESULTS\n")
cat(strrep("=", 60) , "\n")

# Check which files were created
created_files <- list.files(results_dir, pattern = "Mantel_results_.*\\.csv", full.names = TRUE)
cat("\nCreated Mantel result files:\n")
print(created_files)

# Read Mantel results from each niche
if (length(created_files) > 0) {
  # Create a list to store the read data
  mantel_data <- list()
  
  for (file in created_files) {
    niche_name <- gsub("Mantel_results_|\\.csv", "", basename(file))
    cat("\nReading", niche_name, "...\n")
    data <- read.csv(file, header = TRUE)
    mantel_data[[niche_name]] <- data
  }
  
  # Combine all results
  cat("\nCombining results from all niches...\n")
  
  # Add niche label and id for each niche (using correct abbreviations)
  niche_abbrev <- list(
    "ORI" = "MFS",    # Original soil -> MFS
    "BS" = "RZS",     # Bulk soil -> RZS
    "RH" = "RS",      # Rhizosphere -> RS
    "RT" = "RE",      # Root -> RE
    "NOD" = "NE"      # Nodule -> NE
  )
  
  total_list <- list()
  
  for (niche_name in names(mantel_data)) {
    data <- mantel_data[[niche_name]]
    abbrev <- niche_abbrev[[niche_name]]
    
    if (!is.null(abbrev)) {
      data$niche <- abbrev
      data$id <- paste0(data$spec, "_", data$niche)
      total_list[[niche_name]] <- data
    }
  }
  
  if (length(total_list) > 0) {
    total <- do.call(rbind, total_list)
    
    # Separate r values and p values
    total_r <- total[c("id", "niche", "env", "r")]
    total_p <- total[c("id", "niche", "env", "p")]
    
    # Convert p values to significance markers
    b2 <- total_p
    sig <- ifelse(b2$p < 0.001, '***', 
                  ifelse(b2$p >= 0.001 & b2$p < 0.01, '**', 
                         ifelse(b2$p >= 0.01 & b2$p < 0.05, '*', '')))
    b2$p <- sig
    total_p <- b2
    
    # Convert data to wide format
    r_width <- total_r %>% 
      pivot_wider(names_from = env, values_from = r)
    
    p_width <- total_p %>% 
      pivot_wider(names_from = env, values_from = p)
    
    r_width <- as.data.frame(r_width)
    p_width <- as.data.frame(p_width)
    
    # Set row names
    rownames(r_width) <- r_width$id
    rownames(p_width) <- p_width$id
    
    # Arrange rows in specified order
    AMF_NAME <- r_width$id[grepl("AMF", r_width$id)]
    Rhi_NAME <- r_width$id[grepl("Rhizobia", r_width$id)]
    
    if (length(AMF_NAME) > 0 && length(Rhi_NAME) > 0) {
      order <- c(Rhi_NAME, AMF_NAME)
      r_width <- r_width[order, ]
      p_width <- p_width[order, ]
    }
    
    # Extract matrices for plotting
    if (ncol(r_width) > 3) {
      r_plot <- r_width[4:ncol(r_width)]  # Skip first 3 columns (id, niche, env)
      p_plot <- p_width[4:ncol(p_width)]
      
      # Transpose
      r_plot <- t(r_plot)
      p_plot <- t(p_plot)
      
      # Save results
      write.csv(r_plot, file.path(results_dir, "AMF_Rhizobia_mantel_R.csv"))
      write.csv(p_plot, file.path(results_dir, "AMF_Rhizobia_mantel_p.csv"))
      
      cat("\nCombined results saved as:\n")
      cat("1. AMF_Rhizobia_mantel_R.csv\n")
      cat("2. AMF_Rhizobia_mantel_p.csv\n")
      
      # Create summary statistics
      if (nrow(total) > 0) {
        summary_stats <- total %>%
          filter(p < 0.05) %>%
          group_by(niche, spec) %>%
          summarise(
            n_significant = n(),
            mean_r = round(mean(r), 3),
            min_p = min(p),
            max_r = max(r),
            .groups = 'drop'
          )
        
        # Save significant correlations summary
        write.csv(summary_stats, 
                  file.path(results_dir, "Mantel_significant_summary.csv"), 
                  row.names = FALSE)
        
        # Print summary
        cat("\nSummary of significant correlations:\n")
        print(summary_stats)
      }
    } else {
      cat("\nWarning: Not enough columns for pivot_wider\n")
    }
  } else {
    cat("\nNo data to combine.\n")
  }
} else {
  cat("\nNo Mantel result files found.\n")
}

# =============================================================================
# 5. SESSION INFO AND COMPLETION MESSAGE
# =============================================================================

cat("\n" , strrep("=", 60) , "\n")
cat("ANALYSIS COMPLETED\n")
cat(strrep("=", 60) , "\n")
cat(paste("Results saved in:", results_dir, "\n"))

cat("\nFiles generated:\n")
cat("1. Mantel_results_[Niche].csv for each niche\n")
cat("2. AMF_Rhizobia_mantel_R.csv\n")
cat("3. AMF_Rhizobia_mantel_p.csv\n")
cat("4. Mantel_significant_summary.csv\n")

# Session info
cat("\nSession information:\n")
sessionInfo()