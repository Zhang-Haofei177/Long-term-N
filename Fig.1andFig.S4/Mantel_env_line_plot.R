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

if (length(created_files) > 0) {
  mantel_data <- list()
  
  for (file in created_files) {
    niche_name <- gsub("Mantel_results_|\\.csv", "", basename(file))
    cat("\nReading", niche_name, "...\n")
    data <- read.csv(file, header = TRUE)
    mantel_data[[niche_name]] <- data
  }
  
  cat("\nCombining results from all niches...\n")
  
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
    
    total_r <- total[c("id", "niche", "env", "r")]
    total_p <- total[c("id", "niche", "env", "p")]
    
    b2 <- total_p
    sig <- ifelse(b2$p < 0.001, '***', 
                  ifelse(b2$p >= 0.001 & b2$p < 0.01, '**', 
                         ifelse(b2$p >= 0.01 & b2$p < 0.05, '*', '')))
    b2$p <- sig
    total_p <- b2
    
    r_width <- total_r %>% 
      pivot_wider(names_from = env, values_from = r)
    
    p_width <- total_p %>% 
      pivot_wider(names_from = env, values_from = p)
    
    r_width <- as.data.frame(r_width)
    p_width <- as.data.frame(p_width)
    
    rownames(r_width) <- r_width$id
    rownames(p_width) <- p_width$id
    
    AMF_NAME <- r_width$id[grepl("AMF", r_width$id)]
    Rhi_NAME <- r_width$id[grepl("Rhizobia", r_width$id)]
    
    if (length(AMF_NAME) > 0 && length(Rhi_NAME) > 0) {
      order <- c(Rhi_NAME, AMF_NAME)
      r_width <- r_width[order, ]
      p_width <- p_width[order, ]
    }
    
    if (ncol(r_width) > 3) {
      r_plot <- r_width[4:ncol(r_width)]  
      p_plot <- p_width[4:ncol(p_width)]
      
      r_plot <- t(r_plot)
      p_plot <- t(p_plot)
      
      write.csv(r_plot, file.path(results_dir, "AMF_Rhizobia_mantel_R.csv"))
      write.csv(p_plot, file.path(results_dir, "AMF_Rhizobia_mantel_p.csv"))
      
      cat("\nCombined results saved as:\n")
      cat("1. AMF_Rhizobia_mantel_R.csv\n")
      cat("2. AMF_Rhizobia_mantel_p.csv\n")
      
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
        
        write.csv(summary_stats, 
                  file.path(results_dir, "Mantel_significant_summary.csv"), 
                  row.names = FALSE)
        
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
# =============================================================================
# File: 06_Soil_properties_correlation_simple.R
# Description: Calculate Spearman correlations for Original niche soil properties
# =============================================================================

# Clean workspace
rm(list = ls())

# Load required libraries
library(tidyverse)
library(Hmisc)  # For rcorr function with Spearman correlation

# Set paths
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

# =============================================================================
# 1. LOAD AND PREPARE ORIGINAL NICHE DATA
# =============================================================================

cat("Loading and preparing Original niche soil properties data...\n")

# Load soil properties
trait <- read.csv(file.path(data_dir, "Mantel_env.csv"), header = TRUE)

# Prepare soil properties data
trait <- trait[, c(1, 5:25)]
rownames(trait) <- trait$sample

# Filter for Original niche (samples containing ".O.")
original_trait <- trait[grep("\\.O\\.", trait$sample), ]
rownames(original_trait) <- original_trait$sample
original_trait$sample <- NULL

# Check data structure
cat("\nOriginal niche data structure:\n")
cat("Number of samples:", nrow(original_trait), "\n")
cat("Number of variables:", ncol(original_trait), "\n")
cat("\nVariable names:\n")
print(colnames(original_trait))

# =============================================================================
# 2. CALCULATE SPEARMAN CORRELATIONS
# =============================================================================

cat("\nCalculating Spearman correlations...\n")

# Convert to numeric matrix
trait_matrix <- as.matrix(original_trait)

# Calculate Spearman correlation with p-values
cor_result <- rcorr(trait_matrix, type = "spearman")

# Extract correlation coefficients and p-values
cor_r <- cor_result$r
cor_p <- cor_result$P

# Replace NA values (diagonal with 1 for r, 0 for p)
diag(cor_r) <- 1
diag(cor_p) <- 0

# =============================================================================
# 3. RESHAPE AND FORMAT RESULTS
# =============================================================================

cat("\nFormatting results...\n")

# Convert to long format
cor_r_long <- melt(cor_r, varnames = c("tax", "Index"), value.name = "r")
cor_p_long <- melt(cor_p, varnames = c("tax", "Index"), value.name = "P_value")

# Merge r and p values
cor_combined <- merge(cor_r_long, cor_p_long, by = c("tax", "Index"))

# Add significance stars and format p-values
cor_combined <- cor_combined %>%
  mutate(
    P_value_sig = case_when(
      P_value < 0.001 ~ "***",
      P_value >= 0.001 & P_value < 0.01 ~ "**",
      P_value >= 0.01 & P_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Format P-value for display (scientific notation for small values)
    P_value_formatted = ifelse(
      P_value < 0.001,
      sprintf("%.2E", P_value),
      sprintf("%.8f", P_value)
    ),
    # Round r values to 8 decimal places
    r = round(r, 8)
  ) %>%
  # Order by tax and Index
  arrange(tax, Index) %>%
  # Select and rename columns to match requested format
  select(
    tax,
    Index,
    r,
    P_value = P_value_formatted,
    P_value_sig
  )

# =============================================================================
# 4. SAVE RESULTS
# =============================================================================

# Save the correlation table
output_file <- file.path(results_dir, "soil_Spearman_correlation.csv")
write.csv(cor_combined, output_file, row.names = FALSE, quote = FALSE)

cat("\n" , strrep("=", 60) , "\n")
cat("ANALYSIS COMPLETED\n")
cat(strrep("=", 60) , "\n")
cat(paste("\nResults saved as:", output_file, "\n"))

# Display first few rows
cat("\nFirst 10 rows of output:\n")
print(head(cor_combined, 10))

# =============================================================================
# File: 07_Create_mantel_link_table.R
# Description: Create mantel link table from correlation results
# =============================================================================

# Clean workspace
rm(list = ls())

# Load required libraries
library(tidyverse)
library(reshape2)

# Set paths
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

# =============================================================================
# 1. LOAD REQUIRED DATA
# =============================================================================

cat("Loading required data...\n")

# Load environment correlation data (from previous analysis)
env_cor <- read.csv(file.path(results_dir, "soil_Spearman_correlation.csv"), header = TRUE)

# Load mantel results
bio_r <- as.data.frame(t(read.csv(file.path(results_dir, "AMF_Rhizobia_mantel_R.csv"), row.names = 1)))
bio_p <- as.data.frame(t(read.csv(file.path(results_dir, "AMF_Rhizobia_mantel_p.csv"), row.names = 1)))

# =============================================================================
# 2. PREPARE ENVIRONMENT CORRELATION DATA
# =============================================================================

cat("\nPreparing environment correlation data...\n")

# Format environment correlation data to match expected structure
# We only need the diagonal and lower triangle for environment correlations
# But for mantel analysis, we need the correlation between environmental factors

# Create environment correlation matrix format
env_r_long <- env_cor %>%
  select(tax, Index, r) %>%
  mutate(
    .rownames = tax,
    .colnames = Index
  ) %>%
  select(.rownames, .colnames, r)

# Create environment p-value matrix format
env_p_long <- env_cor %>%
  select(tax, Index, P_value_sig) %>%
  mutate(
    .rownames = tax,
    .colnames = Index,
    p = ifelse(P_value_sig == "", "", P_value_sig)
  ) %>%
  select(.rownames, .colnames, p)

# Combine environment correlation data
env_data <- merge(env_r_long, env_p_long, by = c(".rownames", ".colnames"))

# =============================================================================
# 3. PREPARE BIO-ENV CORRELATION DATA (MANTEL RESULTS)
# =============================================================================

cat("\nPreparing bio-environment correlation data...\n")

# Reshape bio_r data from wide to long format
bio_r_long <- bio_r %>%
  rownames_to_column(var = ".rownames") %>%
  melt(id.vars = ".rownames", 
       variable.name = ".colnames", 
       value.name = "r")

# Reshape bio_p data from wide to long format
bio_p_long <- bio_p %>%
  rownames_to_column(var = ".rownames") %>%
  melt(id.vars = ".rownames", 
       variable.name = ".colnames", 
       value.name = "p")

# Combine bio data
bio_data <- merge(bio_r_long, bio_p_long, by = c(".rownames", ".colnames"))

# Ensure r values are numeric
bio_data$r <- as.numeric(bio_data$r)
# =============================================================================
# 4. CREATE MANTEL LINK TABLE
# =============================================================================
AMF_line_site <- read.csv(file.path(results_dir, "AMF_line_site.csv"))
Rhizobia_line_site <- read.csv(file.path(results_dir, "Rhizobia_line_site.csv"))
colnames(bio_p_long) <- c("bio","env","p")
colnames(bio_r_long) <- c("bio","env","r")
AMF_line_site_data <- merge(AMF_line_site, bio_r_long, by = c("bio", "env"), all.x = TRUE)
AMF_line_site_data <- merge(AMF_line_site_data, bio_p_long, by = c("bio", "env"), all.x = TRUE)
AMF_mantel_line <- AMF_line_site_data %>%
  mutate(
    r.sign = ifelse(r > 0, "Positive", "Negative"),
    r.abs = cut(abs(r), 
                breaks = c(0, 0.25, 0.5, 1), 
                labels = c("<0.25", "0.25-0.5", "0.5-1"),
                include.lowest = TRUE,
                right = FALSE),
    p.sign = p
  )

head(AMF_mantel_line)
AMF_mantel_line <- AMF_mantel_line %>%
  select(bio, env, ymin, xmin, xmax, ymax, r, p, r.sign, r.abs, p.sign)
output_file <- file.path(results_dir, "AMF_mantel_line.csv")
write.csv(AMF_mantel_line, output_file, row.names = FALSE)


Rhizobia_line_site_data <- merge(Rhizobia_line_site, bio_r_long, by = c("bio", "env"), all.x = TRUE)
Rhizobia_line_site_data <- merge(Rhizobia_line_site_data, bio_p_long, by = c("bio", "env"), all.x = TRUE)
Rhizobia_mantel_line <- Rhizobia_line_site_data %>%
  mutate(
    r.sign = ifelse(r > 0, "Positive", "Negative"),
    r.abs = cut(abs(r), 
                breaks = c(0, 0.25, 0.5, 1), 
                labels = c("<0.25", "0.25-0.5", "0.5-1"),
                include.lowest = TRUE,
                right = FALSE),
    p.sign = p
  )

head(Rhizobia_mantel_line)
Rhizobia_mantel_line <- Rhizobia_mantel_line %>%
  select(bio, env, ymin, xmin, xmax, ymax, r, p, r.sign, r.abs, p.sign)
output_file <- file.path(results_dir, "Rhizobia_mantel_line.csv")
write.csv(Rhizobia_mantel_line, output_file, row.names = FALSE)
####################draw plot for AMF and Rhizobia########################
# Clean workspace
rm(list = ls())
# Set paths
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

cat("Loading and preparing data...\n")

# Load environment correlation data
env <- read.csv(file.path(data_dir, "soil_Spearman_correlation.csv"), header = TRUE)
env <- env[c(1,2,3,4)]
colnames(env) <- c(".rownames",".colnames","r","p")

# Load mantel results
bio_r <- read.csv(file.path(data_dir, "AMF_Rhizobia_mantel_R.csv"), header = TRUE)
bio_p <- read.csv(file.path(data_dir, "AMF_Rhizobia_mantel_p.csv"), header = TRUE)

# Reshape mantel data
bio_r <- melt(bio_r, id.vars = "X", variable.name = ".colnames", value.name = "r")
colnames(bio_r)[1] <- ".rownames"
bio_p <- melt(bio_p, id.vars = "X", variable.name = ".colnames", value.name = "p")
colnames(bio_p)[1] <- ".rownames"
bio <- cbind(bio_r, bio_p[3])

# Load mantel line data
AMF_line <- read.csv(file.path(data_dir, "AMF_mantel_line.csv"), header = TRUE)
Rhi_line <- read.csv(file.path(data_dir, "Rhizobia_mantel_line.csv"), header = TRUE)

# =============================================================================
# 2. CREATE ENVIRONMENT CORRELATION HEATMAP
# =============================================================================

cat("Creating environment correlation heatmap...\n")

# Prepare correlation matrices
cor2_1 <- env
rownames_unique <- unique(cor2_1$.rownames)
colnames_unique <- unique(cor2_1$.colnames)

r_matrix <- matrix(NA, nrow = length(rownames_unique), ncol = length(colnames_unique),                   
                   dimnames = list(rownames_unique, colnames_unique))
p_matrix <- matrix(NA, nrow = length(rownames_unique), ncol = length(colnames_unique),                   
                   dimnames = list(rownames_unique, colnames_unique))

for (i in 1:nrow(cor2_1)) {  
  r_matrix[cor2_1$.rownames[i], cor2_1$.colnames[i]] <- cor2_1$r[i]  
  p_matrix[cor2_1$.rownames[i], cor2_1$.colnames[i]] <- cor2_1$p[i]
}

# Keep only upper triangle
cor_matrix <- r_matrix
cor_matrix[lower.tri(cor_matrix)] <- NA
p_cor_matrix <- p_matrix
p_cor_matrix[lower.tri(p_cor_matrix)] <- NA

# Convert to long format
long_cor <- melt(cor_matrix, na.rm = TRUE)
long_p <- melt(p_cor_matrix, na.rm = TRUE)

# Add significance stars
long_p <- long_p %>%
  mutate(stars = case_when(
    value < 0.001 ~ "***",
    value < 0.01  ~ "**",
    value < 0.05  ~ "*",
    value >= 0.05 ~ "ns"
  ))

# Create heatmap
p1 <- ggplot() +
  geom_tile(data = long_cor, aes(Var2, Var1), fill = "white", color = "grey") +
  geom_point(data = long_cor, aes(Var2, Var1, size = abs(value), color = value), shape = 15) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -135, hjust = 1),
        axis.text.y = element_text(angle = -135, hjust = 1),
        legend.text = element_text(angle = 0),
        legend.title = element_text(angle = 0)) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdYlBu")), limits = c(-1, 1)) +
  geom_text(data = long_p, aes(x = Var2, y = Var1, label = stars), 
            size = 3, color = "black", family = "sans", angle = -135) +
  coord_fixed(ratio = 1)

# =============================================================================
# 3. CREATE AMF MANTEL LINK PLOT
# =============================================================================

cat("Creating AMF mantel link plot...\n")

# Filter AMF data (excluding MFS if needed)
AMF_line_plot <- AMF_line[grepl("AMF", AMF_line$bio), ]
AMF_line_plot <- AMF_line_plot[!grepl("MFS", AMF_line_plot$bio), ]

# Ensure p.sign has proper values for plotting
AMF_line_plot$p.sign <- ifelse(AMF_line_plot$p.sign == "", "NA", AMF_line_plot$p.sign)

plot_line_AMF <- ggplot() + 
  geom_curve(data = AMF_line_plot, 
             aes(x = xmin, y = ymin, xend = xmax, yend = ymax, 
                 color = p.sign, size = r.abs), 
             curvature = 0.1) + 
  geom_point(data = AMF_line_plot, aes(x = xmin, y = ymin)) + 
  geom_point(data = AMF_line_plot, aes(x = xmax, y = ymax)) + 
  theme_void() + 
  geom_text(data = AMF_line_plot, 
            aes(x = xmin, y = ymin, label = env), 
            hjust = 1.2, vjust = -0.5, size = 3, 
            color = "black", angle = -135) + 
  geom_text(data = AMF_line_plot, 
            aes(x = xmax, y = ymax, label = bio), 
            hjust = -0.2, vjust = -0.5, size = 3, 
            color = "black", angle = -135) + 
  scale_size_manual(values = c("<0.25" = 0.5,                               
                               "0.25-0.5" = 1.5,                               
                               "0.5-1" = 3),
                    name = "|r|") +  
  scale_colour_manual(values = c("***" = "#283c63",                                 
                                 "**" = "#D73027",                     
                                 "*" = "#4575B7",                                 
                                 "NA" = "gray"),
                      name = "Significance") +
  labs(title = "AMF")

# =============================================================================
# 4. CREATE RHIZOBIA MANTEL LINK PLOT
# =============================================================================

cat("Creating Rhizobia mantel link plot...\n")

# Filter Rhizobia data
Rhi_line_plot <- Rhi_line[grepl("Rhizobia", Rhi_line$bio), ]
Rhi_line_plot <- Rhi_line_plot[!grepl("MFS", Rhi_line_plot$bio), ]

# Ensure p.sign has proper values for plotting
Rhi_line_plot$p.sign <- ifelse(Rhi_line_plot$p.sign == "", "NA", Rhi_line_plot$p.sign)

plot_line_Rhi <- ggplot() + 
  geom_curve(data = Rhi_line_plot, 
             aes(x = xmin, y = ymin, xend = xmax, yend = ymax, 
                 color = p.sign, size = r.abs), 
             curvature = -0.1) + 
  geom_point(data = Rhi_line_plot, aes(x = xmin, y = ymin)) + 
  geom_point(data = Rhi_line_plot, aes(x = xmax, y = ymax)) + 
  theme_void() + 
  geom_text(data = Rhi_line_plot, 
            aes(x = xmin, y = ymin, label = env), 
            hjust = 1.2, vjust = -0.5, size = 3, 
            color = "black", angle = -45) + 
  geom_text(data = Rhi_line_plot, 
            aes(x = xmax, y = ymax, label = bio), 
            hjust = -0.2, vjust = -0.5, size = 3, 
            color = "black", angle = -45) + 
  scale_size_manual(values = c("<0.25" = 0.5,                               
                               "0.25-0.5" = 1.5,                               
                               "0.5-1" = 3),
                    name = "|r|") +  
  scale_colour_manual(values = c("***" = "#283c63",                                 
                                 "**" = "#D73027",                     
                                 "*" = "#4575B7",                                 
                                 "NA" = "gray"),
                      name = "Significance") +
  labs(title = "Rhizobia")

# =============================================================================
# 5. COMBINE AND SAVE PLOTS
# =============================================================================

cat("Combining and saving plots...\n")

# 5.1 AMF combined plot
p_total_AMF <- p1 | plot_line_AMF
output_file_AMF <- file.path(results_dir, "AMF_mantel.pdf")
ggsave(plot = p_total_AMF, output_file_AMF, device = "pdf", width = 10, height = 7, dpi = 300)
cat("AMF plot saved as:", output_file_AMF, "\n")

# 5.2 Rhizobia combined plot
p_total_Rhi <- p1 | plot_line_Rhi
output_file_Rhi <- file.path(results_dir, "Rhizobia_mantel.pdf")
ggsave(plot = p_total_Rhi, output_file_Rhi, device = "pdf", width = 10, height = 7, dpi = 300)
cat("Rhizobia plot saved as:", output_file_Rhi, "\n")

# 5.3 Combined all plots
p_combined_all <- (p1 | plot_line_AMF | plot_line_Rhi) +
  plot_layout(widths = c(1, 0.8, 0.8))
output_file_all <- file.path(results_dir, "All_mantel_combined.pdf")
ggsave(plot = p_combined_all, output_file_all, device = "pdf", width = 15, height = 7, dpi = 300)
cat("Combined plot saved as:", output_file_all, "\n")

# 5.4 Create separate plots for each microbe with different curvature
# AMF with positive curvature
plot_line_AMF_pos <- ggplot() + 
  geom_curve(data = AMF_line_plot, 
             aes(x = xmin, y = ymin, xend = xmax, yend = ymax, 
                 color = p.sign, size = r.abs), 
             curvature = 0.1) + 
  geom_point(data = AMF_line_plot, aes(x = xmin, y = ymin)) + 
  geom_point(data = AMF_line_plot, aes(x = xmax, y = ymax)) + 
  theme_void() + 
  geom_text(data = AMF_line_plot, 
            aes(x = xmin, y = ymin, label = env), 
            hjust = 1.2, vjust = -0.5, size = 3, 
            color = "black", angle = -135) + 
  geom_text(data = AMF_line_plot, 
            aes(x = xmax, y = ymax, label = bio), 
            hjust = -0.2, vjust = -0.5, size = 3, 
            color = "black", angle = -135) + 
  scale_size_manual(values = c("<0.25" = 0.5, "0.25-0.5" = 1.5, "0.5-1" = 3)) +  
  scale_colour_manual(values = c("***" = "#283c63", "**" = "#D73027", "*" = "#4575B7", "NA" = "gray")) +
  labs(title = "AMF (Positive Curvature)")

# Rhizobia with negative curvature
plot_line_Rhi_neg <- ggplot() + 
  geom_curve(data = Rhi_line_plot, 
             aes(x = xmin, y = ymin, xend = xmax, yend = ymax, 
                 color = p.sign, size = r.abs), 
             curvature = -0.1) + 
  geom_point(data = Rhi_line_plot, aes(x = xmin, y = ymin)) + 
  geom_point(data = Rhi_line_plot, aes(x = xmax, y = ymax)) + 
  theme_void() + 
  geom_text(data = Rhi_line_plot, 
            aes(x = xmin, y = ymin, label = env), 
            hjust = 1.2, vjust = -0.5, size = 3, 
            color = "black", angle = -45) + 
  geom_text(data = Rhi_line_plot, 
            aes(x = xmax, y = ymax, label = bio), 
            hjust = -0.2, vjust = -0.5, size = 3, 
            color = "black", angle = -45) + 
  scale_size_manual(values = c("<0.25" = 0.5, "0.25-0.5" = 1.5, "0.5-1" = 3)) +  
  scale_colour_manual(values = c("***" = "#283c63", "**" = "#D73027", "*" = "#4575B7", "NA" = "gray")) +
  labs(title = "Rhizobia (Negative Curvature)")

# Save separate plots
ggsave(file.path(results_dir, "AMF_mantel_separate.pdf"), 
       plot_line_AMF_pos, device = "pdf", width = 8, height = 6, dpi = 300)
ggsave(file.path(results_dir, "Rhizobia_mantel_separate.pdf"), 
       plot_line_Rhi_neg, device = "pdf", width = 8, height = 6, dpi = 300)

# =============================================================================
# 6. PRINT SUMMARY
# =============================================================================

cat("\n" , strrep("=", 60) , "\n")
cat("PLOTS GENERATED SUCCESSFULLY\n")
cat(strrep("=", 60) , "\n")

cat("\nFiles created:\n")
cat("1. AMF_mantel.pdf - Combined heatmap + AMF mantel plot\n")
cat("2. Rhizobia_mantel.pdf - Combined heatmap + Rhizobia mantel plot\n")
cat("3. All_mantel_combined.pdf - All three plots combined\n")
cat("4. AMF_mantel_separate.pdf - AMF mantel plot only\n")
cat("5. Rhizobia_mantel_separate.pdf - Rhizobia mantel plot only\n")

cat("\nSummary statistics:\n")
cat("AMF correlations:", nrow(AMF_line_plot), "\n")
cat("Rhizobia correlations:", nrow(Rhi_line_plot), "\n")
cat("Environment factors:", length(unique(c(AMF_line_plot$env, Rhi_line_plot$env))), "\n")
cat("Biological factors (AMF):", length(unique(AMF_line_plot$bio)), "\n")
cat("Biological factors (Rhizobia):", length(unique(Rhi_line_plot$bio)), "\n")
