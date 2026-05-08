# =============================================================================
# File: 05_Mantel_analysis_fixed_FDR.R
# Description: Mantel correlation analysis with FDR correction for multiple testing
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
# 2. MANTEL ANALYSIS FUNCTION (no correction within each niche)
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
  
  pattern <- pattern_map[[niche_name]]
  if (is.null(pattern)) {
    cat("Error: No pattern defined for niche", niche_name, "\n")
    return(NULL)
  }
  
  # Filter trait data for the niche
  trait_niche <- trait_data[grepl(pattern, rownames(trait_data)), ]
  cat("Trait samples found:", nrow(trait_niche), "\n")
  cat("Microbial samples:", nrow(micro_data), "\n")
  
  common_samples <- intersect(rownames(micro_data), rownames(trait_niche))
  cat("Common samples:", length(common_samples), "\n")
  
  if (length(common_samples) == 0) return(NULL)
  
  micro_data <- micro_data[common_samples, ]
  trait_niche <- trait_niche[common_samples, ]
  
  if (nrow(trait_niche) < 2) return(NULL)
  
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
  
  # No correction here; correction will be applied globally later
  return(mantel_result)
}

# =============================================================================
# 3. RUN MANTEL ANALYSIS FOR EACH NICHE AND COLLECT RESULTS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PERFORMING MANTEL ANALYSES (uncorrected p-values)\n")
cat(strrep("=", 60), "\n")

niches <- list(
  "ORI" = ORI,
  "BS" = BS,
  "RH" = RH,
  "RT" = RT,
  "NOD" = NOD
)

all_results_list <- list()

for (niche_name in names(niches)) {
  if (nrow(niches[[niche_name]]) > 0) {
    res <- perform_mantel_analysis(niche_name, niches[[niche_name]], trait)
    if (!is.null(res)) {
      res$niche <- niche_name
      all_results_list[[niche_name]] <- res
    }
  } else {
    cat(paste("\nNo microbial data for", niche_name, "niche\n"))
  }
}

# =============================================================================
# 4. COMBINE ALL RESULTS AND APPLY GLOBAL FDR CORRECTION
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("COMBINING ALL RESULTS AND APPLYING FDR CORRECTION\n")
cat(strrep("=", 60), "\n")

if (length(all_results_list) > 0) {
  # Combine all data frames
  total <- bind_rows(all_results_list, .id = "niche_raw")
  
  # Add standard niche abbreviations
  niche_abbrev <- c(
    "ORI" = "MFS",
    "BS" = "RZS",
    "RH" = "RS",
    "RT" = "RE",
    "NOD" = "NE"
  )
  total$niche <- niche_abbrev[total$niche_raw]
  total$id <- paste0(total$spec, "_", total$niche)
  
  # Apply FDR correction to all p-values (Benjamini-Hochberg)
  total$p_adj <- p.adjust(total$p, method = "fdr")
  
  # Add significance stars based on adjusted p-values
  total$p_sig <- case_when(
    total$p_adj < 0.001 ~ "***",
    total$p_adj < 0.01  ~ "**",
    total$p_adj < 0.05  ~ "*",
    TRUE                ~ ""
  )
  
  # Also keep the original r and env columns
  total_r <- total[, c("id", "niche", "env", "r")]
  total_p <- total[, c("id", "niche", "env", "p_sig")]
  
  # Convert to wide format for heatmap
  r_width <- total_r %>% pivot_wider(names_from = env, values_from = r)
  p_width <- total_p %>% pivot_wider(names_from = env, values_from = p_sig)
  
  r_width <- as.data.frame(r_width)
  p_width <- as.data.frame(p_width)
  rownames(r_width) <- r_width$id
  rownames(p_width) <- p_width$id
  
  # Order rows: Rhizobia first, then AMF
  order_rows <- c(
    rownames(r_width)[grepl("Rhizobia", rownames(r_width))],
    rownames(r_width)[grepl("AMF", rownames(r_width))]
  )
  r_width <- r_width[order_rows, ]
  p_width <- p_width[order_rows, ]
  
  # Remove id, niche, env columns? Actually we already have columns.
  # For matrix output, select only environmental factor columns (after column 3)
  env_cols <- setdiff(names(r_width), c("id", "niche"))
  r_matrix <- r_width[, env_cols]
  p_matrix <- p_width[, env_cols]
  
  # Transpose to have env as rows, microbial groups as columns (common for plotting)
  r_matrix <- t(r_matrix)
  p_matrix <- t(p_matrix)
  
  # Save results with _FDR suffix
  write.csv(r_matrix, file.path(results_dir, "AMF_Rhizobia_mantel_R_FDR.csv"))
  write.csv(p_matrix, file.path(results_dir, "AMF_Rhizobia_mantel_p_FDR.csv"))
  
  cat("\nFDR-corrected results saved:\n")
  cat("1. AMF_Rhizobia_mantel_R_FDR.csv\n")
  cat("2. AMF_Rhizobia_mantel_p_FDR.csv\n")
  
  # Create summary of significant correlations (based on adjusted p < 0.05)
  summary_stats <- total %>%
    filter(p_adj < 0.05) %>%
    group_by(niche, spec) %>%
    summarise(
      n_significant = n(),
      mean_r = round(mean(r), 3),
      min_p_adj = min(p_adj),
      max_r = max(r),
      .groups = 'drop'
    )
  
  write.csv(summary_stats, 
            file.path(results_dir, "Mantel_significant_summary_FDR.csv"), 
            row.names = FALSE)
  cat("3. Mantel_significant_summary_FDR.csv\n")
  
  # Optional: also save the full corrected table
  write.csv(total, file.path(results_dir, "Mantel_all_results_FDR.csv"), row.names = FALSE)
  cat("4. Mantel_all_results_FDR.csv (complete table)\n")
  
} else {
  cat("No Mantel results to combine.\n")
}

# =============================================================================
# 5. COMPLETION
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("ANALYSIS WITH FDR CORRECTION COMPLETED\n")
cat(strrep("=", 60), "\n")
cat(paste("Results saved in:", results_dir, "\n"))
cat("\nAll final output files include '_FDR' suffix to indicate\nthe application of Benjamini-Hochberg false discovery rate correction.\n\n")
sessionInfo()

# =============================================================================
# File: Plot_Mantel_FDR.R
# Description: Plotting Mantel results with FDR-corrected significance
#              Uses *_FDR.csv files and outputs *_FDR.pdf
# =============================================================================

# Clean workspace
rm(list = ls())

# Load required libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(tidyr)

# Set paths (adjust as needed)
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

# =============================================================================
# 1. LOAD AND PREPARE DATA (using FDR-corrected files)
# =============================================================================

cat("Loading FDR-corrected Mantel results and environment correlations...\n")

# Load environment correlation data (unchanged)
env <- read.csv(file.path(results_dir, "soil_Spearman_correlation.csv"), header = TRUE)
env <- env[, c("tax", "Index", "r", "P_value_sig")]
colnames(env) <- c(".rownames", ".colnames", "r", "p")
env$p[env$p == ""] <- "ns"

# Load FDR-corrected Mantel matrices
bio_r <- read.csv(file.path(results_dir, "AMF_Rhizobia_mantel_R_FDR.csv"), header = TRUE, row.names = 1)
bio_p <- read.csv(file.path(results_dir, "AMF_Rhizobia_mantel_p_FDR.csv"), header = TRUE, row.names = 1)

# Reshape mantel data to long format, and rename columns for easy merging
bio_r_long <- as.data.frame(t(bio_r)) %>%
  rownames_to_column(var = "bio") %>%
  pivot_longer(cols = -bio, names_to = "env", values_to = "r")

bio_p_long <- as.data.frame(t(bio_p)) %>%
  rownames_to_column(var = "bio") %>%
  pivot_longer(cols = -bio, names_to = "env", values_to = "p")

# Merge r and p (optional, but we will keep separate for later)
bio <- left_join(bio_r_long, bio_p_long, by = c("bio", "env"))

# Load line site coordinates (unchanged)
AMF_line_site <- read.csv(file.path(results_dir, "AMF_line_site.csv"))
Rhizobia_line_site <- read.csv(file.path(results_dir, "Rhizobia_line_site.csv"))

# =============================================================================
# Environment correlation heatmap (direct Spearman calculation with FDR correction)
# =============================================================================

cat("Calculating Spearman correlations for environmental factors with FDR correction...\n")

# Load raw soil property data
library(Hmisc)
trait_raw <- read.csv(file.path(data_dir, "Mantel_env.csv"), header = TRUE)
# Select the same columns as before
trait <- trait_raw[, c(1, 5:25)]
rownames(trait) <- trait$sample
trait$sample <- NULL
# Subset to original soil niche (.O.)
trait_original <- trait[grep("\\.O\\.", rownames(trait)), ]
trait_original <- as.matrix(trait_original)

# Compute Spearman correlation matrix (including raw p-values)
cor_result <- rcorr(trait_original, type = "spearman")
r_orig <- cor_result$r
p_orig <- cor_result$P

# Apply FDR correction to the upper triangle (pairwise comparisons)
p_upper <- p_orig[upper.tri(p_orig, diag = FALSE)]
# Remove possible NAs
valid <- !is.na(p_upper)
if (sum(valid) > 0) {
  p_adj_upper <- p.adjust(p_upper[valid], method = "fdr")
  p_upper[valid] <- p_adj_upper
}
# Put corrected p-values back to upper triangle and symmetrize
p_adj <- p_orig
p_adj[upper.tri(p_adj, diag = FALSE)] <- p_upper
p_adj[lower.tri(p_adj, diag = FALSE)] <- t(p_adj)[lower.tri(p_adj, diag = FALSE)]
# Set diagonal to 0 (not used in plotting)
diag(p_adj) <- 0

# Build matrices for plotting (keep only upper triangle, lower triangle = NA)
r_plot <- r_orig
r_plot[lower.tri(r_plot)] <- NA
p_plot <- p_adj
p_plot[lower.tri(p_plot)] <- NA

# Convert to long format
long_cor <- melt(r_plot, na.rm = TRUE, varnames = c("Var1", "Var2"), value.name = "r")
long_p <- melt(p_plot, na.rm = TRUE, varnames = c("Var1", "Var2"), value.name = "p")

long_p <- long_p %>%
  mutate(stars = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# Heatmap
p1 <- ggplot() +
  geom_tile(data = long_cor, aes(Var2, Var1), fill = "white", color = "grey") +
  geom_point(data = long_cor, aes(Var2, Var1, size = abs(r), color = r), shape = 15) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = -135, hjust = 1),
        axis.text.y = element_text(angle = -135, hjust = 1)) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(colors = rev(brewer.pal(9, "RdYlBu")), limits = c(-1, 1)) +
  geom_text(data = long_p, aes(x = Var2, y = Var1, label = stars),
            size = 3, color = "black", angle = -135) +
  coord_fixed(ratio = 1)

p1
# =============================================================================
# 3. CREATE MANTEL LINK TABLES FOR AMF AND RHIZOBIA (FDR-based)
# =============================================================================

cat("Creating mantel link tables using FDR-corrected significance...\n")

# Helper function to merge coordinates with Mantel results
create_mantel_line <- function(line_site, bio_r_long, bio_p_long, microbe_name) {
  # Merge r values
  line_data <- merge(line_site, bio_r_long, by = c("bio", "env"), all.x = TRUE)
  # Merge p values
  line_data <- merge(line_data, bio_p_long, by = c("bio", "env"), all.x = TRUE)
  line_data <- line_data %>%
    mutate(
      r.sign = ifelse(r > 0, "Positive", "Negative"),
      r.abs = cut(abs(r),
                  breaks = c(0, 0.25, 0.5, 1),
                  labels = c("<0.25", "0.25-0.5", "0.5-1"),
                  include.lowest = TRUE,
                  right = FALSE),
      p.sign = ifelse(is.na(p) | p == "", "NA", p)
    ) %>%
    select(bio, env, ymin, xmin, xmax, ymax, r, p, r.sign, r.abs, p.sign)
  return(line_data)
}

# Note: bio_r_long and bio_p_long already contain bio, env, r and p
AMF_mantel_line <- create_mantel_line(AMF_line_site, bio_r_long, bio_p_long, "AMF")
Rhizobia_mantel_line <- create_mantel_line(Rhizobia_line_site, bio_r_long, bio_p_long, "Rhizobia")

# Save FDR-based line tables
write.csv(AMF_mantel_line, file.path(results_dir, "AMF_mantel_line_FDR.csv"), row.names = FALSE)
write.csv(Rhizobia_mantel_line, file.path(results_dir, "Rhizobia_mantel_line_FDR.csv"), row.names = FALSE)
cat("Saved: AMF_mantel_line_FDR.csv and Rhizobia_mantel_line_FDR.csv\n")

# =============================================================================
# 4. CREATE MANTEL LINK PLOTS USING FDR-BASED LINES
# =============================================================================

cat("Creating mantel link plots with FDR-corrected significance...\n")

# Filter data (remove MFS if present)
AMF_line_plot <- AMF_mantel_line[grepl("AMF", AMF_mantel_line$bio), ]
AMF_line_plot <- AMF_line_plot[!grepl("MFS", AMF_line_plot$bio), ]
AMF_line_plot$p.sign <- ifelse(AMF_line_plot$p.sign == "", "NA", AMF_line_plot$p.sign)

Rhi_line_plot <- Rhizobia_mantel_line[grepl("Rhizobia", Rhizobia_mantel_line$bio), ]
Rhi_line_plot <- Rhi_line_plot[!grepl("MFS", Rhi_line_plot$bio), ]
Rhi_line_plot$p.sign <- ifelse(Rhi_line_plot$p.sign == "", "NA", Rhi_line_plot$p.sign)

# AMF plot (curvature 0.1, angle -135)
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
  scale_size_manual(values = c("<0.25" = 0.5, "0.25-0.5" = 1.5, "0.5-1" = 3),
                    name = "|r|") +
  scale_colour_manual(values = c("***" = "#283c63", "**" = "#D73027",
                                 "*" = "#4575B7", "NA" = "gray"),
                      name = "Significance") +
  labs(title = "AMF")

# Rhizobia plot (curvature -0.1, angle -45)
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
  scale_size_manual(values = c("<0.25" = 0.5, "0.25-0.5" = 1.5, "0.5-1" = 3),
                    name = "|r|") +
  scale_colour_manual(values = c("***" = "#283c63", "**" = "#D73027",
                                 "*" = "#4575B7", "NA" = "gray"),
                      name = "Significance") +
  labs(title = "Rhizobia")

# =============================================================================
# 5. COMBINE AND SAVE PLOTS (with _FDR suffix)
# =============================================================================

cat("Combining and saving plots (FDR versions)...\n")

# AMF combined
p_total_AMF <- p1 | plot_line_AMF
ggsave(file.path(results_dir, "AMF_mantel_FDR.pdf"),
       plot = p_total_AMF, device = "pdf", width = 10, height = 7, dpi = 300)
cat("Saved: AMF_mantel_FDR.pdf\n")

# Rhizobia combined
p_total_Rhi <- p1 | plot_line_Rhi
ggsave(file.path(results_dir, "Rhizobia_mantel_FDR.pdf"),
       plot = p_total_Rhi, device = "pdf", width = 10, height = 7, dpi = 300)
cat("Saved: Rhizobia_mantel_FDR.pdf\n")

# All three combined
p_all <- (p1 | plot_line_AMF | plot_line_Rhi) +
  plot_layout(widths = c(1, 0.8, 0.8))
ggsave(file.path(results_dir, "All_mantel_combined_FDR.pdf"),
       plot = p_all, device = "pdf", width = 15, height = 7, dpi = 300)
cat("Saved: All_mantel_combined_FDR.pdf\n")

# Also save separate mantel link plots (optional)
ggsave(file.path(results_dir, "AMF_mantel_separate_FDR.pdf"),
       plot = plot_line_AMF, device = "pdf", width = 8, height = 6, dpi = 300)
ggsave(file.path(results_dir, "Rhizobia_mantel_separate_FDR.pdf"),
       plot = plot_line_Rhi, device = "pdf", width = 8, height = 6, dpi = 300)
cat("Saved: AMF_mantel_separate_FDR.pdf and Rhizobia_mantel_separate_FDR.pdf\n")

# =============================================================================
# 6. SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("FDR-BASED MANTEL PLOTS GENERATED SUCCESSFULLY\n")
cat(strrep("=", 60), "\n")
cat("All output files have '_FDR' suffix to indicate use of FDR-corrected p-values.\n")
cat("Files created:\n")
cat("  AMF_mantel_FDR.pdf\n")
cat("  Rhizobia_mantel_FDR.pdf\n")
cat("  All_mantel_combined_FDR.pdf\n")
cat("  AMF_mantel_separate_FDR.pdf\n")
cat("  Rhizobia_mantel_separate_FDR.pdf\n")
cat("  AMF_mantel_line_FDR.csv\n")
cat("  Rhizobia_mantel_line_FDR.csv\n")