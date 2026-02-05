# =============================================================================
# File: 04_PCoA_analysis.R
# Description: PCoA analysis for Bacteria, AMF, and Rhizobia
# =============================================================================

# Clean workspace
rm(list = ls())
library(tidyr)
library(reshape2)
library(ggplot2)
library(vegan)

# Set paths
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

# Load group data once
group <- read.csv(file.path(data_dir, "group.csv"), header = TRUE)

# Function to perform PCoA analysis for a microbe
perform_pcoa_analysis <- function(microbe_name) {
  cat(paste("\n=== Analyzing", microbe_name, "===\n"))
  
  # Load data
  a <- read.csv(file.path(data_dir, paste0(microbe_name, "_ASV_rarefaction.csv")), 
                header = TRUE, row.names = 1)
  
  # Select samples (excluding Original/O samples)
  a2 <- a[grepl("N0|N2|N4|N6", colnames(a))]
  a2 <- a2[!grepl("O", colnames(a2))]
  otu <- t(a2)
  
  # Prepare data frame
  df <- as.data.frame(otu)
  df$sample <- rownames(df)
  
  # Merge with group data
  a_g <- merge(df, group, by = "sample")
  rownames(a_g) <- a_g$sample
  a_g$sample <- NULL
  
  # Remove Original samples from group data
  a_g <- a_g[!grepl("O", rownames(a_g)), ]
  
  # Calculate Bray-Curtis distance
  distance <- vegdist(otu, method = 'bray')
  
  # Perform PCoA
  pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)
  pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
  
  # Prepare sample site data
  sample_site <- data.frame({pcoa$point})[1:2]
  sample_site$sample <- rownames(sample_site)
  names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
  sample_site <- merge(sample_site, group, by = 'sample', all.x = TRUE)
  
  # Set factor levels
  sample_site$treatment <- factor(sample_site$treatment, 
                                  levels = c("N0", "N2", "N4", "N6"))
  sample_site$Niche <- factor(sample_site$Niche,
                              levels = c("Rootzone", "Rhizosphere", "Root", "Nodule"),
                              labels = c("RZS", "RS", "RE", "NE"))
  
  # Perform PERMANOVA tests
  # Test 1: Niche effect
  adonis_niche <- adonis2(distance ~ a_g$Niche, distance = "bray", permutations = 999)
  r2_niche <- round(adonis_niche$R2[1], 3)
  p_niche <- adonis_niche$`Pr(>F)`[1]
  p_niche_formatted <- ifelse(p_niche < 0.001, "<0.001", sprintf("%.3f", p_niche))
  
  # Test 2: Treatment effect
  adonis_treatment <- adonis2(distance ~ a_g$treatment, distance = "bray", permutations = 999)
  r2_treatment <- round(adonis_treatment$R2[1], 3)
  p_treatment <- adonis_treatment$`Pr(>F)`[1]
  p_treatment_formatted <- ifelse(p_treatment < 0.001, "<0.001", sprintf("%.3f", p_treatment))
  
  # Calculate Niche Nitrogen Effects Index
  nne_index <- round(r2_niche / r2_treatment, 3)
  
  # Create annotation text
  annotation_text <- paste0(
    "Niche: R² = ", r2_niche, ", p = ", p_niche_formatted, "\n",
    "Treatment: R² = ", r2_treatment, ", p = ", p_treatment_formatted, "\n",
    "NNE Index = ", nne_index
  )
  
  # Create plot
  p <- ggplot(sample_site, aes(PCoA1, PCoA2)) +
    theme(
      panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
      panel.background = element_rect(color = 'black', fill = 'transparent'), 
      legend.key = element_rect(fill = 'transparent'),
      text = element_text(family = "sans"),  # Changed from "serif" to "sans"
      legend.spacing.y = unit(-0.1, "cm")
    ) +
    geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
    geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
    guides(
      fill = guide_legend(order = 1), 
      shape = guide_legend(order = 2), 
      color = guide_legend(order = 3)
    ) +
    scale_shape_manual(values = c(15, 16, 17, 18, 6)) +
    scale_color_manual(values = c("#3162A0", "#3FBDA7", "#F08080", "#F0C986")) +
    geom_point(aes(color = treatment, shape = Niche), 
               size = 2.5, alpha = 0.8) +
    labs(
      x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'),
      y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%'),
      title = microbe_name
    ) +
    annotate('text', 
             label = annotation_text, 
             x = min(sample_site$PCoA1) + (max(sample_site$PCoA1) - min(sample_site$PCoA1)) * 0.02,
             y = max(sample_site$PCoA2) * 0.95,
             hjust = 0, vjust = 1,
             size = 3, 
             colour = 'black',
             family = "sans")  # Changed from "serif" to "sans"
  
  # Save plot
  ggsave(file.path(results_dir, paste0(microbe_name, "_PCoA.pdf")), 
         plot = p, device = "pdf", dpi = 300, width = 4, height = 3)
  
  # Save data
  write.csv(sample_site, file.path(results_dir, paste0(microbe_name, "_PCoA_coordinates.csv")), 
            row.names = FALSE)
  
  # Return results
  results <- list(
    microbe = microbe_name,
    plot = p,
    data = sample_site,
    stats = data.frame(
      Microbe = microbe_name,
      R2_niche = r2_niche,
      p_niche = p_niche_formatted,
      R2_treatment = r2_treatment,
      p_treatment = p_treatment_formatted,
      NNE_index = nne_index
    )
  )
  
  # Print results
  cat(paste0("Niche effect: R² = ", r2_niche, ", p = ", p_niche_formatted, "\n"))
  cat(paste0("Treatment effect: R² = ", r2_treatment, ", p = ", p_treatment_formatted, "\n"))
  cat(paste0("NNE Index = ", nne_index, "\n"))
  cat(paste("Plot saved as:", paste0(microbe_name, "_PCoA.pdf"), "\n"))
  
  return(results)
}

# Analyze all microbes
microbes <- c("Bacteria", "AMF", "Rhizobia")
all_results <- list()

for (microbe in microbes) {
  results <- perform_pcoa_analysis(microbe)
  all_results[[microbe]] <- results
}

# Combine all statistics into one table
combined_stats <- do.call(rbind, lapply(all_results, function(x) x$stats))
write.csv(combined_stats, file.path(results_dir, "PCoA_combined_statistics.csv"), 
          row.names = FALSE)

# Create summary table
cat(paste0("\n", strrep("=", 50), "\n"))
cat("SUMMARY OF ALL PCoA ANALYSES\n")
cat(strrep("=", 50), "\n\n")

print(combined_stats)

# Create a comparison plot (optional)
if (length(all_results) == 3) {
  # Extract plots
  p_bacteria <- all_results[["Bacteria"]]$plot
  p_amf <- all_results[["AMF"]]$plot
  p_rhizobia <- all_results[["Rhizobia"]]$plot
  
  # Remove titles from individual plots
  p_bacteria <- p_bacteria + labs(title = NULL)
  p_amf <- p_amf + labs(title = NULL)
  p_rhizobia <- p_rhizobia + labs(title = NULL)
  
  # Combine plots
  library(patchwork)
  combined_plot <- (p_bacteria | p_amf | p_rhizobia) + 
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 12, face = "bold", family = "sans"))  # Added sans
  
  # Save combined plot
  ggsave(file.path(results_dir, "All_microbes_PCoA_combined.pdf"), 
         plot = combined_plot, device = "pdf", dpi = 300, 
         width = 12.6, height = 3)  # 3 plots * 4.2 inches each
}

cat("\n", strrep("=", 50), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 50), "\n")
cat(paste("Results saved in:", results_dir, "\n"))
cat("\nFiles generated:\n")
cat("1. Bacteria_PCoA.pdf & Bacteria_PCoA_coordinates.csv\n")
cat("2. AMF_PCoA.pdf & AMF_PCoA_coordinates.csv\n")
cat("3. Rhizobia_PCoA.pdf & Rhizobia_PCoA_coordinates.csv\n")
cat("4. PCoA_combined_statistics.csv\n")
cat("5. All_microbes_PCoA_combined.pdf (optional combined plot)\n")

# Session info
sessionInfo()

