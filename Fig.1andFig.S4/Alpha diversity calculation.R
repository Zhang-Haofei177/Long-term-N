# =============================================================================
# File: 01_Alpha_diversity_analysis.R
# Description: Alpha diversity calculation for Bacteria, AMF, and Rhizobia
# =============================================================================

# Clean workspace
rm(list = ls())

# Load required packages
library(vegan)

# Set random seed for reproducibility
set.seed(123)

# Set paths
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

# Create output directory
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Alpha diversity calculation function
ly.alpha <- function(a, method = "x") {
  names <- c("observed_otus", "shannon", "simpson", "chao1", "ace", "goods_coverage", "pielou.evenness")
  TF <- method %in% names
  if (FALSE %in% TF == TRUE) {
    print("Error: Please input correct alpha diversity index names")
  } else {
    alpha.1 = NULL
    if ("observed_otus" %in% method == TRUE) {
      observed_otus <- rowSums(a > 0)
      alpha.1 <- cbind(alpha.1, observed_otus)
    }
    if ("shannon" %in% method == TRUE) {
      shannon <- {a / as.vector(rowSums(a))} * {log(a / as.vector(rowSums(a)), 2)}
      shannon[is.na(shannon)] <- 0
      shannon <- {-rowSums(shannon)}
      alpha.1 <- cbind(alpha.1, shannon)
    }
    if ("simpson" %in% method == TRUE) {
      simpson <- diversity(a, index = "simpson")
      alpha.1 <- cbind(alpha.1, simpson)
    }
    if ("chao1" %in% method == TRUE) {
      observed_otus <- rowSums(a > 0)
      chao1 <- observed_otus + ((rowSums(a == 1)) * (rowSums(a == 1) - 1)) / (2 * (rowSums(a == 2) + 1))
      alpha.1 <- cbind(alpha.1, chao1)
    }
    if ("ace" %in% method == TRUE) {
      n.rare <- 0
      for (i in 1:10) {
        n.rare = n.rare + i * rowSums(a == i)
      }
      c.ace <- 1 - rowSums(a == 1) / n.rare
      j <- 0
      for (i in 1:10) {
        j = j + i * (i - 1) * rowSums(a == i)
      }
      max <- rowSums(a <= 10 & a > 0) * j / (c.ace * n.rare * (n.rare - 1)) - 1
      max[max < 0] <- 0
      ace <- rowSums(a > 10) + rowSums(a <= 10 & a > 0) / c.ace + max * rowSums(a == 1) / c.ace
      alpha.1 <- cbind(alpha.1, ace)
    }
    if ("goods_coverage" %in% method == TRUE) {
      goods_coverage <- 1 - rowSums(a == 1) / (as.vector(rowSums(a))[1])
      alpha.1 <- cbind(alpha.1, goods_coverage)
    }
    if ("pielou.evenness" %in% method == TRUE) {
      observed_otus <- rowSums(a > 0)
      shannon <- {a / as.vector(rowSums(a))} * {log(a / as.vector(rowSums(a)), 2)}
      shannon[is.na(shannon)] <- 0
      shannon <- {-rowSums(shannon)}
      pielou.evenness <- shannon / log(observed_otus)
      alpha.1 <- cbind(alpha.1, pielou.evenness)
    }
    alpha.1
  }
}

# Calculate alpha diversity indices
diversity_indices <- c("observed_otus", "shannon", "simpson", "chao1", "ace", "goods_coverage", "pielou.evenness")

# 1. Bacteria
bacteria_data <- read.csv(file.path(data_dir, "Bacteria_ASV_rarefaction.csv"), row.names = 1)
bacteria_alpha <- ly.alpha(t(bacteria_data), diversity_indices)
write.csv(bacteria_alpha, file.path(results_dir, "Bacteria_alpha.csv"))

# 2. AMF
amf_data <- read.csv(file.path(data_dir, "AMF_ASV_rarefaction.csv"), row.names = 1)
amf_alpha <- ly.alpha(t(amf_data), diversity_indices)
write.csv(amf_alpha, file.path(results_dir, "AMF_alpha.csv"))

# 3. Rhizobia
rhizobia_data <- read.csv(file.path(data_dir, "Rhizobia_ASV_rarefaction.csv"), row.names = 1)
rhizobia_alpha <- ly.alpha(t(rhizobia_data), diversity_indices)
write.csv(rhizobia_alpha, file.path(results_dir, "Rhizobia_alpha.csv"))

# Session info
sessionInfo()
