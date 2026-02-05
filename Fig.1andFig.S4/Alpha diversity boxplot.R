# =============================================================================
# File: 02_Alpha_diversity_analysis.R
# Description: Statistical analysis and visualization of alpha diversity (observed OTUs only)
# =============================================================================

# Clean workspace
rm(list = ls())

# Load required packages
library(vegan)
library(reshape2)
library(multcomp)
library(ggplot2)
library(dplyr)

# Set random seed for reproducibility
set.seed(123)

# Set paths
data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"

# 1. PROCESS BACTERIA DATA ====================================================

message("Processing Bacteria...")

# Load bacteria alpha diversity data
bacteria_alpha <- read.csv(file.path(results_dir, "Bacteria_alpha.csv"), 
                           header = TRUE, row.names = 1)
bacteria_alpha <- bacteria_alpha[grepl("N0|N2|N4|N6",rownames(bacteria_alpha)),]#######################加上这行########################
bacteria_alpha$sample <- rownames(bacteria_alpha)

# Load group data
group_data <- read.csv(file.path(data_dir, "group.csv"), header = TRUE)

# Merge bacteria data with group data
bacteria_data <- merge(bacteria_alpha, group_data, by = "sample")

# Select only observed_otus and required columns
bacteria_data <- bacteria_data[, c("sample", "group", "Niche", "observed_otus")]

# Add niche abbreviations and Group_Niche
bacteria_data <- bacteria_data %>%
  mutate(
    Niche_abbr = case_when(
      Niche == "Original" ~ "O",
      Niche == "Rootzone" ~ "Z", 
      Niche == "Rhizosphere" ~ "RH",
      Niche == "Root" ~ "R",
      Niche == "Nodule" ~ "N",
      TRUE ~ as.character(Niche)
    ),
    Group_Niche = paste(group, Niche_abbr, sep = ".")
  )

# Define order for plotting
group_niche_order <- c()
for (niche in c("O", "Z", "RH", "RT", "N")) {
  for (treatment in c("N0", "N2", "N4", "N6")) {
    group_niche_order <- c(group_niche_order, paste(treatment, niche, sep = "."))
  }
}

bacteria_data$Group_Niche <- factor(bacteria_data$Group_Niche, levels = group_niche_order)

# Calculate ANOVA and Tukey HSD for each niche
niches <- c("Original", "Rootzone", "Rhizosphere", "Root", "Nodule")
bacteria_stats <- NULL

for (niche in niches) {
  niche_data <- subset(bacteria_data, Niche == niche)
  
  if (nrow(niche_data) > 0) {
    # Prepare data for ANOVA
    dat <- niche_data[, c("group", "observed_otus")]
    names(dat) <- c('class', 'var')
    dat$class <- factor(dat$class)
    
    # One-way ANOVA
    fit <- aov(var ~ class, dat)
    p_value <- summary(fit)[[1]][1, 5]
    
    # Tukey HSD test
    tuk <- cld(glht(fit, alternative = 'two.sided', 
                    linfct = mcp(class = 'Tukey')), 
               decreasing = TRUE)
    
    sig <- data.frame(tuk$mcletters$Letters, stringsAsFactors = FALSE)
    names(sig) <- 'sig'
    sig$class <- rownames(sig)
    
    # Calculate mean and SD
    dat_summary <- dat %>% 
      group_by(class) %>% 
      summarise(mean = mean(var), sd = sd(var))
    
    # Merge with significance
    niche_stats <- merge(dat_summary, sig, by = "class")
    niche_stats$Niche <- niche
    niche_stats$Microbe <- "Bacteria"
    
    bacteria_stats <- rbind(bacteria_stats, niche_stats)
  }
}

# Get max values for significance positioning
bacteria_max <- bacteria_data %>%
  group_by(group) %>%
  summarise(max_value = max(observed_otus, na.rm = TRUE))

bacteria_overall_max <- max(bacteria_data$observed_otus, na.rm = TRUE)

# Merge stats with max values
bacteria_stats <- bacteria_stats %>%
  left_join(bacteria_max, by = c("class" = "group"))

bacteria_stats$overall_max <- bacteria_overall_max

# 2. PROCESS AMF DATA =========================================================

message("Processing AMF...")

# Load AMF alpha diversity data
amf_alpha <- read.csv(file.path(results_dir, "AMF_alpha.csv"), 
                      header = TRUE, row.names = 1)
amf_alpha <- amf_alpha[grepl("N0|N2|N4|N6",rownames(amf_alpha)),]
amf_alpha$sample <- rownames(amf_alpha)

# Merge AMF data with group data
amf_data <- merge(amf_alpha, group_data, by = "sample")

# Select only observed_otus and required columns
amf_data <- amf_data[grepl("N0|N2|N4|N6",amf_alpha$sample), c("sample", "group", "Niche", "observed_otus")]

# Add niche abbreviations and Group_Niche
amf_data <- amf_data %>%
  mutate(
    Niche_abbr = case_when(
      Niche == "Original" ~ "O",
      Niche == "Rootzone" ~ "Z", 
      Niche == "Rhizosphere" ~ "RH",
      Niche == "Root" ~ "RT",
      Niche == "Nodule" ~ "N",
      TRUE ~ as.character(Niche)
    ),
    Group_Niche = paste(group, Niche_abbr, sep = ".")
  )

amf_data$Group_Niche <- factor(amf_data$Group_Niche, levels = group_niche_order)

# Calculate ANOVA and Tukey HSD for each niche
amf_stats <- NULL

for (niche in niches) {
  niche_data <- subset(amf_data, Niche == niche)
  
  if (nrow(niche_data) > 0) {
    # Prepare data for ANOVA
    dat <- niche_data[, c("group", "observed_otus")]
    names(dat) <- c('class', 'var')
    dat$class <- factor(dat$class)
    
    # One-way ANOVA
    fit <- aov(var ~ class, dat)
    p_value <- summary(fit)[[1]][1, 5]
    
    # Tukey HSD test
    tuk <- cld(glht(fit, alternative = 'two.sided', 
                    linfct = mcp(class = 'Tukey')), 
               decreasing = TRUE)
    
    sig <- data.frame(tuk$mcletters$Letters, stringsAsFactors = FALSE)
    names(sig) <- 'sig'
    sig$class <- rownames(sig)
    
    # Calculate mean and SD
    dat_summary <- dat %>% 
      group_by(class) %>% 
      summarise(mean = mean(var), sd = sd(var))
    
    # Merge with significance
    niche_stats <- merge(dat_summary, sig, by = "class")
    niche_stats$Niche <- niche
    niche_stats$Microbe <- "AMF"
    
    amf_stats <- rbind(amf_stats, niche_stats)
  }
}

# Get max values for significance positioning
amf_max <- amf_data %>%
  group_by(group) %>%
  summarise(max_value = max(observed_otus, na.rm = TRUE))

amf_overall_max <- max(amf_data$observed_otus, na.rm = TRUE)

# Merge stats with max values
amf_stats <- amf_stats %>%
  left_join(amf_max, by = c("class" = "group"))

amf_stats$overall_max <- amf_overall_max

# 3. PROCESS RHIZOBIA DATA ====================================================

message("Processing Rhizobia...")

# Load Rhizobia alpha diversity data
rhizobia_alpha <- read.csv(file.path(results_dir, "Rhizobia_alpha.csv"), 
                           header = TRUE, row.names = 1)
rhizobia_alpha <- rhizobia_alpha[grepl("N0|N2|N4|N6",rownames(rhizobia_alpha)),]
rhizobia_alpha$sample <- rownames(rhizobia_alpha)

# Merge Rhizobia data with group data
rhizobia_data <- merge(rhizobia_alpha, group_data, by = "sample")

# Select only observed_otus and required columns
rhizobia_data <- rhizobia_data[grepl("N0|N2|N4|N6",rhizobia_alpha$sample), c("sample", "group", "Niche", "observed_otus")]

# Add niche abbreviations and Group_Niche
rhizobia_data <- rhizobia_data %>%
  mutate(
    Niche_abbr = case_when(
      Niche == "Original" ~ "O",
      Niche == "Rootzone" ~ "Z", 
      Niche == "Rhizosphere" ~ "RH",
      Niche == "Root" ~ "RT",
      Niche == "Nodule" ~ "N",
      TRUE ~ as.character(Niche)
    ),
    Group_Niche = paste(group, Niche_abbr, sep = ".")
  )

rhizobia_data$Group_Niche <- factor(rhizobia_data$Group_Niche, levels = group_niche_order)

# Calculate ANOVA and Tukey HSD for each niche
rhizobia_stats <- NULL

for (niche in niches) {
  niche_data <- subset(rhizobia_data, Niche == niche)
  
  if (nrow(niche_data) > 0) {
    # Prepare data for ANOVA
    dat <- niche_data[, c("group", "observed_otus")]
    names(dat) <- c('class', 'var')
    dat$class <- factor(dat$class)
    
    # One-way ANOVA
    fit <- aov(var ~ class, dat)
    p_value <- summary(fit)[[1]][1, 5]
    
    # Tukey HSD test
    tuk <- cld(glht(fit, alternative = 'two.sided', 
                    linfct = mcp(class = 'Tukey')), 
               decreasing = TRUE)
    
    sig <- data.frame(tuk$mcletters$Letters, stringsAsFactors = FALSE)
    names(sig) <- 'sig'
    sig$class <- rownames(sig)
    
    # Calculate mean and SD
    dat_summary <- dat %>% 
      group_by(class) %>% 
      summarise(mean = mean(var), sd = sd(var))
    
    # Merge with significance
    niche_stats <- merge(dat_summary, sig, by = "class")
    niche_stats$Niche <- niche
    niche_stats$Microbe <- "Rhizobia"
    
    rhizobia_stats <- rbind(rhizobia_stats, niche_stats)
  }
}

# Get max values for significance positioning
rhizobia_max <- rhizobia_data %>%
  group_by(group) %>%
  summarise(max_value = max(observed_otus, na.rm = TRUE))

rhizobia_overall_max <- max(rhizobia_data$observed_otus, na.rm = TRUE)

# Merge stats with max values
rhizobia_stats <- rhizobia_stats %>%
  left_join(rhizobia_max, by = c("class" = "group"))

rhizobia_stats$overall_max <- rhizobia_overall_max

# 4. CREATE COMBINED STATS TABLE ==============================================

combined_stats <- rbind(bacteria_stats, amf_stats, rhizobia_stats)
write.csv(combined_stats, file.path(results_dir, "combined_alpha_stats.csv"), 
          row.names = FALSE)

# 5. CREATE PLOTS FOR EACH MICROBE ============================================

# Colors for treatments (20 colors as in original code)
colors <- c("#3162A0", "#3FBDA7", "#F08080", "#F0C986",
            "#3162A0", "#3FBDA7", "#F08080", "#F0C986",
            "#3162A0", "#3FBDA7", "#F08080", "#F0C986",
            "#3162A0", "#3FBDA7", "#F08080", "#F0C986",
            "#3162A0", "#3FBDA7", "#F08080", "#F0C986")

# 5.1 Bacteria plot
message("Creating Bacteria plot...")
bacteria_data$group <- factor(bacteria_data$group,levels = c("N0.O" ,"N2.O", "N4.O" ,"N6.O",
                                                             "N0.Z" ,"N2.Z", "N4.Z" ,"N6.Z",
                                                             "N0.RH" ,"N2.RH", "N4.RH" ,"N6.RH",
                                                             "N0.RT" ,"N2.RT", "N4.RT" ,"N6.RT",
                                                             "N0.N" ,"N2.N", "N4.N" ,"N6.N"))
p_bacteria <- ggplot(bacteria_data, aes(x = group, y = observed_otus)) +
  # Gray background rectangles
  geom_rect(aes(xmin = "N0.Z", xmax = "N6.Z", ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  geom_rect(aes(xmin = 12.5, xmax = 16.5, ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  
  # Boxplot and jitter points
  geom_boxplot(aes(color = group), width = 0.6) +
  geom_jitter(aes(x = group, y = observed_otus, color = group), 
              width = 0.01, size = 1.5) +
  
  # Significance labels
  geom_text(data = bacteria_stats, 
            aes(x = class, label = sig, y = max_value + overall_max * 0.1), 
            size = 4) +
  
  # Color scales
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  
  # Labels and theme
  labs(x = NULL, y = "Richness") +
  ggtitle("Bacteria") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", 
                              face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold",
                               color = "#333333", size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif", size = 11)
  )
p_bacteria
 ggsave(file.path(results_dir, "Bacteria_richness_boxplot.pdf"), 
       plot = p_bacteria, device = "pdf", width = 9, height = 4, dpi = 300)

# 5.2 AMF plot
message("Creating AMF plot...")
amf_data$group <- factor(amf_data$group,levels = c("N0.O" ,"N2.O", "N4.O" ,"N6.O",
                                                              "N0.Z" ,"N2.Z", "N4.Z" ,"N6.Z",
                                                              "N0.RH" ,"N2.RH", "N4.RH" ,"N6.RH",
                                                              "N0.RT" ,"N2.RT", "N4.RT" ,"N6.RT",
                                                              "N0.N" ,"N2.N", "N4.N" ,"N6.N"))
p_amf <- ggplot(amf_data, aes(x = group, y = observed_otus)) +
  # Gray background rectangles
  geom_rect(aes(xmin = "N0.Z", xmax = "N6.Z", ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  geom_rect(aes(xmin = 12.5, xmax = 16.5, ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  
  # Boxplot and jitter points
  geom_boxplot(aes(color = group), width = 0.6) +
  geom_jitter(aes(x = group, y = observed_otus, color = group), 
              width = 0.01, size = 1.5) +
  
  # Significance labels
  geom_text(data = amf_stats, 
            aes(x = class, label = sig, y = max_value + overall_max * 0.1), 
            size = 4) +
  
  # Color scales
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  
  # Labels and theme
  labs(x = NULL, y = "Richness") +
  ggtitle("AMF") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", 
                              face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold",
                               color = "#333333", size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif", size = 11)
  )
p_amf
ggsave(file.path(results_dir, "AMF_Richness_boxplot.pdf"), 
       plot = p_amf, device = "pdf", width = 6, height = 4, dpi = 300)

# 5.3 Rhizobia plot
message("Creating Rhizobia plot...")
rhizobia_data$group <- factor(rhizobia_data$group,levels = c("N0.O" ,"N2.O", "N4.O" ,"N6.O",
                                                   "N0.Z" ,"N2.Z", "N4.Z" ,"N6.Z",
                                                   "N0.RH" ,"N2.RH", "N4.RH" ,"N6.RH",
                                                   "N0.RT" ,"N2.RT", "N4.RT" ,"N6.RT",
                                                   "N0.N" ,"N2.N", "N4.N" ,"N6.N"))
p_rhizobia <- ggplot(rhizobia_data, aes(x = group, y = observed_otus)) +
  # Gray background rectangles
  geom_rect(aes(xmin = "N0.Z", xmax = "N6.Z", ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  geom_rect(aes(xmin = 12.5, xmax = 16.5, ymin = -Inf, ymax = Inf), 
            fill = "#ECECEC", alpha = 0.3) +
  
  # Boxplot and jitter points
  geom_boxplot(aes(color = group), width = 0.6) +
  geom_jitter(aes(x = group, y = observed_otus, color = group), 
              width = 0.01, size = 1.5) +
  
  # Significance labels
  geom_text(data = rhizobia_stats, 
            aes(x = class, label = sig, y = max_value + overall_max * 0.1), 
            size = 4) +
  
  # Color scales
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  
  # Labels and theme
  labs(x = NULL, y = "Richness") +
  ggtitle("Rhizobia") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", 
                              face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold",
                               color = "#333333", size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif", size = 11)
  )
p_rhizobia
ggsave(file.path(results_dir, "Rhizobia_Richness_boxplot.pdf"), 
       plot = p_rhizobia, device = "pdf", width = 6, height = 4, dpi = 300)

message("\nAnalysis completed successfully!")
message(paste("All plots saved in:", results_dir))

# Session info
sessionInfo()
