setup_environment <- function() {
  # Clean workspace
  rm(list = ls())
  
  # List of required packages
  required_packages <- c(
    "linkET",      # For Mantel test
    "dplyr",       # For data manipulation
    "tidyverse",   # For data wrangling and visualization
    "reshape2",    # For data reshaping
    "Hmisc",       # For Spearman correlation
    "ggplot2",     # For plotting
    "RColorBrewer",# For color palettes
    "patchwork"    # For combining plots
  )
  
  # Install missing packages
  missing_packages <- required_packages[!required_packages %in% installed.packages()]
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, dependencies = TRUE)
  }
  
  # Load packages
  suppressPackageStartupMessages({
    library(linkET)
    library(dplyr)
    library(tidyverse)
    library(reshape2)
    library(Hmisc)
    library(ggplot2)
    library(RColorBrewer)
    library(patchwork)
  })
  
  # Set paths (modify these according to your directory structure)
  data_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
  results_dir <- "C:/Users/ADMIN/Desktop/Longterm N code/Fig.1andFig.S4/"
  
  # Create results directory if it doesn't exist
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Return paths
  return(list(data_dir = data_dir, results_dir = results_dir))
}

# =============================================================================
# SECTION 2: DATA PREPARATION
# =============================================================================

#' @title Load and Prepare Data
#' @description Load microbial and environmental data, and prepare for analysis
load_and_prepare_data <- function(data_dir) {
  cat("Step 1: Loading and preparing data...\n")
  
  # Load soil properties data
  trait <- read.csv(file.path(data_dir, "Mantel_env.csv"), header = TRUE)
  
  # Load microbial data
  AMF <- read.csv(file.path(data_dir, "AMF_ASV_rarefaction.csv"), 
                  header = TRUE, row.names = 1)
  Rhizobia <- read.csv(file.path(data_dir, "Rhizobia_ASV_rarefaction.csv"), 
                       header = TRUE, row.names = 1)
  
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
  
  # Split data by niche
  niches <- list(
    ORI = TASV[grepl("\\.O\\.", rownames(TASV)), ],      # Original soil
    BS  = TASV[grepl("\\.Z\\.", rownames(TASV)), ],      # Bulk soil (Rootzone)
    RH  = TASV[grepl("\\.RH\\.", rownames(TASV)), ],     # Rhizosphere
    RT  = TASV[grepl("\\.RT\\.", rownames(TASV)), ],     # Root
    NOD = TASV[grepl("\\.N\\.", rownames(TASV)), ]       # Nodule
  )
  
  # Print sample counts
  cat("\nSample counts by niche:\n")
  for (niche_name in names(niches)) {
    cat(paste0(niche_name, ": ", nrow(niches[[niche_name]]), " samples\n"))
  }
  
  return(list(
    trait = trait,
    TAMF = TAMF,
    TRhizobia = TRhizobia,
    niches = niches
  ))
}

# =============================================================================
# SECTION 3: MANTEL ANALYSIS
# =============================================================================

#' @title Perform Mantel Analysis
#' @description Calculate Mantel correlations between microbial communities and environmental factors
perform_mantel_analysis <- function(niche_name, micro_data, trait_data, TAMF, TRhizobia, results_dir) {
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
  
  # Filter trait data for the niche
  trait_niche <- trait_data[grepl(pattern, rownames(trait_data)), ]
  
  # Ensure samples match
  common_samples <- intersect(rownames(micro_data), rownames(trait_niche))
  
  if (length(common_samples) == 0) {
    cat("No common samples found.\n")
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
  
  cat("Mantel results saved as:", results_file, "\n")
  
  # Print summary
  sig_results <- mantel_result[mantel_result$p < 0.05, ]
  if (nrow(sig_results) > 0) {
    cat("\nSignificant correlations (p < 0.05):\n")
    for (i in 1:nrow(sig_results)) {
      cat(paste0(
        sig_results$spec[i], " - ", sig_results$env[i], 
        ": r = ", round(sig_results$r[i], 3), 
        ", p = ", ifelse(sig_results$p[i] < 0.001, "<0.001", round(sig_results$p[i], 3)), "\n"
      ))
    }
  } else {
    cat("\nNo significant correlations found.\n")
  }
  
  return(mantel_result)
}

#' @title Run Mantel Analysis for All Niches
run_all_mantel_analyses <- function(data_prep, results_dir) {
  cat("\n", strrep("=", 60), "\n")
  cat("PERFORMING MANTEL ANALYSES FOR ALL NICHES\n")
  cat(strrep("=", 60), "\n")
  
  all_results <- list()
  
  for (niche_name in names(data_prep$niches)) {
    if (nrow(data_prep$niches[[niche_name]]) > 0) {
      results <- perform_mantel_analysis(
        niche_name = niche_name,
        micro_data = data_prep$niches[[niche_name]],
        trait_data = data_prep$trait,
        TAMF = data_prep$TAMF,
        TRhizobia = data_prep$TRhizobia,
        results_dir = results_dir
      )
      if (!is.null(results)) {
        all_results[[niche_name]] <- results
      }
    } else {
      cat(paste("\nNo microbial data for", niche_name, "niche\n"))
    }
  }
  
  return(all_results)
}

# =============================================================================
# SECTION 4: COMBINE AND PROCESS MANTEL RESULTS
# =============================================================================

#' @title Combine Mantel Results
combine_mantel_results <- function(results_dir) {
  cat("\n", strrep("=", 60), "\n")
  cat("COMBINING MANTEL RESULTS\n")
  cat(strrep("=", 60), "\n")
  
  # Check which files were created
  created_files <- list.files(results_dir, pattern = "Mantel_results_.*\\.csv", full.names = TRUE)
  
  if (length(created_files) == 0) {
    cat("No Mantel result files found.\n")
    return(NULL)
  }
  
  # Read all Mantel result files
  mantel_data <- list()
  for (file in created_files) {
    niche_name <- gsub("Mantel_results_|\\.csv", "", basename(file))
    cat("Reading", niche_name, "...\n")
    mantel_data[[niche_name]] <- read.csv(file, header = TRUE)
  }
  
  # Define niche abbreviations
  niche_abbrev <- list(
    "ORI" = "MFS",    # Original soil -> MFS
    "BS" = "RZS",     # Bulk soil -> RZS
    "RH" = "RS",      # Rhizosphere -> RS
    "RT" = "RE",      # Root -> RE
    "NOD" = "NE"      # Nodule -> NE
  )
  
  # Combine all results
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
  
  if (length(total_list) == 0) {
    cat("No data to combine.\n")
    return(NULL)
  }
  
  total <- do.call(rbind, total_list)
  
  # Separate r and p values
  total_r <- total[c("id", "niche", "env", "r")]
  total_p <- total[c("id", "niche", "env", "p")]
  
  # Convert p-values to significance symbols
  total_p <- total_p %>%
    mutate(p = ifelse(p < 0.001, '***',
                      ifelse(p >= 0.001 & p < 0.01, '**',
                             ifelse(p >= 0.01 & p < 0.05, '*', ''))))
  
  # Convert to wide format
  r_width <- total_r %>% 
    pivot_wider(names_from = env, values_from = r)
  
  p_width <- total_p %>% 
    pivot_wider(names_from = env, values_from = p)
  
  r_width <- as.data.frame(r_width)
  p_width <- as.data.frame(p_width)
  
  # Set row names
  rownames(r_width) <- r_width$id
  rownames(p_width) <- p_width$id
  
  # Order rows by microbial group
  AMF_NAME <- r_width$id[grepl("AMF", r_width$id)]
  Rhi_NAME <- r_width$id[grepl("Rhizobia", r_width$id)]
  
  if (length(AMF_NAME) > 0 && length(Rhi_NAME) > 0) {
    order <- c(Rhi_NAME, AMF_NAME)
    r_width <- r_width[order, ]
    p_width <- p_width[order, ]
  }
  
  # Extract matrices for plotting
  if (ncol(r_width) > 3) {
    r_plot <- r_width[4:ncol(r_width)]
    p_plot <- p_width[4:ncol(p_width)]
    
    # Transpose for final format
    r_plot <- t(r_plot)
    p_plot <- t(p_plot)
    
    # Save results
    write.csv(r_plot, file.path(results_dir, "AMF_Rhizobia_mantel_R.csv"))
    write.csv(p_plot, file.path(results_dir, "AMF_Rhizobia_mantel_p.csv"))
    
    cat("\nCombined results saved:\n")
    cat("1. AMF_Rhizobia_mantel_R.csv\n")
    cat("2. AMF_Rhizobia_mantel_p.csv\n")
    
    # Create summary statistics
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
    
    cat("3. Mantel_significant_summary.csv\n")
    
    return(list(
      r_matrix = r_plot,
      p_matrix = p_plot,
      summary = summary_stats
    ))
  }
  
  return(NULL)
}

# =============================================================================
# SECTION 5: SOIL PROPERTIES CORRELATION ANALYSIS
# =============================================================================

#' @title Calculate Soil Properties Correlations
analyze_soil_correlations <- function(data_dir, results_dir) {
  cat("\n", strrep("=", 60), "\n")
  cat("ANALYZING SOIL PROPERTIES CORRELATIONS\n")
  cat(strrep("=", 60), "\n")
  
  # Load soil properties data
  trait <- read.csv(file.path(data_dir, "Mantel_env.csv"), header = TRUE)
  
  # Prepare soil properties data
  trait <- trait[, c(1, 5:25)]
  rownames(trait) <- trait$sample
  
  # Filter for Original niche (samples containing ".O.")
  original_trait <- trait[grep("\\.O\\.", trait$sample), ]
  rownames(original_trait) <- original_trait$sample
  original_trait$sample <- NULL
  
  # Calculate Spearman correlation with p-values
  cor_result <- rcorr(as.matrix(original_trait), type = "spearman")
  
  # Extract correlation coefficients and p-values
  cor_r <- cor_result$r
  cor_p <- cor_result$P
  
  # Replace NA values
  diag(cor_r) <- 1
  diag(cor_p) <- 0
  
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
      P_value_formatted = ifelse(
        P_value < 0.001,
        sprintf("%.2E", P_value),
        sprintf("%.8f", P_value)
      ),
      r = round(r, 8)
    ) %>%
    arrange(tax, Index) %>%
    select(tax, Index, r, P_value = P_value_formatted, P_value_sig)
  
  # Save results
  output_file <- file.path(results_dir, "soil_Spearman_correlation.csv")
  write.csv(cor_combined, output_file, row.names = FALSE, quote = FALSE)
  
  cat("Soil properties correlation saved as:", output_file, "\n")
  
  return(cor_combined)
}

# =============================================================================
# SECTION 6: CREATE MANTEL LINK TABLES
# =============================================================================

#' @title Create Mantel Link Tables
create_mantel_link_tables <- function(results_dir) {
  cat("\n", strrep("=", 60), "\n")
  cat("CREATING MANTEL LINK TABLES\n")
  cat(strrep("=", 60), "\n")
  
  # Load required data
  bio_r <- read.csv(file.path(results_dir, "AMF_Rhizobia_mantel_R.csv"), header = TRUE)
  bio_p <- read.csv(file.path(results_dir, "AMF_Rhizobia_mantel_p.csv"), header = TRUE)
  
  # Load line site data (coordinates for plotting)
  # Note: These files need to be created separately with the coordinate information
  AMF_line_site <- read.csv(file.path(results_dir, "AMF_line_site.csv"))
  Rhizobia_line_site <- read.csv(file.path(results_dir, "Rhizobia_line_site.csv"))
  
  # Reshape bio data to long format
  bio_r_long <- melt(bio_r, id.vars = "X", variable.name = "env", value.name = "r")
  colnames(bio_r_long)[1] <- "bio"
  
  bio_p_long <- melt(bio_p, id.vars = "X", variable.name = "env", value.name = "p")
  colnames(bio_p_long)[1] <- "bio"
  
  # Merge with line site data
  create_mantel_line <- function(line_site, bio_r_long, bio_p_long, microbe_name) {
    line_data <- merge(line_site, bio_r_long, by = c("bio", "env"), all.x = TRUE)
    line_data <- merge(line_data, bio_p_long, by = c("bio", "env"), all.x = TRUE)
    
    line_data <- line_data %>%
      mutate(
        r.sign = ifelse(r > 0, "Positive", "Negative"),
        r.abs = cut(abs(r), 
                    breaks = c(0, 0.25, 0.5, 1), 
                    labels = c("<0.25", "0.25-0.5", "0.5-1"),
                    include.lowest = TRUE,
                    right = FALSE),
        p.sign = p
      ) %>%
      select(bio, env, ymin, xmin, xmax, ymax, r, p, r.sign, r.abs, p.sign)
    
    output_file <- file.path(results_dir, paste0(microbe_name, "_mantel_line.csv"))
    write.csv(line_data, output_file, row.names = FALSE)
    
    cat(paste0(microbe_name, " mantel line table saved as: ", output_file, "\n"))
    
    return(line_data)
  }
  
  # Create mantel line tables for both microbes
  AMF_mantel_line <- create_mantel_line(AMF_line_site, bio_r_long, bio_p_long, "AMF")
  Rhizobia_mantel_line <- create_mantel_line(Rhizobia_line_site, bio_r_long, bio_p_long, "Rhizobia")
  
  return(list(
    AMF = AMF_mantel_line,
    Rhizobia = Rhizobia_mantel_line
  ))
}

# =============================================================================
# SECTION 7: VISUALIZATION
# =============================================================================

#' @title Create Environment Correlation Heatmap
create_environment_heatmap <- function(env_data) {
  # Prepare environment correlation data
  env <- env_data[c(1,2,3,4)]
  colnames(env) <- c(".rownames", ".colnames", "r", "p")
  
  # Create correlation matrices
  rownames_unique <- unique(env$.rownames)
  colnames_unique <- unique(env$.colnames)
  
  r_matrix <- matrix(NA, nrow = length(rownames_unique), ncol = length(colnames_unique),
                     dimnames = list(rownames_unique, colnames_unique))
  p_matrix <- matrix(NA, nrow = length(rownames_unique), ncol = length(colnames_unique),
                     dimnames = list(rownames_unique, colnames_unique))
  
  for (i in 1:nrow(env)) {
    r_matrix[env$.rownames[i], env$.colnames[i]] <- env$r[i]
    p_matrix[env$.rownames[i], env$.colnames[i]] <- env$p[i]
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
      TRUE ~ "ns"
    ))
  
  # Create heatmap
  p1 <- ggplot() +
    geom_tile(data = long_cor, aes(Var2, Var1), fill = "white", color = "grey") +
    geom_point(data = long_cor, aes(Var2, Var1, size = abs(value), color = value), shape = 15) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_text(angle = -135, hjust = 1),
      axis.text.y = element_text(angle = -135, hjust = 1),
      legend.text = element_text(angle = 0),
      legend.title = element_text(angle = 0)
    ) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL) +
    scale_color_gradientn(
      colors = rev(brewer.pal(n = 9, name = "RdYlBu")),
      limits = c(-1, 1)
    ) +
    geom_text(
      data = long_p,
      aes(x = Var2, y = Var1, label = stars),
      size = 3, color = "black", angle = -135
    ) +
    coord_fixed(ratio = 1)
  
  return(p1)
}

#' @title Create Mantel Link Plot
create_mantel_link_plot <- function(line_data, title, curvature = 0.1, text_angle = -135) {
  # Ensure p.sign has proper values for plotting
  line_data$p.sign <- ifelse(line_data$p.sign == "", "NA", line_data$p.sign)
  
  plot <- ggplot() +
    geom_curve(
      data = line_data,
      aes(x = xmin, y = ymin, xend = xmax, yend = ymax,
          color = p.sign, size = r.abs),
      curvature = curvature
    ) +
    geom_point(data = line_data, aes(x = xmin, y = ymin)) +
    geom_point(data = line_data, aes(x = xmax, y = ymax)) +
    theme_void() +
    geom_text(
      data = line_data,
      aes(x = xmin, y = ymin, label = env),
      hjust = 1.2, vjust = -0.5, size = 3,
      color = "black", angle = text_angle
    ) +
    geom_text(
      data = line_data,
      aes(x = xmax, y = ymax, label = bio),
      hjust = -0.2, vjust = -0.5, size = 3,
      color = "black", angle = text_angle
    ) +
    scale_size_manual(
      values = c("<0.25" = 0.5, "0.25-0.5" = 1.5, "0.5-1" = 3),
      name = "|r|"
    ) +
    scale_colour_manual(
      values = c("***" = "#283c63", "**" = "#D73027", 
                 "*" = "#4575B7", "NA" = "gray"),
      name = "Significance"
    ) +
    labs(title = title)
  
  return(plot)
}

#' @title Generate All Plots
generate_all_plots <- function(results_dir, env_cor) {
  cat("\n", strrep("=", 60), "\n")
  cat("GENERATING PLOTS\n")
  cat(strrep("=", 60), "\n")
  
  # Load mantel line data
  AMF_line <- read.csv(file.path(results_dir, "AMF_mantel_line.csv"), header = TRUE)
  Rhi_line <- read.csv(file.path(results_dir, "Rhizobia_mantel_line.csv"), header = TRUE)
  
  # Filter data (excluding MFS if needed)
  AMF_line_plot <- AMF_line[grepl("AMF", AMF_line$bio), ]
  AMF_line_plot <- AMF_line_plot[!grepl("MFS", AMF_line_plot$bio), ]
  
  Rhi_line_plot <- Rhi_line[grepl("Rhizobia", Rhi_line$bio), ]
  Rhi_line_plot <- Rhi_line_plot[!grepl("MFS", Rhi_line_plot$bio), ]
  
  # Create environment correlation heatmap
  p1 <- create_environment_heatmap(env_cor)
  
  # Create mantel link plots
  plot_AMF <- create_mantel_link_plot(AMF_line_plot, "AMF", curvature = 0.1, text_angle = -135)
  plot_Rhi <- create_mantel_link_plot(Rhi_line_plot, "Rhizobia", curvature = -0.1, text_angle = -45)
  
  # Combine and save plots
  # 1. AMF combined plot
  p_total_AMF <- p1 | plot_AMF
  ggsave(
    file.path(results_dir, "AMF_mantel.pdf"),
    p_total_AMF, device = "pdf", width = 10, height = 7, dpi = 300
  )
  cat("AMF plot saved as: AMF_mantel.pdf\n")
  
  # 2. Rhizobia combined plot
  p_total_Rhi <- p1 | plot_Rhi
  ggsave(
    file.path(results_dir, "Rhizobia_mantel.pdf"),
    p_total_Rhi, device = "pdf", width = 10, height = 7, dpi = 300
  )
  cat("Rhizobia plot saved as: Rhizobia_mantel.pdf\n")
  
  # 3. Combined all plots
  p_combined_all <- (p1 | plot_AMF | plot_Rhi) +
    plot_layout(widths = c(1, 0.8, 0.8))
  ggsave(
    file.path(results_dir, "All_mantel_combined.pdf"),
    p_combined_all, device = "pdf", width = 15, height = 7, dpi = 300
  )
  cat("Combined plot saved as: All_mantel_combined.pdf\n")
  
  # Print summary
  cat("\nSummary statistics:\n")
  cat("AMF correlations:", nrow(AMF_line_plot), "\n")
  cat("Rhizobia correlations:", nrow(Rhi_line_plot), "\n")
  cat("Environment factors:", length(unique(c(AMF_line_plot$env, Rhi_line_plot$env))), "\n")
  
  return(list(
    heatmap = p1,
    AMF_plot = plot_AMF,
    Rhizobia_plot = plot_Rhi
  ))
}

# =============================================================================
# SECTION 8: MAIN PIPELINE FUNCTION
# =============================================================================

#' @title Main Analysis Pipeline
#' @description Run the complete Mantel analysis pipeline
main_analysis_pipeline <- function() {
  # Record start time
  start_time <- Sys.time()
  cat("Starting Mantel Analysis Pipeline\n")
  cat("Start time:", as.character(start_time), "\n")
  
  # Step 1: Setup environment
  paths <- setup_environment()
  
  # Step 2: Load and prepare data
  data_prep <- load_and_prepare_data(paths$data_dir)
  
  # Step 3: Run Mantel analysis for all niches
  mantel_results <- run_all_mantel_analyses(data_prep, paths$results_dir)
  
  # Step 4: Combine Mantel results
  combined_results <- combine_mantel_results(paths$results_dir)
  
  # Step 5: Analyze soil properties correlations
  soil_cor <- analyze_soil_correlations(paths$data_dir, paths$results_dir)
  
  # Step 6: Create mantel link tables
  # Note: This step requires AMF_line_site.csv and Rhizobia_line_site.csv files
  # These files contain the coordinate information for plotting
  link_tables <- create_mantel_link_tables(paths$results_dir)
  
  # Step 7: Generate plots
  plots <- generate_all_plots(paths$results_dir, soil_cor)
  
  # Record end time and calculate duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  cat("\n", strrep("=", 60), "\n")
  cat("ANALYSIS PIPELINE COMPLETED SUCCESSFULLY\n")
  cat(strrep("=", 60), "\n")
  cat("End time:", as.character(end_time), "\n")
  cat("Total duration:", round(duration, 2), "minutes\n")
  
  # Print session info for reproducibility
  cat("\nSession information:\n")
  print(sessionInfo())
  
  # Return all results
  return(list(
    data_prep = data_prep,
    mantel_results = mantel_results,
    combined_results = combined_results,
    soil_cor = soil_cor,
    link_tables = link_tables,
    plots = plots
  ))
}

# =============================================================================
# SECTION 9: RUN THE PIPELINE
# =============================================================================

# Uncomment the line below to run the complete pipeline
# results <- main_analysis_pipeline()

# =============================================================================
# SECTION 10: HELPER FUNCTIONS FOR GITHUB REPOSITORY
# =============================================================================

#' @title Create README for GitHub Repository
#' @description Generate a README file with instructions for using the pipeline
create_readme <- function() {
  readme_content <- '
# Mantel Analysis Pipeline for AMF and Rhizobia

## Overview
This repository contains R code for performing Mantel correlation analysis between microbial communities (AMF and Rhizobia) and soil properties across different niches. The analysis pipeline includes data preparation, statistical analysis, and visualization.

## Requirements

### R Packages
- linkET
- dplyr
- tidyverse
- reshape2
- Hmisc
- ggplot2
- RColorBrewer
- patchwork

### Data Files
1. `Mantel_env.csv` - Soil properties data
2. `AMF_ASV_rarefaction.csv` - AMF abundance data
3. `Rhizobia_ASV_rarefaction.csv` - Rhizobia abundance data
4. `AMF_line_site.csv` - AMF coordinate data for plotting
5. `Rhizobia_line_site.csv` - Rhizobia coordinate data for plotting

## Usage

1. **Set up directory structure**:
   - Place data files in the appropriate directory
   - Update paths in the `setup_environment()` function

2. **Run the complete pipeline**:
   ```r
   source("mantel_analysis_pipeline.R")
   results <- main_analysis_pipeline()
   