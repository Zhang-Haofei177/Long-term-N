###### Rhizobia_Niche_neutral community model #############
rm(list = ls())
library(Hmisc)
library(minpack.lm)
library(stats4)
library(grid)
library(ggplot2)
library(gridExtra)
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S6/")
spp <- read.csv('Rhizobia_ASV_rarefaction.csv', head = TRUE, stringsAsFactors = FALSE, row.names = 1)
spp <- spp[, grep("N0|N2|N4|N6", colnames(spp))]

# Define data for five compartments
ori <- spp[, grep("O", colnames(spp))]
rz <- spp[, grep("Z", colnames(spp))]
rh <- spp[, grep("RH", colnames(spp))]
rt <- spp[, grep("RT", colnames(spp))]
nod <- spp[, grep("N.1|N.2|N.3|N.4|N.5", colnames(spp))]

# Create lists to store results
plot_list <- list()
results_list <- list()

# Define NCM analysis function
run_ncm_analysis <- function(spp_data, compartment_name) {
  spp <- t(spp_data)
  
  # Calculate mean total relative abundance
  N <- mean(apply(spp, 1, sum))
  
  # Calculate mean relative abundance for each species
  p.m <- apply(spp, 2, mean)
  
  # Remove species with zero mean abundance
  p.m <- p.m[p.m != 0]
  
  # Calculate relative abundance for each species
  p <- p.m / N
  
  # Binarize original data (presence/absence)
  spp.bi <- 1 * (spp > 0)
  
  # Calculate occurrence frequency for each species
  freq <- apply(spp.bi, 2, mean)
  
  # Remove species with zero frequency
  freq <- freq[freq != 0]
  
  # Merge relative abundance and frequency data
  C <- merge(p, freq, by = 0)
  
  # Sort by frequency
  C <- C[order(C[, 2]), ]
  
  # Convert to data frame
  C <- as.data.frame(C)
  
  # Remove rows containing zeros
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
  
  # Extract relative abundance and frequency data
  p <- C.0[, 2]
  freq <- C.0[, 3]
  
  # Assign names to data
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]
  
  # Calculate d value
  d = 1 / N
  
  ## Fit model parameter m (or Nm) using nonlinear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
  
  # Output fitting results
  cat(paste("\n====", compartment_name, "compartment analysis results ====\n"))
  cat("m =", round(coef(m.fit), 4), "\n")
  cat("Nm =", round(coef(m.fit) * N, 2), "\n")
  
  # Calculate confidence interval for m
  m.ci <- confint(m.fit, 'm', level = 0.95)
  cat("95% confidence interval for m:", round(m.ci[1], 4), "-", round(m.ci[2], 4), "\n")
  
  # Calculate predicted frequency values
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
  
  # Calculate confidence intervals using Wilson method
  pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)
  
  # Calculate R-squared value
  Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
  cat("R-squared =", round(Rsqr, 3), "\n")
  
  # Define bacnlsALL data frame containing p, freq, freq.pred, and pred.ci data
  bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])
  
  # Define point colors
  inter.col <- rep('black', nrow(bacnlsALL))
  inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#00a4ac'  # Points below confidence interval
  inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#e06a5d'  # Points above confidence interval
  
  # Create plot
  grid.newpage()
  
  # Set viewport
  pushViewport(viewport(h = 0.6, w = 0.6))
  
  # Set data viewport range
  pushViewport(dataViewport(xData = range(log10(bacnlsALL$p)), yData = c(0, 1.02), extension = c(0.02, 0)))
  
  # Draw rectangle
  grid.rect()
  
  # Draw scatter plot
  grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))
  
  # Add y-axis and x-axis
  grid.yaxis()
  grid.xaxis()
  
  # Draw predicted curve
  grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = '#0a71b4', lwd = 2), default = 'native')
  
  # Draw confidence intervals
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = '#0a71b4', lwd = 2, lty = 2), default = 'native')
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = '#0a71b4', lwd = 2, lty = 2), default = 'native')
  
  # Add text labels
  grid.text(y = unit(0, 'npc') - unit(2.5, 'lines'), label = 'Mean Relative Abundance (log10)', gp = gpar(fontface = 2))
  grid.text(x = unit(0, 'npc') - unit(3, 'lines'), label = 'Frequency of Occurrence', gp = gpar(fontface = 2), rot = 90)
  
  # Add compartment name and statistical results
  grid.text(paste(compartment_name, "\nRsqr=", round(Rsqr, 3), "\nm=", round(coef(m.fit), 4), "\nNm=", round(coef(m.fit) * N, 2)), 
            x = unit(0.8, "npc"), y = unit(0.8, "npc"), just = c("centre", "centre"), gp = gpar(cex = 0.8, fontface = "bold"))
  
  # Capture current plot
  current_plot <- grid.grab()
  
  # Store results
  results <- list(
    compartment = compartment_name,
    m = coef(m.fit),
    Nm = coef(m.fit) * N,
    Rsqr = Rsqr,
    m_ci = m.ci
  )
  
  return(list(plot = current_plot, results = results))
}

# Analyze five compartments separately
compartments <- list(
  ori = ori,
  rz = rz,
  rh = rh,
  rt = rt,
  nod = nod
)

compartment_names <- c("ori", "rz", "rh", "rt", "nod")

# Run analysis and collect results
for (i in 1:length(compartment_names)) {
  name <- compartment_names[i]
  data <- compartments[[name]]
  analysis_result <- run_ncm_analysis(data, name)
  plot_list[[name]] <- analysis_result$plot
  results_list[[name]] <- analysis_result$results
}

# Create combined plot - five plots arranged vertically
grid.newpage()
grid.arrange(
  plot_list$rz,
  plot_list$rh,
  plot_list$rt,
  plot_list$nod,
  ncol = 1,
  heights = c(1, 1, 1, 1, 1)
)

# Save combined plot
ggsave("NCM_Rhizobia_Niche.pdf", 
       arrangeGrob(plot_list$rz, plot_list$rh, plot_list$rt, plot_list$nod, ncol = 1),
       width = 5, height = 15)

###### Rhizobia_Nitrogen_neutral community model #############
rm(list = ls())
library(Hmisc)
library(minpack.lm)
library(stats4)
library(grid)
library(ggplot2)
library(gridExtra)
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S6/")
spp <- read.csv('Rhizobia_ASV_rarefaction.csv', head = TRUE, stringsAsFactors = FALSE, row.names = 1)
spp <- spp[, grep("N0|N2|N4|N6", colnames(spp))]

# Define data for four treatments
N0 <- spp[, grep("N0", colnames(spp))]
N2 <- spp[, grep("N2", colnames(spp))]
N4 <- spp[, grep("N4", colnames(spp))]
N6 <- spp[, grep("N6", colnames(spp))]

# Create lists to store results
plot_list <- list()
results_list <- list()

# Define NCM analysis function
run_ncm_analysis <- function(spp_data, compartment_name) {
  spp <- t(spp_data)
  
  # Calculate mean total relative abundance
  N <- mean(apply(spp, 1, sum))
  
  # Calculate mean relative abundance for each species
  p.m <- apply(spp, 2, mean)
  
  # Remove species with zero mean abundance
  p.m <- p.m[p.m != 0]
  
  # Calculate relative abundance for each species
  p <- p.m / N
  
  # Binarize original data (presence/absence)
  spp.bi <- 1 * (spp > 0)
  
  # Calculate occurrence frequency for each species
  freq <- apply(spp.bi, 2, mean)
  
  # Remove species with zero frequency
  freq <- freq[freq != 0]
  
  # Merge relative abundance and frequency data
  C <- merge(p, freq, by = 0)
  
  # Sort by frequency
  C <- C[order(C[, 2]), ]
  
  # Convert to data frame
  C <- as.data.frame(C)
  
  # Remove rows containing zeros
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
  
  # Extract relative abundance and frequency data
  p <- C.0[, 2]
  freq <- C.0[, 3]
  
  # Assign names to data
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]
  
  # Calculate d value
  d = 1 / N
  
  ## Fit model parameter m (or Nm) using nonlinear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
  
  # Output fitting results
  cat(paste("\n====", compartment_name, "compartment analysis results ====\n"))
  cat("m =", round(coef(m.fit), 4), "\n")
  cat("Nm =", round(coef(m.fit) * N, 2), "\n")
  
  # Calculate confidence interval for m
  m.ci <- confint(m.fit, 'm', level = 0.95)
  cat("95% confidence interval for m:", round(m.ci[1], 4), "-", round(m.ci[2], 4), "\n")
  
  # Calculate predicted frequency values
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
  
  # Calculate confidence intervals using Wilson method
  pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)
  
  # Calculate R-squared value
  Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
  cat("R-squared =", round(Rsqr, 3), "\n")
  
  # Define bacnlsALL data frame containing p, freq, freq.pred, and pred.ci data
  bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])
  
  # Define point colors
  inter.col <- rep('black', nrow(bacnlsALL))
  inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#00a4ac'  # Points below confidence interval
  inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#e06a5d'  # Points above confidence interval
  
  # Create plot
  grid.newpage()
  
  # Set viewport
  pushViewport(viewport(h = 0.6, w = 0.6))
  
  # Set data viewport range
  pushViewport(dataViewport(xData = range(log10(bacnlsALL$p)), yData = c(0, 1.02), extension = c(0.02, 0)))
  
  # Draw rectangle
  grid.rect()
  
  # Draw scatter plot
  grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))
  
  # Add y-axis and x-axis
  grid.yaxis()
  grid.xaxis()
  
  # Draw predicted curve
  grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = '#0a71b4', lwd = 2), default = 'native')
  
  # Draw confidence intervals
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = '#0a71b4', lwd = 2, lty = 2), default = 'native')
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = '#0a71b4', lwd = 2, lty = 2), default = 'native')
  
  # Add text labels
  grid.text(y = unit(0, 'npc') - unit(2.5, 'lines'), label = 'Mean Relative Abundance (log10)', gp = gpar(fontface = 2))
  grid.text(x = unit(0, 'npc') - unit(3, 'lines'), label = 'Frequency of Occurrence', gp = gpar(fontface = 2), rot = 90)
  
  # Add compartment name and statistical results
  grid.text(paste(compartment_name, "\nRsqr=", round(Rsqr, 3), "\nm=", round(coef(m.fit), 4), "\nNm=", round(coef(m.fit) * N, 2)), 
            x = unit(0.8, "npc"), y = unit(0.8, "npc"), just = c("centre", "centre"), gp = gpar(cex = 0.8, fontface = "bold"))
  
  # Capture current plot
  current_plot <- grid.grab()
  
  # Store results
  results <- list(
    compartment = compartment_name,
    m = coef(m.fit),
    Nm = coef(m.fit) * N,
    Rsqr = Rsqr,
    m_ci = m.ci
  )
  
  return(list(plot = current_plot, results = results))
}

# Analyze five compartments separately
compartments <- list(
  N0 = N0,
  N2 = N2,
  N4 = N4,
  N6 = N6
)

compartment_names <- c("N0", "N2", "N4", "N6")

# Run analysis and collect results
for (i in 1:length(compartment_names)) {
  name <- compartment_names[i]
  data <- compartments[[name]]
  analysis_result <- run_ncm_analysis(data, name)
  plot_list[[name]] <- analysis_result$plot
  results_list[[name]] <- analysis_result$results
}

# Create combined plot - five plots arranged vertically
grid.newpage()
grid.arrange(
  plot_list$N0,
  plot_list$N2,
  plot_list$N4,
  plot_list$N6,
  ncol = 1,
  heights = c(1, 1, 1, 1)
)

# Save combined plot
ggsave("NCM_Rhizobia_Nitrogen.pdf", 
       arrangeGrob(plot_list$N0, plot_list$N2, plot_list$N4, plot_list$N6, ncol = 1),
       width = 5, height = 15)
