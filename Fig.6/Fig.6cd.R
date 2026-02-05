rm(list = ls())
library(piecewiseSEM) 
library(readxl)
library(nlme)
library(lme4)
library(QuantPsyc)
library(reshape2)
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.6/")

# Load main data
df <- read.csv("SEM input data.csv", header = TRUE, row.names = 1)

# Nutrient accumulation index calculation
EFC <- read.csv("SEM_difference in soil factors between RZS and MFS.csv", row.names = 1)
mydata <- EFC
mydata <- mydata[c("TN", "NH4_N", "NO3_N", "MBN", "TC", "OM", "MBC", "TP", "AP", "TK", "AK")]
mydata <- scale(mydata)
head(mydata)
Zcore <- (mydata - mean(mydata)) / sd(mydata)
Zcore <- as.data.frame(Zcore)
Zcore$Nutrition <- rowSums(Zcore) / 11
Zcore <- Zcore[rownames(df), ]
df$emf_deltaNutrition <- Zcore$Nutrition

# Nitrification function index calculation
mydata <- df[13:24]
mydata <- mydata[c("hzsB", "amoA_archaea", "amoA_becteria", "nxrA", "comammoxA")]
mydata <- scale(mydata)
head(mydata)
Zcore <- (mydata - mean(mydata)) / sd(mydata)
Zcore <- as.data.frame(Zcore)
Zcore$emf_Nitrification <- rowSums(Zcore) / 5
Zcore <- Zcore[rownames(df), ]
df$emf_Nitrification <- Zcore$emf_Nitrification

# Denitrification function index calculation
mydata <- df[13:24]
mydata <- mydata[c("napA", "nirS", "nirK", "qnorB", "nosZ")]
mydata <- scale(mydata)
head(mydata)
Zcore <- (mydata - mean(mydata)) / sd(mydata)
Zcore <- as.data.frame(Zcore)
Zcore$emf_Denitrification <- rowSums(Zcore) / 5
Zcore <- Zcore[rownames(df), ]
df$emf_Denitrification <- Zcore$emf_Denitrification

# Nitrogen fixation function index calculation
mydata <- df[13:24]
mydata <- mydata[c("nifH")]
mydata <- scale(mydata)
head(mydata)
Zcore <- (mydata - mean(mydata)) / sd(mydata)
Zcore <- as.data.frame(Zcore)
Zcore$emf_Nfix <- rowSums(Zcore) / 1
Zcore <- Zcore[rownames(df), ]
df$emf_Nfix <- Zcore$emf_Nfix

# Build piecewise SEM model
model <- psem(
  lm(emf_deltaNutrition ~ richness_AMF + richness_RHI + Nitrogen + emf_Nitrification + emf_Denitrification + emf_Nfix, df),
  lm(emf_Nitrification ~ richness_AMF + richness_RHI + Nitrogen + NicheSelc, df), 
  lm(emf_Denitrification ~ richness_AMF + richness_RHI + Nitrogen + NicheSelc, df), 
  lm(emf_Nfix ~ richness_AMF + richness_RHI + Nitrogen + NicheSelc, df), 
  lm(richness_AMF ~ Nitrogen + NicheSelc, df),
  lm(richness_RHI ~ Nitrogen + NicheSelc, df)
)

# Plot SEM diagram
plot(model, layout = "tree")

# Model summary
summary(model, .progressBar = FALSE)

# Save summary to text file
sink("SEM_summary.txt")
summary(model, .progressBar = FALSE)
sink()

# Model evaluation using lavaan
library(lavaan)

# Define lavaan model syntax
model_lavaan <- '
  # Regression equations
  emf_deltaNutrition ~ richness_AMF + richness_RHI + Nitrogen + emf_Nitrification + emf_Denitrification + emf_Nfix
  emf_Nitrification ~ richness_AMF + richness_RHI + Nitrogen + NicheSelc
  emf_Denitrification ~ richness_AMF + richness_RHI + Nitrogen + NicheSelc
  emf_Nfix ~ richness_AMF + richness_RHI + Nitrogen + NicheSelc
  richness_AMF ~ Nitrogen + NicheSelc
  richness_RHI ~ Nitrogen + NicheSelc
'

# Fit lavaan model
fit_lavaan <- sem(model_lavaan, data = df)

# Standardize data to address scale differences
df_scaled <- as.data.frame(scale(df[c(1:3, 7:49)]))

# Fit lavaan model with standardized data
fit_lavaan <- sem(model_lavaan, data = df_scaled)

# Output model fit indices
summary(fit_lavaan, fit.measures = TRUE)

# Save fit measures to text file
sink("SEM_fit_measures.txt")
summary(fit_lavaan, fit.measures = TRUE)
sink()

# Calculate direct, indirect, and total effects
coefs_df <- coefs(model)
print("=== Raw coefs_df ===")
print(coefs_df)

# Check if there are valid unidirectional regression paths after filtering
if (nrow(coefs_df) == 0) {  
  stop("No valid unidirectional regression paths after filtering. Please check model or data.")
}

coefs_df$effect <- coefs_df$Std.Estimate

# Function to build graph from coefficient data
build_graph <- function(coefs_df) {  
  graph <- list()  
  for (i in 1:nrow(coefs_df)) {    
    from <- coefs_df$Predictor[i]    
    to <- coefs_df$Response[i]    
    effect <- coefs_df$effect[i]    
    if (!is.null(graph[[from]])) {      
      graph[[from]] <- rbind(graph[[from]], data.frame(to = to, effect = effect, stringsAsFactors = FALSE))    
    } else {      
      graph[[from]] <- data.frame(to = to, effect = effect, stringsAsFactors = FALSE)    
    }  
  }  
  return(graph)
}

graph <- build_graph(coefs_df)
print("=== Constructed graph ===")
print(graph)

# Function to find all paths between nodes
find_all_paths <- function(graph, start, end, visited = character()) {  
  if (start == end) {    
    return(list(c(end)))  
  }  
  if (!start %in% names(graph)) {    
    return(list())  
  }  
  visited <- c(visited, start)  
  paths <- list()  
  for (i in 1:nrow(graph[[start]])) {    
    next_node <- graph[[start]]$to[i]    
    if (!(next_node %in% visited)) {      
      sub_paths <- find_all_paths(graph, next_node, end, visited)      
      for (sp in sub_paths) {        
        paths <- c(paths, list(c(start, sp)))      
      }    
    }  
  }  
  return(paths)
}

# Function to calculate path effect
path_effect <- function(path, graph) {  
  eff <- 1  
  for (i in 1:(length(path) - 1)) {    
    from <- path[i]    
    to <- path[i + 1]    
    edge <- graph[[from]]    
    eff <- eff * edge$effect[edge$to == to]  
  }  
  return(eff)
}

# Function to compute all effects (direct, indirect, total)
compute_all_effects <- function(graph, variables) {  
  results <- data.frame(    
    predictor = character(),    
    outcome = character(),    
    direct_effect = numeric(),    
    indirect_effect = numeric(),    
    total_effect = numeric(),    
    stringsAsFactors = FALSE  
  )
  
  for (pred in variables) {    
    for (outc in variables) {      
      if (pred != outc) {        
        paths <- find_all_paths(graph, pred, outc)
        if (length(paths) > 0) {          
          direct_effect <- 0          
          indirect_effect <- 0          
          for (p in paths) {            
            eff <- path_effect(p, graph)            
            if (length(p) == 2) {              
              direct_effect <- direct_effect + eff            
            } else {              
              indirect_effect <- indirect_effect + eff            
            }          
          }          
          total_effect <- direct_effect + indirect_effect          
          results <- rbind(results, data.frame(            
            predictor = pred,            
            outcome = outc,            
            direct_effect = direct_effect,            
            indirect_effect = indirect_effect,            
            total_effect = total_effect,            
            stringsAsFactors = FALSE          
          ))        
        }     
      }    
    } 
  }  
  return(results)
}

all_vars <- unique(c(coefs_df$Response, coefs_df$Predictor))
effects_df <- compute_all_effects(graph, all_vars)
print("=== Calculated effect results effects_df ===")
print(effects_df)

if (nrow(effects_df) == 0) {  
  stop("No effect values calculated. Please check if paths exist.")
}

# Custom order for displaying three types of effects
metrics_to_show <- c("direct_effect", "indirect_effect", "total_effect")
plot_df <- melt(effects_df,                 
                id.vars = c("predictor", "outcome"),                 
                measure.vars = metrics_to_show,                
                variable.name = "effect_type",                 
                value.name = "effect_size")

# Custom order for X-axis display and select metrics to show
plot_df$effect_type <- factor(plot_df$effect_type, levels = metrics_to_show)
desired_order <- c("Nitrogen", "NicheSelc", "richness_AMF", "richness_RHI", "emf_Nitrification", "emf_Denitrification", "emf_Nfix")
plot_df <- plot_df %>% filter(predictor %in% desired_order)
plot_df$predictor <- factor(plot_df$predictor, levels = desired_order)

# Filter for specific outcome variable (adjust based on final variable names)
plot_df_filtered <- plot_df %>% filter(outcome == "emf_deltaNutrition")
plot_df_filtered

# Plotting
library(ggplot2)
plot_effects <- ggplot(plot_df_filtered, aes(x = predictor, y = effect_size, fill = effect_type)) +  
  geom_bar(stat = "identity", position = "dodge") +   
  labs(x = "Predictor", y = "Effect Size (Standardized)",        
       title = "Effects on TRAD: Direct, Indirect, and Total") +  
  scale_fill_manual(values = c("direct_effect" = "#004EA2",                                
                               "indirect_effect" = "#009D8A",                                
                               "total_effect" = "#F5B71A")) +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),        
        legend.title = element_blank()) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +  # Add horizontal line  
  ylim(-0.5, 0.5)  # Set Y-axis range

plot_effects

# Save plot
ggsave(plot = plot_effects, filename = "SEM effect size.pdf", device = "pdf", width = 10, height = 8, dpi = 300)

# Save effect size data
write.csv(plot_df_filtered, "SEM effect size.csv")