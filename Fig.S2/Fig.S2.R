# =============================================================================
# File: 02_Phylum_level_analysis.R
# Description: Bacterial phylum (Phylum) level relative abundance analysis and stacked bar chart plotting
# Strictly follows the original code logic, only modularized and path adjusted
# =============================================================================

# Environment setup -------------------------------------------------------------------
# Clear workspace
rm(list = ls())

# Load necessary packages (exactly the same as original code)
library(dplyr)
library(reshape2)
library(ggplot2)

# Set random seed for reproducible results
set.seed(123)

# Set working directory
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S2/")

# Load data -------------------------------------------------------------------
# Read OTU table (exactly the same as original code, only path modified)
otu <- read.csv("Bacteria_ASV_rarefaction.csv", header = TRUE)

# Read taxonomy annotation file
tax <- read.csv("Bacteria_tax.csv", header = TRUE)

# Read group information
group <- read.csv("group.csv", header = TRUE)

# Data preprocessing: Convert to relative abundance (exactly the same as original code)-----------------------
df <- otu
df_rel <- df
df_rel[, 2:ncol(df)] <- 100 * df[, 2:ncol(df)] / colSums(df[, 2:ncol(df)])
otu <- df_rel

Order <- merge(otu, tax, by = "ASV")
Order <- Order[, grepl("N0|N2|N4|N6|Order", colnames(Order))]

sum_Order <- group_by(Order, Order) %>%
  summarise_all(sum)
sum_Order$Order[is.na(sum_Order$Order) | sum_Order$Order == ""] <- "Others"

### Sort by total abundance size, set the order of each order
sum_Order$row_sum <- rowSums(sum_Order[, -1], na.rm = TRUE)  # Calculate sum for each row
sorted_data <- sum_Order[order(-sum_Order$row_sum), ]  # Sort dataframe by row sum size
sorted_data$row_sum <- NULL

#### Calculate top twenty at order level table
Ord_20 <- sorted_data[c(1:21), ]
Ord_20 <- subset(Ord_20, Ord_20$Order != "Others")
name_20 <- Ord_20$Order

Ord_other <- subset(sorted_data, !grepl(paste(name_20, collapse = "|"), Order))
Ord_other <- colSums(Ord_other[2:ncol(Ord_other)])
Ord_other <- c("Others", Ord_other)

### Merge top 20 rows and the 21st row
Ord_20_other <- rbind(Ord_20, Ord_other)
Ord_21_name <- Ord_20_other$Order

### Wide to long format, add grouping information
Ord_long <- melt(Ord_20_other, id.vars = "Order", variable.name = "sample", value.name = "Value")
Ord_long_group <- merge(Ord_long, group, by = "sample")

### Factorize, set display order
Ord_long_group$treatment <- factor(Ord_long_group$treatment, levels = c("N0", "N2", "N4", "N6"))
Ord_long_group$location <- factor(Ord_long_group$Niche,
                                  levels = c("Original", "Rootzone", "Rhizosphere", "Root", "Nodule"),
                                  labels = c("MFS", "RZS", "RS", "RE", "NE"))  # Use abbreviations

Ord_long_group$Order <- factor(Ord_long_group$Order, levels = Ord_21_name)

### Set colors for all orders (exactly the same as original code)
CS_cols <- c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3',
             '#fdb462', '#b3de69', '#fccde5', '#ccebc5', '#bcf4f5',
             '#ffb6c1', '#00ffff', '#9400d3', '#228B22', '#FF6347',
             '#ffed6f', '#00FF7F', '#EEE8AA', '#4682B4', '#FFA500',
             'gray')

names(CS_cols) <- c(Ord_21_name)
Ord_long_group$Value <- as.numeric(Ord_long_group$Value)

# Plotting (exact same aesthetic mappings as original code)--------------------------------------------------
p <- ggplot(Ord_long_group, aes(treatment, Value/5, fill = Order)) +
  geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~location, scales = 'free_x', ncol = 5) +
  scale_fill_manual(values = CS_cols) +
  labs(x = NULL, y = 'Bacteria Relative Abundance(%)') +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        strip.text = element_text(family = "serif", size = 12), 
        axis.text = element_text(family = "serif", size = 11), 
        axis.title = element_text(family = "serif", size = 13), 
        legend.text = element_text(family = "serif", size = 11, face = "italic"),  # Order level in italics
        legend.position = "bottom",
        legend.key.size = unit(0.7, "lines"))

# Display plot
print(p)

# Save plot (exactly same parameters as original code)
ggsave("16S_Order_relative_abundance.pdf", plot = p, device = "pdf", width = 7.5, height = 4)
#######################AMF Genus####################################################
# =============================================================================
# File: 03_AMF_Genus_level_analysis.R
# Description: AMF fungal genus (Genus) level relative abundance analysis and stacked bar chart plotting
# Strictly follows the original code logic, only data source and classification level modified
# =============================================================================
rm(list = ls())

# Load necessary packages (exactly the same as original code)
library(dplyr)
library(reshape2)
library(ggplot2)

# Set random seed for reproducible results
set.seed(123)

# Set working directory
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S2/")

# Load data -------------------------------------------------------------------
# Read OTU table (exactly the same as original code, only path modified)
otu <- read.csv("AMF_ASV_rarefaction.csv", header = TRUE)

# Read taxonomy annotation file
tax <- read.csv("AMF_tax.csv", header = TRUE)

# Read group information
group <- read.csv("group.csv", header = TRUE)

# Data preprocessing: Convert to relative abundance (exactly the same as original code)-----------------------
df <- otu
df_rel <- df
df_rel[, 2:ncol(df)] <- 100 * df[, 2:ncol(df)] / colSums(df[, 2:ncol(df)])
otu <- df_rel

# Check column name for genus level in tax file
genus_col <- "Genus"  # Assuming column name is "Genus", modify if not

# Check if Genus column exists in taxonomy file
if (!genus_col %in% colnames(tax)) {
  # Try to find other possible names for Genus column
  possible_names <- grep("genus|Genus", colnames(tax), ignore.case = TRUE, value = TRUE)
  if (length(possible_names) > 0) {
    genus_col <- possible_names[1]
    message(paste("Using column name:", genus_col, "as Genus classification column"))
  } else {
    stop("Genus classification column not found, please check AMF taxonomy file")
  }
}

# Merge OTU table and taxonomy annotation
Genus <- merge(otu, tax, by = "ASV")

# Extract sample columns and Genus column
# Assuming sample columns contain "N0|N2|N4|N6"
Genus <- Genus[, grepl("N0|N2|N4|N6|Genus", colnames(Genus))]

# If column name is not "Genus", rename to "Genus"
colnames(Genus)[colnames(Genus) == genus_col] <- "Genus"

# Group by Genus and sum
sum_Genus <- group_by(Genus, Genus) %>%
  summarise_all(sum)

# Handle blank classifications
blank_count <- sum(is.na(sum_Genus$Genus) | 
                     trimws(as.character(sum_Genus$Genus)) == "" |
                     sum_Genus$Genus %in% c(" ", "NA", "N/A", "unknown", "unclassified"))

if (blank_count > 0) {
  message(paste("Found", blank_count, "blank classifications, filled as 'Others'"))
  sum_Genus$Genus[is.na(sum_Genus$Genus) | 
                    trimws(as.character(sum_Genus$Genus)) == "" |
                    sum_Genus$Genus %in% c(" ", "NA", "N/A", "unknown", "unclassified")] <- "Others"
} else {
  message("No blank classifications found")
}

### Sort by total abundance size, set the order of each genus
sum_Genus$row_sum <- rowSums(sum_Genus[, -1], na.rm = TRUE)  # Calculate sum for each row
sorted_data <- sum_Genus[order(-sum_Genus$row_sum), ]  # Sort dataframe by row sum size
sorted_data$row_sum <- NULL

#### Calculate top twenty at genus level table
Genus_20 <- sorted_data[c(1:21), ]
Genus_20 <- subset(Genus_20, Genus_20$Genus != "Others")
name_20 <- Genus_20$Genus

Genus_other <- subset(sorted_data, !grepl(paste(name_20, collapse = "|"), Genus))
Genus_other <- colSums(Genus_other[2:ncol(Genus_other)])
Genus_other <- c("Others", Genus_other)

### Merge top 20 rows and the 21st row
Genus_20_other <- rbind(Genus_20, Genus_other)
Genus_21_name <- Genus_20_other$Genus

### Wide to long format, add grouping information
Genus_long <- melt(Genus_20_other, id.vars = "Genus", variable.name = "sample", value.name = "Value")
Genus_long_group <- merge(Genus_long, group, by = "sample")

### Factorize, set display order
# Note: AMF group file may have different column name structure, adjust according to actual situation
# Assuming AMF grouping information is same as 16S
Genus_long_group$treatment <- factor(Genus_long_group$treatment, levels = c("N0", "N2", "N4", "N6"))
Genus_long_group$location <- factor(Genus_long_group$Niche,
                                    levels = c("Original", "Rootzone", "Rhizosphere", "Root", "Nodule"),
                                    labels = c("MFS", "RZS", "RS", "RE", "NE"))  # Use abbreviations

Genus_long_group$Genus <- factor(Genus_long_group$Genus, levels = Genus_21_name)

### Set colors for all genera (exactly the same as original code)
CS_cols <- c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3',
             '#fdb462', '#b3de69', '#fccde5', 'gray')

names(CS_cols) <- c(Genus_21_name)
Genus_long_group$Value <- as.numeric(Genus_long_group$Value)


p <- ggplot(Genus_long_group, aes(treatment, Value/5, fill = Genus)) +
  geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~location, scales = 'free_x', ncol = 5) +
  scale_fill_manual(values = CS_cols) +
  labs(x = NULL, y = 'AMF Relative Abundance(%)') +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        strip.text = element_text(family = "serif", size = 12), 
        axis.text = element_text(family = "serif", size = 11), 
        axis.title = element_text(family = "serif", size = 13), 
        legend.text = element_text(family = "serif", size = 11, face = "italic"),  # Genus level in italics
        legend.position = "bottom",
        legend.key.size = unit(0.7, "lines"))
print(p)

# Save plot
ggsave("AMF_Genus_relative_abundance.pdf", plot = p, device = "pdf", width = 7.5, height = 4)
########################Rhizobia Genus###################################
# =============================================================================
# File: 04_Rhizobia_Genus_level_analysis.R
# Description: Rhizobia genus level relative abundance analysis and stacked bar chart plotting
# =============================================================================

rm(list = ls())
library(dplyr)
library(reshape2)
library(ggplot2)
# Set random seed for reproducible results
set.seed(123)

# Set working directory
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S2/")

# Load data -------------------------------------------------------------------
# Read OTU table (exactly the same as original code, only path modified)
otu <- read.csv("Rhizobia_ASV_rarefaction.csv", header = TRUE)

# Read taxonomy annotation file
tax <- read.csv("Rhizobia_tax.csv", header = TRUE)

# Read group information
group <- read.csv("group.csv", header = TRUE)


# Data preprocessing: Convert to relative abundance (exactly the same as original code)-----------------------
df <- otu
df_rel <- df
df_rel[, 2:ncol(df)] <- 100 * df[, 2:ncol(df)] / colSums(df[, 2:ncol(df)])
otu <- df_rel


# Check column name for genus level in tax file
genus_col <- "Genus"  # Assuming column name is "Genus", modify if not

# Check if Genus column exists in taxonomy file
if (!genus_col %in% colnames(tax)) {
  # Try to find other possible names for Genus column
  possible_names <- grep("genus|Genus", colnames(tax), ignore.case = TRUE, value = TRUE)
  if (length(possible_names) > 0) {
    genus_col <- possible_names[1]
    message(paste("Using column name:", genus_col, "as Genus classification column"))
  } else {
    # If still not found, Rhizobia annotation file may have different column names, try other common names
    possible_names <- grep("属|genus|Genus", colnames(tax), value = TRUE)
    if (length(possible_names) > 0) {
      genus_col <- possible_names[1]
      message(paste("Using column name:", genus_col, "as Genus classification column"))
    } else {
      stop("Genus classification column not found, please check Rhizobia taxonomy file")
    }
  }
}

# Merge OTU table and taxonomy annotation
Genus <- merge(otu, tax, by = "ASV")

# Extract sample columns and Genus column
# Assuming sample columns contain "N0|N2|N4|N6"
Genus <- Genus[, grepl("N0|N2|N4|N6|Genus", colnames(Genus))]

# If column name is not "Genus", rename to "Genus"
colnames(Genus)[colnames(Genus) == genus_col] <- "Genus"

# Group by Genus and sum
sum_Genus <- group_by(Genus, Genus) %>%
  summarise_all(sum)

# Handle blank classifications
blank_count <- sum(is.na(sum_Genus$Genus) | 
                     trimws(as.character(sum_Genus$Genus)) == "" |
                     sum_Genus$Genus %in% c(" ", "NA", "N/A", "unknown", "unclassified", "Unclassified"))

if (blank_count > 0) {
  message(paste("Found", blank_count, "blank classifications, filled as 'Others'"))
  sum_Genus$Genus[is.na(sum_Genus$Genus) | 
                    trimws(as.character(sum_Genus$Genus)) == "" |
                    sum_Genus$Genus %in% c(" ", "NA", "N/A", "unknown", "unclassified", "Unclassified")] <- "Others"
} else {
  message("No blank classifications found")
}

### Sort by total abundance size, set the order of each genus
sum_Genus$row_sum <- rowSums(sum_Genus[, -1], na.rm = TRUE)  # Calculate sum for each row
sorted_data <- sum_Genus[order(-sum_Genus$row_sum), ]  # Sort dataframe by row sum size
sorted_data$row_sum <- NULL

#### Calculate top twenty at genus level table
Genus_20 <- sorted_data[c(1:21), ]
Genus_20 <- subset(Genus_20, Genus_20$Genus != "Others")
name_20 <- Genus_20$Genus

Genus_other <- subset(sorted_data, !grepl(paste(name_20, collapse = "|"), Genus))
Genus_other <- colSums(Genus_other[2:ncol(Genus_other)])
Genus_other <- c("Others", Genus_other)

### Merge top 20 rows and the 21st row
Genus_20_other <- rbind(Genus_20, Genus_other)
Genus_21_name <- Genus_20_other$Genus

### Wide to long format, add grouping information
Genus_long <- melt(Genus_20_other, id.vars = "Genus", variable.name = "sample", value.name = "Value")
Genus_long_group <- merge(Genus_long, group, by = "sample")

### Factorize, set display order
# Note: Rhizobia group file may have different column name structure, adjust according to actual situation
# Assuming Rhizobia grouping information is same as 16S, using treatment and Niche columns
# If column names are different, modify accordingly

# Check if treatment column exists in group file, if not, use other column names
if (!"treatment" %in% colnames(group)) {
  # Try to find other possible names for treatment column
  possible_treatment <- grep("treatment|Treatment|fert|Fert", colnames(group), value = TRUE)
  if (length(possible_treatment) > 0) {
    treatment_col <- possible_treatment[1]
    # Temporarily create treatment column for unified subsequent code
    group$treatment <- group[[treatment_col]]
    message(paste("Using column name", treatment_col, "as treatment column"))
  } else {
    stop("Treatment column not found, please check Rhizobia group file")
  }
}

# Check if Niche column exists in group file, if not, use other column names
if (!"Niche" %in% colnames(group)) {
  # Try to find other possible names for Niche column
  possible_niche <- grep("niche|Niche|location|Location|sample_type", colnames(group), value = TRUE)
  if (length(possible_niche) > 0) {
    niche_col <- possible_niche[1]
    # Temporarily create Niche column for unified subsequent code
    group$Niche <- group[[niche_col]]
    message(paste("Using column name", niche_col, "as Niche column"))
  } else {
    stop("Niche column not found, please check Rhizobia group file")
  }
}

Genus_long_group$treatment <- factor(Genus_long_group$treatment, levels = c("N0", "N2", "N4", "N6"))
Genus_long_group$location <- factor(Genus_long_group$Niche,
                                    levels = c("Original", "Rootzone", "Rhizosphere", "Root", "Nodule"),
                                    labels = c("MFS", "RZS", "RS", "RE", "NE"))  # Use abbreviations

Genus_long_group$Genus <- factor(Genus_long_group$Genus, levels = Genus_21_name)

### Set colors for all genera (exactly the same as original code)
CS_cols <- c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3',
             '#fdb462', '#b3de69', '#fccde5', '#ccebc5', '#bcf4f5',
             '#ffb6c1', '#00ffff', '#9400d3', '#228B22', '#FF6347',
             '#ffed6f', '#00FF7F', '#EEE8AA', '#4682B4', '#FFA500',
             'gray')

names(CS_cols) <- c(Genus_21_name)
Genus_long_group$Value <- as.numeric(Genus_long_group$Value)


p <- ggplot(Genus_long_group, aes(treatment, Value/5, fill = Genus)) +
  geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~location, scales = 'free_x', ncol = 5) +
  scale_fill_manual(values = CS_cols) +
  labs(x = NULL, y = 'Rhizobia Relative Abundance(%)') +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        strip.text = element_text(family = "serif", size = 12), 
        axis.text = element_text(family = "serif", size = 11), 
        axis.title = element_text(family = "serif", size = 13), 
        legend.text = element_text(family = "serif", size = 11, face = "italic"),  # Genus level in italics
        legend.position = "bottom",
        legend.key.size = unit(0.7, "lines"))

# Display plot
print(p)

# Save plot (exactly same parameters as original code)
ggsave("Rhizobia_Genus_relative_abundance.pdf", plot = p, device = "pdf", width = 7.5, height = 4)